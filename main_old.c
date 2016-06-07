#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memset */
#include <sys/time.h>
#include <mpi.h>
#include <omp.h>
#include <unistd.h> /* getopt? */
#include <ctype.h>  /* getopt? */
#include <math.h>

#define END_READ_RESERVE 1024 // Bytes
#define DEBUG 1
// #define RAW_FINAL_OUTPUT 1
 
#ifdef DEBUG
    FILE *fdebug;
#endif

char *in_path = "/Users/larry/Desktop/tests/";
char *in_file_name = "IT_Met_1a.chr22.fa";
char *out_path = "/Users/larry/Desktop/outs/IT_Met_1a.chr22/";

 // HELP FCNS
 size_t getFileSize(FILE *f)
{
    fseek(f, 0, SEEK_END);
    size_t fs = ftell(f);
    fseek(f, 0, SEEK_SET);
    return fs;
}

size_t getEndOfSequenceForward(char *buf, size_t length) 
{
    size_t i;
    size_t end_offset = 0;
    for (i = 0; i < length; i++) {
        end_offset++;
        
        if (buf[i] == '\n') {
            return end_offset;
        }
    }
    return 0;
}

typedef struct KmerCountPair { // [16B]
    uint64_t kmer;
    uint64_t count;
} KmerCountPair;
 
// COUNT PARTS
typedef struct CountPart {
    int id;
    size_t start_nmer;
    size_t end_nmer;
    size_t kmers_count;
    int files_count;
    size_t *start_offsets; // v kazdom subore kde je zaciatok (kolko k-tic preskocit)
    size_t *end_offsets; // v kazdom subore kde je zaciatok (kolko k-tic preskocit - kde konci)
    size_t *files_size;
} CountPart;
 
// CORE PARTS
typedef struct CorePart {
    int id; // inner id
    size_t start_index;
    size_t end_index;
    size_t size;
    // char file_full_path[1024];
} CorePart;
 
 // THREAD PARTS
 
 typedef struct ThreadPart {
	int id;
	size_t start_index; // local for buffer (not for file)
	size_t end_index;   // local for buffer (not for file)
	size_t size;
} ThreadPart;

void getThreadPartsStartsAndEnds(char *buffer, CorePart corePart, size_t reserve, int nthreads, ThreadPart *threadParts) {
    size_t part_size = (corePart.size / nthreads) + reserve;

    int i;
    for (i = 0; i < nthreads; i++) {
        size_t start = i * part_size;
        size_t end = start + part_size - reserve;

        size_t move_start_by = 0;
        size_t move_end_by = 0;
        
        char *buf = (char*)malloc(sizeof(char) * reserve);
        memset(buf, '\0', reserve);
        
        if (start == 0) { // zaciatok // na zaciatku neposuvat zaciatocny pointer
            strncpy(buf, &buffer[end], reserve);
            move_end_by = getEndOfSequenceForward(buf, reserve);
        } else if (end < corePart.end_index) { // medzi
            start -= reserve;
            move_start_by = getEndOfSequenceForward(buf, reserve) + 1;

            memset(buf, '\0', reserve);
            
            if (end > corePart.size) { end = corePart.size; } // => lebo segmentation fault
            strncpy(buf, &buffer[end], reserve);
            move_end_by = getEndOfSequenceForward(buf, reserve);
        } else {
            start -= reserve;
            strncpy(buf, &buffer[start], reserve);
            move_start_by = getEndOfSequenceForward(buf, reserve) + 1;
        }
        start += move_start_by;
        end += move_end_by;
        
        if (end > corePart.size) {
            end = corePart.size;
        }
        threadParts[i].id = i;
        threadParts[i].start_index = start;
        threadParts[i].end_index = end;
        threadParts[i].size = end - start;
        free(buf);
    }
}
 
// FILE PARTS
// file is splitted according to number of cores(processors)
typedef struct FilePart 
{
    int pid;
    size_t start_index;
    size_t end_index;
    size_t size;
    CorePart *coreParts;
    int num_of_core_parts;
} FilePart;

// splits input file into parts where reads can be splitted only at their ends (end of lines) => no chopping reads somewhere in the middle
void getFilePartsStartsAndEnds(FILE *f, size_t filesize, size_t reserve, FilePart *fileParts, int file_parts_count, size_t memory_for_kmers, size_t memory_for_raw_part) 
{
    size_t core_part_size = (filesize / file_parts_count) + reserve;
    char *buf = (char*)malloc(sizeof(char) * reserve);
    int i;
    for (i = 0; i < file_parts_count; i++) {
        size_t start = i * core_part_size;
        size_t end = start + core_part_size - reserve;
        size_t move_start_by = 0;
        size_t move_end_by = 0;
        memset(buf, '\0', reserve);
        
        if (start == 0) { // beginning of input file
            fseek(f, end, SEEK_SET);
            fread(buf, sizeof(char), reserve, f);
            move_end_by = getEndOfSequenceForward(buf, reserve);
        } else if (end < filesize) { // betweeen
            start -= reserve;
            fseek(f, start, SEEK_SET);
            fread(buf, sizeof(char), reserve, f);
            move_start_by = getEndOfSequenceForward(buf, reserve);

            memset(buf, '\0', reserve);

            fseek(f, end, SEEK_SET);
            fread(buf, sizeof(char), reserve, f);
            move_end_by = getEndOfSequenceForward(buf, reserve);
        } else { // end
            start -= reserve;
            fseek(f, start, SEEK_SET);
            fread(buf, sizeof(char), reserve, f);
            move_start_by = getEndOfSequenceForward(buf, reserve);
        }
        start += move_start_by;
        end += move_end_by;
        
        if (end > filesize) {
            end = filesize;
        }
        fileParts[i].pid = i;
        fileParts[i].start_index = start;
        fileParts[i].end_index = end;
        fileParts[i].size = end - start;
        
        // kolko parts bude mat filePart = <velkost vstupneho retazca> / max pamat pre string
        int num_of_core_parts = 1 + (end - start) / memory_for_raw_part;

        fileParts[i].coreParts = (CorePart*)malloc(sizeof(CorePart) * num_of_core_parts);
        fileParts[i].num_of_core_parts = num_of_core_parts;
        
        size_t core_part_size = ((end - start) / num_of_core_parts) + reserve;
        int j;
        char *p_buf = (char*)malloc(sizeof(char) * reserve);
        for (j = 0; j < num_of_core_parts; j++) {
            size_t p_start = fileParts[i].start_index + j * core_part_size;
            size_t p_end = p_start + core_part_size - reserve;
            
            size_t p_move_start_by = 0;
            size_t p_move_end_by = 0;
            memset(p_buf, '\0', reserve);
            
            if (p_start == fileParts[i].start_index) {
                fseek(f, p_end, SEEK_SET);
                fread(p_buf, sizeof(char), reserve, f);
                p_move_end_by = getEndOfSequenceForward(p_buf, reserve);
            } else if (p_end < fileParts[i].end_index) {
                p_start -= reserve;
                fseek(f, p_start, SEEK_SET);
                fread(p_buf, sizeof(char), reserve, f);
                p_move_start_by = getEndOfSequenceForward(p_buf, reserve);

                memset(p_buf, '\0', reserve);
                
                fseek(f, p_end, SEEK_SET);
                fread(p_buf, sizeof(char), reserve, f);
                p_move_end_by = getEndOfSequenceForward(p_buf, reserve);
            } else {
                p_start -= reserve;
                fseek(f, p_start, SEEK_SET);
                fread(p_buf, sizeof(char), reserve, f);
                p_move_start_by = getEndOfSequenceForward(p_buf, reserve);
            }
            p_start += p_move_start_by;
            p_end += p_move_end_by;
            
            if (p_end > fileParts[i].end_index) {
                p_end = fileParts[i].end_index;
            }
            fileParts[i].coreParts[j].id = j;
            fileParts[i].coreParts[j].start_index = p_start;
            fileParts[i].coreParts[j].end_index = p_end;
            fileParts[i].coreParts[j].size = p_end - p_start;
        }
        free(p_buf);
        
// #ifdef DEBUG
//         fprintf(fdebug, "\ncore_part_size = %d\n", core_part_size);
//         fprintf(fdebug, "memory_for_raw_part = %lu\n", memory_for_raw_part);
// #endif
        
        
    }
    free(buf);
    
    fseek(f, 0, SEEK_SET);
}

void printFilePartData(FilePart filePart) 
{
    printf("pid = %d\n", filePart.pid);
    printf("start_index = %lu\n", filePart.start_index);
    printf("end_index = %lu\n", filePart.end_index);
    printf("size = %lu (req. memory: %lu)\n", filePart.size, filePart.size * sizeof(char));
    printf("%d coreParts:\n", filePart.num_of_core_parts);
    int i;
    for (i = 0; i < filePart.num_of_core_parts; i++) {
        printf("\tpart id = %d\n", filePart.coreParts[i].id);
        printf("\tstart_index = %lu\n", filePart.coreParts[i].start_index);
        printf("\tend_index = %lu\n", filePart.coreParts[i].end_index);
        printf("\tsize = %d (req. memory: %lu)\n", filePart.coreParts[i].size, filePart.coreParts[i].size * sizeof(char));
    }
    printf("-----\n");
}

void printCorePart(CorePart corePart) {
    printf("COREPART %d: start_index = %lu\tend_index = %lu\tsize = %lu\n", corePart.id, corePart.start_index, corePart.end_index, corePart.size);
}

void printCountPart(CountPart countPart) {
    printf("CountPart ID = %d\tstart_nmer = %lu\tend_nmer = %lu\tkmers_count = %lu\n", countPart.id, countPart.start_nmer, countPart.end_nmer, countPart.kmers_count);
    int i;
    for (i = 0; i < countPart.files_count; i++) {
        printf("\tfile #%d: %lu => %lu / %lu (bytes: %lu)\n", i, countPart.start_offsets[i], countPart.end_offsets[i], countPart.files_size[i] / 8, countPart.files_size[i]);
    }
}

void printThreadPart(ThreadPart threadPart) {
    printf("THREAD %d: start_index = %lu\t end_index = %lu\t size = %lu\n", threadPart.id, threadPart.start_index, threadPart.end_index, threadPart.size);
}

#ifdef DEBUG
void printFilePartDataToFile(FilePart filePart, FILE *f) {
    fprintf(f, "pid = %d\n", filePart.pid);
    fprintf(f, "start_index = %lu\n", filePart.start_index);
    fprintf(f, "end_index = %lu\n", filePart.end_index);
    fprintf(f, "size = %lu (req. memory: %lu)\n", filePart.size, filePart.size * sizeof(char));
    fprintf(f, "%d coreParts:\n", filePart.num_of_core_parts);
    int i;
    for (i = 0; i < filePart.num_of_core_parts; i++) {
        fprintf(f, "\tpart id = %d\n", filePart.coreParts[i].id);
        fprintf(f, "\tstart_index = %lu\n", filePart.coreParts[i].start_index);
        fprintf(f, "\tend_index = %lu\n", filePart.coreParts[i].end_index);
        fprintf(f, "\tsize = %d (req. memory: %lu)\n", filePart.coreParts[i].size, filePart.coreParts[i].size * sizeof(char));
    }
    fprintf(f, "-----\n\n");
    fflush(f);
} 

void printCorePartToFile(CorePart corePart, FILE *f) {
    fprintf(f, "COREPART %d: start_index = %lu\tend_index = %lu\tsize = %lu\n", corePart.id, corePart.start_index, corePart.end_index, corePart.size);
    fflush(f);
}

void printCountPartToFile(CountPart countPart, FILE *f) {
    fprintf(f, "CountPart ID = %d\tstart_nmer = %lu\tend_nmer = %lu\tkmers_count = %lu\n", countPart.id, countPart.start_nmer, countPart.end_nmer, countPart.kmers_count);
    int i;
    for (i = 0; i < countPart.files_count; i++) {
        fprintf(f, "\tfile #%d: %lu => %lu / %lu (bytes: %lu)\n", i, countPart.start_offsets[i], countPart.end_offsets[i], countPart.files_size[i] / 8, countPart.files_size[i]);
    }
    fflush(f);
}

void printThreadPartToFile(ThreadPart threadPart, FILE *f) {
    fprintf(f, "THREAD %d: start_index = %lu\t end_index = %lu\t size = %lu\n", threadPart.id, threadPart.start_index, threadPart.end_index, threadPart.size);
    fflush(f);
}

#endif

// KMERS
void createKmers(char *buf, size_t buf_size, int KMER_LENGTH, size_t kmers_limit, uint64_t *kmers, size_t *kmers_count, int *nmers) {
	// printf("buf_size = %lu\tkmers_limit = %d\n", buf_size, kmers_limit);
    char BASES[255];
	BASES[65]=0; BASES[67]=1; BASES[71]=2; BASES[84]=3;

	int i;
 	char BASES_PARSING[255];
 	for (i = 0; i < 85; i++) { BASES_PARSING[i]=0; }
 	BASES_PARSING[65]=1; BASES_PARSING[67]=1; BASES_PARSING[71]=1; BASES_PARSING[84]=1;

	uint64_t BASES_SHIFTED[64][85];
	for (i = 0; i < 64; i += 2) {
		BASES_SHIFTED[i][65] = (uint64_t)0;
		BASES_SHIFTED[i][67] = (uint64_t)1 << i;
		BASES_SHIFTED[i][71] = (uint64_t)2 << i;
		BASES_SHIFTED[i][84] = (uint64_t)3 << i;
	}

    for (i = 0; i < 256; i++) {
        nmers[i] = 0; 
    }
	uint64_t nmers_mask = ((uint64_t)255 << 56);
    
	(*kmers_count) = 0;
	int kmer_length_counter = 0;
	uint64_t kmer = 0;
	for (i = 0; i < buf_size; i++) {
			
			if ((BASES_PARSING[ buf[i] ]) && (*kmers_count) < kmers_limit) {

					if (kmer_length_counter < KMER_LENGTH) {
//						kmer |= (((uint64_t)BASES[buf[i]]) << (62 - kmer_length_counter * 2)); // TODO: ked bude ine ako 31
                        kmer |= BASES_SHIFTED[ 62 - kmer_length_counter * 2 ][buf[i]];
                        kmer_length_counter++;

                        if (kmer_length_counter == KMER_LENGTH) {
                            // printf("kmers_count = %d\n", (*kmers_count));
                            kmers[(*kmers_count)] = kmer;
                            nmers[(int)(kmer >> 56)]++;
                            (*kmers_count)++;
                        }
					} else {
//  					kmer = (uint64_t)((kmer << 2) | ((uint64_t)BASES[buf[i]] << 2)); // TODO: ked bude ine ako 31
                        kmer = (uint64_t)((kmer << 2) | BASES_SHIFTED[2][buf[i]]);
                        // printf("kmers_count = %d\n", (*kmers_count));
                        kmers[(*kmers_count)] = kmer;
                        nmers[(int)(kmer >> 56)]++;
                        (*kmers_count)++;
					}
			} else {
                kmer = 0;
                kmer_length_counter = 0;
			}
	
			if ((*kmers_count) >= kmers_limit) {
// 				(*buf_offset) = i;
				break;
			}
	}
}

void saveKmers(uint64_t *kmers, int kmers_count, int a, int b) {
    char filename_out[1024];
    // sprintf(filename_out, "/Users/larry/Desktop/outs/kmers_%d_%d_%lu_%lu", myid, round * nthreads + tid, tp.start_index, tp.end_index);
    sprintf(filename_out, "%skmers_%d_%d", out_path, a, b);
    FILE *fw = fopen(filename_out, "wb");
    fwrite(kmers, sizeof(uint64_t), kmers_count, fw);
    fclose(fw);
}

void saveKmersRaw(uint64_t *kmers, int kmers_count, int a, int b) {
    int i;
    char filename_kmers[1024];
    sprintf(filename_kmers, "%skmers_raw_%d_%d.txt", out_path, a, b);
    FILE *fw_kmers = fopen(filename_kmers, "w");
    for (i = 0; i < kmers_count; i++) {
        fprintf(fw_kmers, "%llu\n", kmers[i]);
    }                        
    fclose(fw_kmers);
}

void saveNmersCounts(int *nmers_counts, int count, int a, int b) {
    char filename_nmers_bin[1024];
    sprintf(filename_nmers_bin, "%snmers_%d_%d", out_path, a, b);
    FILE *fw_nmers_bin = fopen(filename_nmers_bin, "wb");
    fwrite(nmers_counts, sizeof(int), count, fw_nmers_bin);
    fclose(fw_nmers_bin);
}

void saveNmersCountsRaw(int *nmers_counts, int count, int a, int b) {
    int i;
    // NMERS
    char filename_nmers[1024];
    sprintf(filename_nmers, "%snmers_%d_%d.txt", out_path, a, b);
    FILE *fw_nmers = fopen(filename_nmers, "w");
    for (i = 0; i < count; i++) {
        fprintf(fw_nmers, "%d => %d\n", i, nmers_counts[i]);
    }                        
    fclose(fw_nmers);
}

// SORT
void insertion_sort_uint64(uint64_t *a, int n) {
    int i, j;
    uint64_t value;
    for (i = 1; i < n; i++) {
        value = a[i];
        for (j = i; j > 0 && value < a[j - 1]; j--) {
            a[j] = a[j - 1];
        }
        a[j] = value;
    }
}

void quick_sort_uint64 (uint64_t *data, int *pivInxs, int n) {
    int i, j, k;
    uint64_t pivot, auxn;
    int npivsL = 0;
    int npivsR = 0;

    if (n < 32){
        insertion_sort_uint64(data, n);
        return;
    }

    pivot = data[n / 2];
    i = 0;
    j = n - 1;
    for (; ; i++, j--) {
        while (data[i] < pivot){
            i++;
        }
        while (pivot < data[j]){
            j--;
        }
        if (i >= j)
            break;
       
        if(data[j]==pivot){
            pivInxs[npivsL] = i;
            npivsL++;
        }
        if(data[i]==pivot){
            npivsR++;
            pivInxs[n-npivsR] = j;
        }

        auxn = data[i];
        data[i] = data[j];
        data[j] = auxn;
    }

    int midInx = i-1;
    for(npivsL--;npivsL>=0;npivsL--){
        if(data[midInx]!=pivot){
            auxn = data[midInx];
            data[midInx] = data[pivInxs[npivsL]];
            data[pivInxs[npivsL]] = auxn;
        }
        midInx--;
    }
    int leftEnd = midInx+1;


    midInx = i;
    for(npivsR--;npivsR>=0;npivsR--){
        if(data[midInx]!=pivot){
            auxn = data[midInx];
            data[midInx] = data[pivInxs[n-npivsR-1]];
            data[pivInxs[n-npivsR-1]] = auxn;
        }
        midInx++;
    }
    int rightStart = midInx;

    quick_sort_uint64(data, pivInxs, leftEnd);
    quick_sort_uint64(data + rightStart, pivInxs + rightStart, n - rightStart);
}


void bucketSortSerial_uint64(uint64_t *data, uint64_t *aux, int dataCount, int baseBitShift, int lgNBUCKETS, short beNested) { //vstup v aux, vystup v data

    int i;
    
    int NBUCKETS = (1 << lgNBUCKETS);
    int bitShift = 64 - lgNBUCKETS;

    int bCounts[NBUCKETS];       memset(bCounts, 0, NBUCKETS*sizeof(int));
    int cumBCounts[NBUCKETS];    memset(cumBCounts, 0, NBUCKETS*sizeof(int));//aka bucket start

   
    ////////////////////////////////////////////////////1st pass (zisti pocetnosti v aux)
    {
        for(i=0; i<dataCount; i++) {
            int inx = (aux[i] << baseBitShift) >> (bitShift);
            bCounts[inx]++;
        }
        cumBCounts[0] = 0;
        for(i=1; i<NBUCKETS; i++) {
            cumBCounts[i] = cumBCounts[i-1] + bCounts[i-1] ;
        }
    }

    ////////////////////////////////////////////////////2nd pass (prekopiruj z aux na prislusne miesto v data)
    {
        for(i=0; i<dataCount; i++) {
            int bIdx = (aux[i] << baseBitShift) >> (bitShift);
            data[cumBCounts[bIdx]] = aux[i];
            cumBCounts[bIdx]++;
        }
        for(i=NBUCKETS-1; i>0; i--) {
            cumBCounts[i] = cumBCounts[i-1];
        }
        cumBCounts[0] = 0;
    }


    ////////////////////////////////////////////////////bucket sort
    {
        if(beNested == 0){
            for(i=0; i<NBUCKETS; i++) {
               quick_sort_uint64(&data[cumBCounts[i]], (int*)&aux[cumBCounts[i]], bCounts[i]);
            }
        }else{
            for(i=0; i<NBUCKETS; i++) {
              int n = floor(log(bCounts[i] / 25) / log(2) );
              int nextShift = n < 11 ? n:12; 
              
              if(nextShift<4){
                quick_sort_uint64(&data[cumBCounts[i]], (int*)&aux[cumBCounts[i]], bCounts[i]);
                memcpy(&aux[cumBCounts[i]], &data[cumBCounts[i]], 8*bCounts[i]); //lebo vysledok cakam v globalnom data, co je teraz aux
              }
              else{
                bucketSortSerial_uint64(&aux[cumBCounts[i]], &data[cumBCounts[i]], bCounts[i], baseBitShift+64-bitShift, nextShift, 0);
              }
            }
        }
    }
}

void createSortSave(char *buffer, ThreadPart *threadParts, size_t memory_for_kmers, int myid, int core_part_num, int nthreads, int tid) {
    struct timeval t_start, t_end;
    
    size_t kmers_limit = memory_for_kmers / sizeof(uint64_t);  // max k-mers for core
    // printf("kmers_limit = %lu\n", kmers_limit);
    size_t kmers_limit_thread = kmers_limit / nthreads;        // max k-mers for thread // TODO: deleno 2 -> pre aux?
    size_t kmers_count = 0;
    uint64_t *kmers = (uint64_t*)malloc(sizeof(uint64_t) * kmers_limit_thread);
    uint64_t *aux = (uint64_t*)malloc(sizeof(uint64_t) * kmers_limit_thread);

    int k;
    int *nmers_counts = (int*)malloc(sizeof(int) * 256); // TODO 256...
    int *nmers_counted = (int*)malloc(sizeof(int) * 256); // TODO 256...
    for (k = 0; k < 256; k++) {nmers_counts[k] = 0; nmers_counted[k] = 0;}

#ifdef DEBUG
    fprintf(fdebug, "CREATE START myid=%d tid=%d\n", myid, tid);
    fflush(fdebug);
    gettimeofday(&t_start, NULL);
#endif

    createKmers(&buffer[threadParts[tid].start_index], threadParts[tid].size, 31, kmers_limit_thread, kmers, &kmers_count, nmers_counts);

#ifdef DEBUG
    gettimeofday(&t_end, NULL);
    double time_createKmers = (t_end.tv_sec - t_start.tv_sec) * 1000 + (t_end.tv_usec - t_start.tv_usec)/1000.0;  
    fprintf(fdebug, "CREATE DONE myid=%d tid=%d: kmers_limit_thread = %lu\tkmers_count = %lu\ttime = %lf\n", myid, tid, kmers_limit_thread, kmers_count, time_createKmers);
    fflush(fdebug);
#endif

#ifdef DEBUG
    fprintf(fdebug, "SORT START myid=%d tid=%d\n", myid, tid);
    fflush(fdebug);
    gettimeofday(&t_start, NULL);
#endif

    bucketSortSerial_uint64(aux, kmers, kmers_count, 0, 12, 1);
    
#ifdef DEBUG
    gettimeofday(&t_end, NULL);
    double time_sortKmers = (t_end.tv_sec - t_start.tv_sec) * 1000 + (t_end.tv_usec - t_start.tv_usec)/1000.0;  
    fprintf(fdebug, "SORT DONE myid=%d tid=%d\ttime=%lf\n", myid, tid, time_sortKmers);
    fflush(fdebug);
#endif

    saveKmers(kmers, kmers_count, myid, core_part_num * nthreads + tid);    // TODO: filename as parameter
    saveNmersCounts(nmers_counts, 256, myid, core_part_num * nthreads + tid);    // TODO: filename as parameter

    // saveKmersRaw(kmers, kmers_count, myid, core_part_num * nthreads + tid);
    // saveNmersCountsRaw(nmers_counts, 256, myid, core_part_num * nthreads + tid);    // TODO: filename as parameter
    free(nmers_counts);
    free(nmers_counted);

    free(aux);
    free(kmers);
}

int main(int argc, char *argv[]) 
{
    struct timeval t_start, t_end;
    struct timeval t_start_createSortSave, t_end_createSortSave;
    gettimeofday(&t_start, NULL);
    
    int nps, myid;
    int i, j, c, t;
    
    int ARG_CORE_THREADS = 2;
    int ARG_CORE_MEMORY = 2048; 

    opterr = 0;
    while ((c = getopt (argc, argv, "t:m:")) != -1) {
        switch (c) {
        case 't':
            //printf("SETTING THREADS!\n"); // optarg
            ARG_CORE_THREADS = atoi(optarg);
            break;
        case 'm':
            ARG_CORE_MEMORY = atoi(optarg);
            break;
        case '?':
            if (optopt == 'c')
            fprintf (stderr, "Option -%c requires an argument.\n", optopt);
            else if (isprint (optopt))
            fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else
            fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
            return 1;
        default:
            ;
        }
    }
    // Open input file
    char input_file_path_full[1024];
    sprintf(input_file_path_full, "%s%s", in_path, in_file_name);
    FILE *input_file = fopen(input_file_path_full, "r");
    if (input_file == NULL) {
        fprintf(stderr, "Error: input file missing\n");
        exit(1);
    }
    size_t input_file_size = getFileSize(input_file);

    size_t MEMORY_MB = ARG_CORE_MEMORY;
    MEMORY_MB -= 16; // nejaka rezerva // debug uncomment
    int NUM_THREADS_PER_CORE = ARG_CORE_THREADS;

    // MPI Init
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nps);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

#ifdef DEBUG
    char file_debug_path[1024];
    sprintf(file_debug_path, "%sdebug_np%d_nt%d_m%d_p%d.txt", out_path, nps, ARG_CORE_THREADS, ARG_CORE_MEMORY, myid);
    fdebug = fopen(file_debug_path, "w");
#endif

    // core required memories
    size_t max_memory = (size_t)MEMORY_MB*1024*1024;
    size_t memory_for_kmers = max_memory / 2;
    size_t memory_for_raw_part = memory_for_kmers / 8;
    
#ifdef DEBUG
    fprintf(fdebug, "max_memory = %lu\nmemory_for_kmers %lu\nmemory_for_raw_part = %lu\n", max_memory, memory_for_kmers, memory_for_raw_part);
    fflush(fdebug);
#endif
    
    FilePart *fileParts = (FilePart*)malloc(sizeof(FilePart) * nps);
    getFilePartsStartsAndEnds(input_file, input_file_size, END_READ_RESERVE, fileParts, nps, memory_for_kmers, memory_for_raw_part);
    
    fclose(input_file);
    
    // debug fileParts
    // if (myid == 0) for (i = 0; i < nps; i++) {printFilePartData(fileParts[i]);}
#ifdef DEBUG
    printFilePartDataToFile(fileParts[myid], fdebug);
    fflush(fdebug);
#endif

    // debug vypis casti
    /*
        i = myid;
        char file_test_out_path[1024];
        sprintf(file_test_out_path, "%stest_out_%d.txt", out_path, myid);
        FILE *file_test = fopen(file_test_out_path, "w");
        char *buf = (char*)malloc(sizeof(char) * fileParts[i].size); 
        fseek(input_file, fileParts[i].start_index, SEEK_SET);
        size_t r = fread(buf, sizeof(char), fileParts[i].size, input_file);
        fwrite(buf, sizeof(char), fileParts[i].size, file_test);
        fclose(file_test);
    */
    int nthreads = NUM_THREADS_PER_CORE;
    ThreadPart *threadParts;
    char *buffer = (char*)malloc(sizeof(char) * memory_for_raw_part);

    #pragma omp parallel num_threads(nthreads) shared(nthreads, threadParts, buffer, input_file_path_full) private(i)
    {
        nthreads = omp_get_num_threads();
        int tid = omp_get_thread_num();
        
        FILE *input_file = fopen(input_file_path_full, "r");
        
        for (i = 0; i < fileParts[myid].num_of_core_parts; i++) {
            CorePart corePart = fileParts[myid].coreParts[i];
            printCorePart(corePart);

#ifdef DEBUG
    printCorePartToFile(corePart, fdebug);
#endif
            threadParts = (ThreadPart*)malloc(sizeof(ThreadPart) * nthreads);  
   
            #pragma omp master
            {
                // buffer = (char*)malloc(sizeof(char) * corePart.size);
                fseek(input_file, corePart.start_index, SEEK_SET);
                size_t buf_size = fread(buffer, sizeof(char), corePart.size, input_file);
                getThreadPartsStartsAndEnds(buffer, corePart, END_READ_RESERVE, nthreads, threadParts);
                
                for (t = 0 ; t < nthreads; t++) {
                    printThreadPart(threadParts[t]);
#ifdef DEBUG
    printThreadPartToFile(threadParts[t], fdebug);
#endif
                }                
            }
            
            #ifdef DEBUG
                fprintf(fdebug, "createSortSave START myid=%d tid=%d\n", myid, tid);
                fflush(fdebug);
                gettimeofday(&t_start_createSortSave, NULL);
            #endif

            #pragma omp barrier
            createSortSave(buffer, threadParts, memory_for_kmers, myid, i, nthreads, tid);
            #pragma omp barrier // !!
            
            #ifdef DEBUG
                gettimeofday(&t_end_createSortSave, NULL);
                double time_createSortSave = (t_end_createSortSave.tv_sec - t_start_createSortSave.tv_sec) * 1000 + (t_end_createSortSave.tv_usec - t_start_createSortSave.tv_usec)/1000.0;
                fprintf(fdebug, "createSortSave DONE myid=%d tid=%d\ttime = %lf\n", myid, tid, time_createSortSave);
                fflush(fdebug);
            #endif
            
            // free(threadParts);
            // free(buffer);
        }
    }
    
    #pragma omp barrier // nanic?
    MPI_Barrier(MPI_COMM_WORLD);
    free(buffer);
    free(threadParts);
    
    int bits_in_int = 256; // TODO
    int **nmers_parts = (int**)malloc(sizeof(int*) * nps * nthreads * fileParts[myid].num_of_core_parts);     // zapamatat si pre jednotlive casti (Part <=> subory)
    // printf("|nmers_parts| = %d\n", nps * nthreads * coreParts[myid].num_of_parts);
    int *nmers_total = (int*)malloc(sizeof(int) * 256);  // zapamatat si kolko je vo vsetkych partoch (suboroch)
    for (i = 0; i < bits_in_int; i++) {nmers_total[i] = 0;}
    
    for (i = 0; i < fileParts[myid].num_of_core_parts * nps * nthreads; i++) {   // TODO: skontrolovat
        nmers_parts[i] = (int*)malloc(sizeof(int) * bits_in_int);
        for (j = 0; j < bits_in_int; j++) {
            nmers_parts[i][j] = 0;
        }
    }
    // prejst subory s nmermi
    size_t total_kmers_count = 0;
    int *nmers = (int*)malloc(sizeof(int) * bits_in_int);
    for (c = 0; c < nps; c++) {
        for (t = 0; t < nthreads * fileParts[c].num_of_core_parts; t++) {
            char filename_nmers[1024];
            sprintf(filename_nmers, "%snmers_%d_%d", out_path, c, t);
            
#ifdef DEBUG    
            fprintf(fdebug, "nmer file = %s\n", filename_nmers);
            fflush(fdebug);
#endif
            FILE *fr_nmers = fopen(filename_nmers, "r");
            
            if (fr_nmers == NULL) {
                fprintf(stderr, "FILE '%s' DOES NOT EXIST!\n", filename_nmers);
                exit(1);
            }
            fread(nmers, sizeof(int), bits_in_int, fr_nmers);
            fclose(fr_nmers);
            
            for (i = 0; i < bits_in_int; i++) {
                // printf("%d => %d\n", i, nmers[i]);
                // printf("c * nps + t = %d, i = %d\n", c * nps + t, i);
                // printf("c * fileParts[c].num_of_core_parts + t = %d, i = %d\n", c * fileParts[c].num_of_core_parts + t, i);
                // printf("===> %d\n", c * nthreads * fileParts[c].num_of_core_parts + t);
                nmers_parts[c * nthreads * fileParts[c].num_of_core_parts + t][i] = nmers[i];
                nmers_total[i] += nmers[i];
                total_kmers_count += nmers[i];
            }
        }
    }
    // COUNT PARTS
    size_t max_kmers_per_count_part = (memory_for_kmers * nthreads) / sizeof(uint64_t);
    size_t count_parts_cca = (total_kmers_count / max_kmers_per_count_part) + 1;
    
    // count parts count dividable by nps ! 
    size_t count_parts_tmp = 0;
    while (count_parts_tmp < count_parts_cca) {
        count_parts_tmp += nps;
    }
    count_parts_cca = count_parts_tmp;
    size_t balanced_kmers_per_count_part = total_kmers_count / count_parts_cca;
    // postupne pridavat 1% k balanced kmers poctu
    size_t one_percent_kmers_count = (max_kmers_per_count_part - balanced_kmers_per_count_part) / 100;
    
    while (balanced_kmers_per_count_part * count_parts_cca < total_kmers_count) {
        balanced_kmers_per_count_part += one_percent_kmers_count;
    }
    CountPart *countParts = (CountPart*)malloc(sizeof(CountPart) * count_parts_cca * 2); // *2 == rezerva
    
    int count_part_id = 0;
    size_t count_kmers = 0;
    size_t total_count_kmers = 0;
    // MPI_Finalize();
    // exit(0);

    countParts[count_part_id].id = count_part_id;
    countParts[count_part_id].start_nmer = 0;
    countParts[count_part_id].files_count = nps * nthreads * fileParts[myid].num_of_core_parts;
    countParts[count_part_id].start_offsets = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
    countParts[count_part_id].end_offsets = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
    countParts[count_part_id].files_size = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
    
    for (i = 0; i < bits_in_int; i++) {

        if (i == bits_in_int - 1) { // ked posledny n-mer je moc velky 
            // printf("[%d]: %d / %d\n", myid, count_part_id, count_parts_cca);
            
            if (count_kmers + nmers_total[i] > max_kmers_per_count_part) { // je vacsi -> treba osetrit
                // ukoncit predchadzajuci
                total_count_kmers += count_kmers;
                countParts[count_part_id].id = count_part_id;
                countParts[count_part_id].kmers_count = count_kmers;
                countParts[count_part_id].end_nmer = i - 1;
                countParts[count_part_id].files_count = nps * nthreads * fileParts[myid].num_of_core_parts;
                countParts[count_part_id].start_offsets = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
                countParts[count_part_id].end_offsets = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
                countParts[count_part_id].files_size = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
                
                // spravit posledny
                count_kmers = nmers_total[i];
                total_count_kmers += count_kmers;
                count_part_id++;
                countParts[count_part_id].start_nmer = i;
                countParts[count_part_id].end_nmer = i;
                countParts[count_part_id].id = count_part_id;
                countParts[count_part_id].kmers_count = count_kmers;
                countParts[count_part_id].files_count = nps * nthreads * fileParts[myid].num_of_core_parts;
                countParts[count_part_id].start_offsets = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
                countParts[count_part_id].end_offsets = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
                countParts[count_part_id].files_size = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
            } else { // nie je vacsia - ukoncit
                count_kmers += nmers_total[i];
                total_count_kmers += count_kmers;
                countParts[count_part_id].id = count_part_id;
                countParts[count_part_id].kmers_count = count_kmers;
                // countParts[count_part_id].end_nmer = i - 1;
                countParts[count_part_id].files_count = nps * nthreads * fileParts[myid].num_of_core_parts;
                countParts[count_part_id].start_offsets = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
                countParts[count_part_id].end_offsets = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
                countParts[count_part_id].files_size = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
                countParts[count_part_id].end_nmer = i;
            }
            
        } else if (count_kmers + nmers_total[i] > balanced_kmers_per_count_part/* || i == bits_in_int - 1 */) {
            // printf("[%d]: %d / %d\n", myid, count_part_id, count_parts_cca);
            total_count_kmers += count_kmers;
            countParts[count_part_id].id = count_part_id;
            countParts[count_part_id].kmers_count = count_kmers;
            countParts[count_part_id].end_nmer = i - 1;
            countParts[count_part_id].files_count = nps * nthreads * fileParts[myid].num_of_core_parts;
            countParts[count_part_id].start_offsets = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
            countParts[count_part_id].end_offsets = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
            countParts[count_part_id].files_size = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
            
            count_part_id++;
            count_kmers = nmers_total[i];
            countParts[count_part_id].start_nmer = i;
        } else {
            count_kmers += nmers_total[i];
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    // MPI_Finalize();
    // exit(0);

// FILES SIZE 11.4
    for (i = 0; i <= count_part_id; i++) {
        int cp_id = i;
        CountPart cp = countParts[cp_id];
        size_t offset = 0;
        
        for (c = 0; c < nps; c++) {
            int num_of_core_parts = fileParts[myid].num_of_core_parts;
            for (t = 0; t < num_of_core_parts * nthreads; t++) {
                char filename_kmers[1024];
                sprintf(filename_kmers, "%skmers_%d_%d", out_path, c, t);
                FILE *fr_kmers = fopen(filename_kmers, "r");
                
                if (fr_kmers == NULL) {
                    fprintf(stderr, "FILE '%s' DOES NOT EXIST!\n", filename_kmers);
                    exit(1);
                }
                // printf("fileSize of '%s' is %lu\n", filename_kmers, getFileSize(fr_kmers));
                // fprintf(fdebug, "fileSize of '%s' is %lu\n", filename_kmers, getFileSize(fr_kmers));
                countParts[cp_id].files_size[c * num_of_core_parts * nthreads + t] = getFileSize(fr_kmers);// TODO; getFileSize
                fclose(fr_kmers);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

//


    
    // finding cummulative counts of nmers in each file (where those nmers are, start_offsets and end_offsets)
    //int count_parts_count = count_part_id + 1; // ?
    for (i = 0; i <= count_part_id; i++) { // 3.2
        CountPart cp = countParts[i];
        
#ifdef DEBUG
        fprintf(fdebug, "count_part_id = %d: start_nmer = %lu\tend_nmer = %lu\n", i, cp.start_nmer, cp.end_nmer);
        fflush(fdebug);
#endif

        for (c = 0; c < nps; c++) {
            int num_of_core_parts = fileParts[c].num_of_core_parts;
            for (t = 0; t < nthreads * num_of_core_parts; t++) {
                size_t start_nmer = cp.start_nmer;
                size_t end_nmer = cp.end_nmer;

                cp.start_offsets[c * nthreads * num_of_core_parts + t] = 0;
                cp.end_offsets[c * nthreads * num_of_core_parts + t] = 0;
                int n;

                for (n = 0; n <= end_nmer; n++) {
                    
                    if (n < start_nmer) {
                        cp.start_offsets[c * nthreads * num_of_core_parts + t] += nmers_parts[c * nthreads * num_of_core_parts + t][n];
                        cp.end_offsets[c * nthreads * num_of_core_parts + t] += nmers_parts[c * nthreads * num_of_core_parts + t][n];
                    } else if (n <= end_nmer) {
                        cp.end_offsets[c * nthreads * num_of_core_parts + t] += nmers_parts[c * nthreads * num_of_core_parts + t][n];
                    }
                }
            }
        }   
#ifdef DEBUG
        // printCountPart(cp);
        printCountPartToFile(cp, fdebug);
#endif 
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    // MPI_Finalize();
    // exit(0);
    
    int count_parts_count = count_part_id + 1;
#ifdef DEBUG
    fprintf(fdebug, "count_parts_count = %d\n", count_parts_count);
    fflush(fdebug);
#endif
    uint64_t *kmers = (uint64_t*)malloc(max_kmers_per_count_part * sizeof(uint64_t));
    uint64_t *aux = (uint64_t*)malloc(max_kmers_per_count_part * sizeof(uint64_t));

    for (i = 0; i < count_parts_count; i++) {
        int cp_id = i;
        
        printf("[%d]: TRYING cp_id = %d\n", myid, cp_id);
#ifdef DEBUG
        fprintf(fdebug, "[%d]: TRYING cp_id = %d\n", myid, cp_id);
        fflush(fdebug);
#endif

        if (cp_id % nps != myid) { // kazdy svoju cast
            continue;            
        }
        aux = (uint64_t*)aux;
        
        printf("[%d]: DOING cp_id = %d\n", myid, cp_id);
#ifdef DEBUG
        fprintf(fdebug, "[%d]: DOING cp_id = %d\n", myid, cp_id);
        fflush(fdebug);
#endif
        
        CountPart cp = countParts[cp_id];
        size_t offset = 0;
        
        for (c = 0; c < nps; c++) {
            int num_of_core_parts = fileParts[myid].num_of_core_parts;
            for (t = 0; t < num_of_core_parts * nthreads; t++) {
                char filename_kmers[1024];
                sprintf(filename_kmers, "%skmers_%d_%d", out_path, c, t);
                FILE *fr_kmers = fopen(filename_kmers, "r");
                
                if (fr_kmers == NULL) {
                    fprintf(stderr, "FILE '%s' DOES NOT EXIST!\n", filename_kmers);
                    exit(1);
                }
                int offset_id = c * num_of_core_parts * nthreads + t; // 11.4
                fseek(fr_kmers, cp.start_offsets[offset_id] * sizeof(uint64_t), SEEK_SET);
                size_t count_to_read =  cp.end_offsets[offset_id] - cp.start_offsets[offset_id];
#ifdef DEBUG
                fprintf(fdebug, "[%d, %d]: '%s' idx = %d\tstart_offset = %lu\tend_offset = %lu\tcount_to_read = %lu\toffset_id = %d\t offset = %lu\n", c, t, filename_kmers, i, cp.start_offsets[offset_id], cp.end_offsets[offset_id], count_to_read, offset_id, offset);
                fflush(fdebug);
#endif
                // printf("[%d]: count_to_read = %d\n", myid, count_to_read);
                size_t really_read_kmers = fread(&kmers[offset], sizeof(uint64_t), count_to_read, fr_kmers);
#ifdef DEBUG
                fprintf(fdebug, "filename_kmers = '%s'\tc = %d, t = %d\t really_read_kmers = %lu\n", filename_kmers, c, t, really_read_kmers);
                fflush(fdebug);
#endif                    
                // offset += count_to_read;
                offset += really_read_kmers;
                fclose(fr_kmers);
            }
        }
        printf("[%d]: start_nmer = %d, end_nmer = %d => count_of_kmers = %lu, max kmers = %lu\n", myid, cp.start_nmer, cp.end_nmer, offset, max_kmers_per_count_part);
#ifdef DEBUG
        fprintf(fdebug, "\n[%d]: start_nmer = %d, end_nmer = %d => count_of_kmers = %lu, max kmers = %lu\n", myid, cp.start_nmer, cp.end_nmer, offset, max_kmers_per_count_part);
        fflush(fdebug);
#endif        
/*
        if (offset > max_kmers_per_count_part) {
            printf("CHYBA offset > max_kmers_per_count_part!\n >>> TODO <<<"); // TODO !!! => s balanced by sa to nemalo stat
            offset = max_kmers_per_count_part;  // TODO: inak (pridat dalsiu cast?)
        }
*/
        
        printf("bucketSort START (%d -> %d)\n", cp.start_nmer, cp.end_nmer);
#ifdef DEBUG
        fprintf(fdebug, "bucketSort START (%d -> %d)\n", cp.start_nmer, cp.end_nmer);
#endif
        bucketSortSerial_uint64(aux, kmers, offset, 0, 12, 1);
#ifdef DEBUG        
        fprintf(fdebug, "bucketSort END (%d -> %d)\n", cp.start_nmer, cp.end_nmer);
#endif
        printf("bucketSort END (%d -> %d)\n", cp.start_nmer, cp.end_nmer);


        
        printf("counting START (%d -> %d)\n", cp.start_nmer, cp.end_nmer);
#ifdef DEBUG
        fprintf(fdebug, "counting START (%d -> %d)\n", cp.start_nmer, cp.end_nmer);
#endif

// iba test ci to tu padne
size_t ii;
uint64_t kmer = kmers[0];
for (ii = 1; ii < offset; ii++) {
   kmer = kmers[ii]; 
}

/*
        char filename_kmers_counted[1024];
        sprintf(filename_kmers_counted, "%skmers_counted_%d_%d", out_path, cp.start_nmer, cp.end_nmer);
        FILE *fw_kmers_counted = fopen(filename_kmers_counted, "wb");

#ifdef RAW_FINAL_OUTPUT
        char filename_kmers_counted_raw[1024];
        sprintf(filename_kmers_counted_raw, "%skmers_counted_raw_%d_%d.txt", out_path, cp.start_nmer, cp.end_nmer);
        FILE *fw_kmers_counted_raw = fopen(filename_kmers_counted_raw, "w");
#endif
        uint64_t kmer = kmers[0];
        size_t kmer_count = 1;
        size_t max_pairs = (max_kmers_per_count_part * sizeof(uint64_t)) / sizeof(KmerCountPair); // sizeof(aux) / sizeof(KmerCountPair)

        KmerCountPair *kcp = (KmerCountPair*)aux;
        size_t pair_count = 0;
        int flushed = 0;
                
        size_t ii;
        
        for (ii = 1; ii < offset; ii++) {
                    
            if (kmers[ii] != kmer) {
                kcp[pair_count].kmer = kmer;
                kcp[pair_count].count = kmer_count;
                kmer_count = 0;
                kmer = kmers[ii];
                pair_count++;
                flushed = 0;
                        
                if (pair_count == max_pairs) {
                    fwrite(kcp, sizeof(KmerCountPair), pair_count, fw_kmers_counted);
                            
#ifdef RAW_FINAL_OUTPUT
                    for (j = 0; j < pair_count; j++) {
                        fprintf(fw_kmers_counted_raw, "%lu %d\n", kcp[j].kmer, kcp[j].count);
                    }
#endif
                    pair_count = 0;
                    flushed = 1;
                }             
            }
            kmer_count++;
        } 
                
        if (flushed == 0 && pair_count > 0) {
            fwrite(kcp, sizeof(KmerCountPair), pair_count, fw_kmers_counted);
#ifdef RAW_FINAL_OUTPUT
            for (j = 0; j < pair_count; j++) {
                fprintf(fw_kmers_counted_raw, "%lu %d\n", kcp[j].kmer, kcp[j].count);
            }
#endif
        }
                
        fclose(fw_kmers_counted);
        
#ifdef RAW_FINAL_OUTPUT
        fclose(fw_kmers_counted_raw);
#endif
*/


        printf("counting END (%d -> %d)\n", cp.start_nmer, cp.end_nmer);
#ifdef DEBUG
        fprintf(fdebug, "counting END (%d -> %d)\n", cp.start_nmer, cp.end_nmer);
#endif


    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    // free each corePart in fileParts?
    free(fileParts);
    free(countParts);
    free(kmers);
    free(aux);
    
    gettimeofday(&t_end, NULL);
    double time_total = (t_end.tv_sec - t_start.tv_sec) * 1000 + (t_end.tv_usec - t_start.tv_usec)/1000.0;  
    printf("[%d]: total time = %lf\n", myid, time_total);
#ifdef DEBUG
    fprintf(fdebug, "total time = %lf\n", myid, time_total);
    fflush(fdebug);
#endif
    fclose(fdebug);
    MPI_Finalize();
    return 0;
}