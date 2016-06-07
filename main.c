#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memset */
#include <sys/time.h>
#include <mpi.h>
#include <omp.h>
#include <unistd.h> /* getopt? */
#include <ctype.h>  /* getopt? */
#include <math.h>
#include <inttypes.h>
#include <omp.h>

#include "defs.h"
#include "sort.h"
#include "kmers.h"
#include "partition.h"

#ifdef DEBUG
    FILE *fdebug;
#endif

void createSortSave(char *buffer, ThreadPart *threadParts, size_t memory_for_kmers, int kmer_length, char out_path[1024], int myid, int core_part_num, int nthreads, int tid, int nmers_int) {
    struct timeval t_start, t_end;
    
    size_t kmers_limit = memory_for_kmers / sizeof(uint64_t);  // max k-mers for core
    size_t kmers_limit_thread = kmers_limit / nthreads;        // max k-mers for thread // TODO: deleno 2 -> pre aux?
    size_t kmers_count = 0;
    uint64_t *kmers = (uint64_t*)malloc(sizeof(uint64_t) * kmers_limit_thread);
    uint64_t *aux = (uint64_t*)malloc(sizeof(uint64_t) * kmers_limit_thread);

    int k;
    int *nmers_counts = (int*)malloc(sizeof(int) * nmers_int);
    int *nmers_counted = (int*)malloc(sizeof(int) * nmers_int);
    for (k = 0; k < nmers_int; k++) {nmers_counts[k] = 0; nmers_counted[k] = 0;}

#ifdef DEBUG
    fprintf(fdebug, "CREATE START myid=%d tid=%d\n", myid, tid);
    fflush(fdebug);
    gettimeofday(&t_start, NULL);
#endif

    createKmers(&buffer[threadParts[tid].start_index], threadParts[tid].size, kmer_length, kmers_limit_thread, kmers, &kmers_count, nmers_counts, nmers_int);

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

    saveKmers(kmers, kmers_count, out_path, myid, core_part_num * nthreads + tid);
    saveNmersCounts(nmers_counts, nmers_int, out_path, myid, core_part_num * nthreads + tid);

#ifdef RAW_KMERS_NMERS
    // saveKmersRaw(kmers, kmers_count, out_path, myid, core_part_num * nthreads + tid);
    saveNmersCountsRaw(nmers_counts, nmers_int, out_path, myid, core_part_num * nthreads + tid);
#endif
    free(nmers_counts);
    free(nmers_counted);

    free(aux);
    free(kmers);
}

int main(int argc, char *argv[]) 
{
    // loadCountedParams("/Users/larry/Documents/PhD/2015-2016/program/outs/raw_part_aa/counted_params");
    // exit(0);
    
    struct timeval t_start, t_end;
    struct timeval t_start_createSortSave, t_end_createSortSave;
    gettimeofday(&t_start, NULL);
    
    int i, j, c, t;
    int nps, myid;

    // MPI Init
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nps);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    int ARG_CORE_THREADS = 2;
    int ARG_CORE_MEMORY = 2048; 
    int kmer_length = 31;
    int nmers_bits = 16;
    int nmers_int = (int)pow((double)2, (double)nmers_bits);

    char in_path[1024] = "/Users/larry/Documents/PhD/2015-2016/program/tests/raw_part_aa";
    char out_path[1024] = "/Users/larry/Documents/PhD/2015-2016/program/outs/raw_part_aa/";
    char out_filename_base[1024] = "kmers_counted";

    int min_threshold = 5; // min amount o k-mers to be saved

    opterr = 0;
    while ((c = getopt (argc, argv, "t:m:o:k:l:h")) != -1) {
        switch (c) {
        case 't':
            ARG_CORE_THREADS = atoi(optarg);
            break;
        case 'm':
            ARG_CORE_MEMORY = atoi(optarg);
            break;
        case 'o':
            sprintf(out_path, "%s", optarg);
            break;
        case 'k':
            kmer_length = atoi(optarg);
            break;
        case 'l':
            min_threshold = atoi(optarg);
            break;
        case 'h':
            if (myid == 0) {
                printf("Usage: mpirun -n <num_of_processors> [options] <input_file_path>\n\n");
                printf("Options (default value in (), *required):\n");
                printf(" -k=uint32\t *Length of k-mer (31)\n");
                printf(" -t=uint32\t Number of threads (2)\n");
                printf(" -m=uint32\t Available memory for core/processor in MB (2048)\n");
                printf(" -o=string\t Output path\n");
                printf(" -l=uint32\t Don't output k-mers with count lower thanthis value (5)\n");
            }
            MPI_Finalize();
            return 0;
            break;
        case '?':
        //     if (optopt == 'c')
        //     fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        //     else if (isprint (optopt))
        //     fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        //     else
        //     fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
            return 1;
        default:
            ;
        }
    }
    // k-mer length vs nmers_int
    // TODO: checkDirectoryExists
    // TODO: checkFileExists
    // TODO: fixPath ? => add slash /

    // nastavenie in_path    
    sprintf(in_path, "%s", argv[optind]);
    FILE *input_file = fopen(in_path, "r");
    if (input_file == NULL) {
        fprintf(stderr, "Error: input file missing\n");
        exit(1);
    }
    size_t input_file_size = getFileSize(input_file);

    size_t MEMORY_MB = ARG_CORE_MEMORY;
    MEMORY_MB -= 16; // nejaka rezerva // debug uncomment
    int NUM_THREADS_PER_CORE = ARG_CORE_THREADS;

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
    
#ifdef DEBUG
    printFilePartDataToFile(fileParts[myid], fdebug);
    fflush(fdebug);
#endif

#ifdef DEBUG
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
#endif
    int nthreads = NUM_THREADS_PER_CORE;
    ThreadPart *threadParts;
    char *buffer = (char*)malloc(sizeof(char) * memory_for_raw_part);

    #pragma omp parallel num_threads(nthreads) shared(nthreads, threadParts, buffer, in_path) private(i)
    {
        nthreads = omp_get_num_threads();
        int tid = omp_get_thread_num();
        
        FILE *input_file = fopen(in_path, "r");
        
        for (i = 0; i < fileParts[myid].num_of_core_parts; i++) {
            CorePart corePart = fileParts[myid].coreParts[i];

#ifdef DEBUG
    printCorePart(corePart);
    printCorePartToFile(corePart, fdebug);
#endif
            threadParts = (ThreadPart*)malloc(sizeof(ThreadPart) * nthreads);  
   
            #pragma omp master
            {
                // buffer = (char*)malloc(sizeof(char) * corePart.size);
                fseek(input_file, corePart.start_index, SEEK_SET);
                size_t buf_size = fread(buffer, sizeof(char), corePart.size, input_file);
                getThreadPartsStartsAndEnds(buffer, corePart, END_READ_RESERVE, nthreads, threadParts);
                
#ifdef DEBUG
    for (t = 0 ; t < nthreads; t++) {
        printThreadPart(threadParts[t]);
        printThreadPartToFile(threadParts[t], fdebug);
        }                
#endif
            }
            
            #ifdef DEBUG
                fprintf(fdebug, "createSortSave START myid=%d tid=%d\n", myid, tid);
                fflush(fdebug);
                gettimeofday(&t_start_createSortSave, NULL);
            #endif

            #pragma omp barrier
            createSortSave(buffer, threadParts, memory_for_kmers, kmer_length, out_path, myid, i, nthreads, tid, nmers_int);
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
    
    int **nmers_parts = (int**)malloc(sizeof(int*) * nps * nthreads * fileParts[myid].num_of_core_parts);     // zapamatat si pre jednotlive casti (Part <=> subory)
    // printf("|nmers_parts| = %d\n", nps * nthreads * coreParts[myid].num_of_parts);
    int *nmers_total = (int*)malloc(sizeof(int) * nmers_int);  // zapamatat si kolko je vo vsetkych partoch (suboroch)
    for (i = 0; i < nmers_int; i++) {nmers_total[i] = 0;}
    
    for (i = 0; i < fileParts[myid].num_of_core_parts * nps * nthreads; i++) {   // TODO: skontrolovat
        nmers_parts[i] = (int*)malloc(sizeof(int) * nmers_int);
        for (j = 0; j < nmers_int; j++) {
            nmers_parts[i][j] = 0;
        }
    }
    
    // prejst subory s nmermi
    // prerobit na MPI_Allreduce // kedy? // uz! // uz to prerob! // ne!
    size_t total_kmers_count = 0;
    int *nmers = (int*)malloc(sizeof(int) * nmers_int);
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
            fread(nmers, sizeof(int), nmers_int, fr_nmers);
            fclose(fr_nmers);
            
            for (i = 0; i < nmers_int; i++) {
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
    CountPart *countParts = (CountPart*)malloc(sizeof(CountPart) * count_parts_cca * 2); // aj nejaka rezerva
    
    int count_part_id = 0;
    size_t count_kmers = 0;
    size_t total_count_kmers = 0;

    countParts[count_part_id].id = count_part_id;
    countParts[count_part_id].start_nmer = 0;
    countParts[count_part_id].files_count = nps * nthreads * fileParts[myid].num_of_core_parts;
    countParts[count_part_id].start_offsets = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
    countParts[count_part_id].end_offsets = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
    countParts[count_part_id].files_size = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
    
    for (i = 0; i < nmers_int; i++) {

        if (i == nmers_int - 1) { // ked posledny n-mer je moc velky 
            count_kmers += nmers_total[i];
            total_count_kmers += count_kmers;
            countParts[count_part_id].id = count_part_id;
            countParts[count_part_id].kmers_count = count_kmers;
            countParts[count_part_id].files_count = nps * nthreads * fileParts[myid].num_of_core_parts;
            countParts[count_part_id].start_offsets = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
            countParts[count_part_id].end_offsets = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
            countParts[count_part_id].files_size = (size_t*)malloc(sizeof(size_t) * nps * nthreads * fileParts[myid].num_of_core_parts);
            countParts[count_part_id].end_nmer = i;
        } else if (count_kmers + nmers_total[i] > balanced_kmers_per_count_part/* || i == nmers_int - 1 */) {
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
    
    // 1.6
    // [CountedParams init]
    CountedParams countedParams;
    countedParams.nmers_bits = nmers_bits;
    countedParams.counted_parts_count = count_part_id + 1;
    countedParams.start_nmers = (int*)malloc(sizeof(int) * (count_part_id + 1));
    countedParams.end_nmers = (int*)malloc(sizeof(int) * (count_part_id + 1));
    // [END CountedParams init]
    
    for (i = 0; i <= count_part_id; i++) {
        CountPart cp = countParts[i];
        size_t offset = 0;

        countedParams.start_nmers[i] = cp.start_nmer;
        countedParams.end_nmers[i] = cp.end_nmer;
        
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
                // countParts[i].files_size[c * num_of_core_parts * nthreads + t] = getFileSize(fr_kmers);
                fclose(fr_kmers);
            }
        }
    }
    
    if (myid == 0) {
        saveCountedParams(countedParams, out_path);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    // finding cummulative counts of nmers in each file (where those nmers are, start_offsets and end_offsets)
    
    int count_parts_count = count_part_id + 1;
    for (i = 0; i < count_parts_count; i++) {
        CountPart cp = countParts[i];
        
#ifdef DEBUG
        fprintf(fdebug, "count_part_id = %d: start_nmer=%d\tend_nmer=%d\n", myid, i, cp.start_nmer, cp.end_nmer);
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
        printCountPartToFile(cp, fdebug);
#endif 
    }
    
#ifdef DEBUG
    gettimeofday(&t_end, NULL);
    double time_total_wc = (t_end.tv_sec - t_start.tv_sec) * 1000 + (t_end.tv_usec - t_start.tv_usec)/1000.0;  
    fprintf(fdebug, "total time without count = %lf\n", myid, time_total_wc);
    fprintf(fdebug, "count_parts_count = %d\n", count_parts_count);
    fflush(fdebug);
#endif

    uint64_t *kmers = (uint64_t*)malloc(max_kmers_per_count_part * sizeof(uint64_t));
    uint64_t *aux = (uint64_t*)malloc(max_kmers_per_count_part * sizeof(uint64_t));
    for (i = 0; i < count_parts_count; i++) {
//     	printf("[%d]: count_part_id = %d\n", myid, i);
        int cp_id = i;

//         printf("[%d]: TRYING cp_id = %d\n", myid, cp_id);
// #ifdef DEBUG
//         fprintf(fdebug, "[%d]: TRYING cp_id = %d\n", myid, cp_id);
//         fflush(fdebug);
// #endif
        if (cp_id % nps != myid) { // kazdy svoju cast
            continue;            
        }
        aux = (uint64_t*)aux;
        
//         printf("[%d]: DOING cp_id = %d\n", myid, cp_id);
// #ifdef DEBUG
//         fprintf(fdebug, "[%d]: DOING cp_id = %d\n", myid, cp_id);
//         fflush(fdebug);
// #endif

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
//                 offset += count_to_read;
								offset += really_read_kmers;
                fclose(fr_kmers);
            }
        }
#ifdef DEBUG
        printf("[%d]: start_nmer = %d, end_nmer = %d => count_of_kmers = %lu, max kmers = %lu\n", myid, cp.start_nmer, cp.end_nmer, offset, max_kmers_per_count_part);
        fprintf(fdebug, "\n[%d]: start_nmer = %d, end_nmer = %d => count_of_kmers = %lu, max kmers = %lu\n", myid, cp.start_nmer, cp.end_nmer, offset, max_kmers_per_count_part);
        fflush(fdebug);
#endif        
        if (offset > max_kmers_per_count_part) {
            printf("CHYBA offset > max_kmers_per_count_part!\n >>> TODO <<<"); // TODO !!! => s balanced by sa to nemalo stat
            offset = max_kmers_per_count_part;  // TODO: inak (pridat dalsiu cast?)
        }

#ifdef DEBUG
        printf("bucketSort START (%d -> %d)\n", cp.start_nmer, cp.end_nmer);
        fprintf(fdebug, "bucketSort START (%d -> %d)\n", cp.start_nmer, cp.end_nmer);
#endif
        bucketSortSerial_uint64(aux, kmers, offset, 0, 12, 1);
#ifdef DEBUG        
        printf("bucketSort END (%d -> %d)\n", cp.start_nmer, cp.end_nmer);
        fprintf(fdebug, "bucketSort END (%d -> %d)\n", cp.start_nmer, cp.end_nmer);
#endif

#ifdef DEBUG
        printf("counting START (%d -> %d)\n", cp.start_nmer, cp.end_nmer);
        fprintf(fdebug, "counting START (%d -> %d)\n", cp.start_nmer, cp.end_nmer);
#endif

//////////
        char filename_kmers_counted[1024];
        sprintf(filename_kmers_counted, "%s%s_%d_%d", out_path, out_filename_base, cp.start_nmer, cp.end_nmer);
        FILE *fw_kmers_counted = fopen(filename_kmers_counted, "wb");

#ifdef RAW_FINAL_OUTPUT
        char filename_kmers_counted_raw[1024];
        sprintf(filename_kmers_counted_raw, "%s%s_raw_%d_%d.txt", out_path, out_filename_base, cp.start_nmer, cp.end_nmer);
        FILE *fw_kmers_counted_raw = fopen(filename_kmers_counted_raw, "w");
#endif

				// printf("last kmer = %llu\n", kmers[offset - 1]);
        uint64_t kmer = kmers[0];
        size_t diff_kmers_count = 0;
        size_t total_kmers_counted = 0;
        size_t threshold_kmers_count = 0;
        
        size_t kmer_count = 1;
        size_t max_pairs = (max_kmers_per_count_part * sizeof(uint64_t)) / sizeof(KmerCountPair); // sizeof(aux) / sizeof(KmerCountPair)
//         max_pairs /= 2; // test 4.4
// 				size_t max_pairs = 10000000;

        KmerCountPair *kcp = (KmerCountPair*)aux;
// 				KmerCountPair *kcp = (KmerCountPair*)malloc(sizeof(KmerCountPair) * (max_pairs + 1000));
        size_t pair_count = 0;
        size_t flushed = 0;
                
        size_t ii;
       
        // printf("[%d]: offset = %lu\n", myid, offset);
				
				// size_t ttt_count = 0;				
        for (ii = 1; ii <= offset; ii++) {

						// if (kmers[ii] == ((uint64_t)18446744073709551612)) {
						// 	ttt_count++;
						// }

            if (kmers[ii] != kmer) {
                diff_kmers_count++;
                
                if (kmer_count >= min_threshold) {
                    // printf("%llu %d\n", kmer, kmer_count);
                    kcp[pair_count].kmer = kmer;
                    kcp[pair_count].count = kmer_count;
                    total_kmers_counted += kmer_count;
                    kmer_count = 0;
                    kmer = kmers[ii];
                    pair_count++;
                    flushed = 0;
                } else {
                    threshold_kmers_count++;
                }
                        
                if (pair_count == max_pairs) {
                    fwrite(kcp, sizeof(KmerCountPair), pair_count, fw_kmers_counted);                            
#ifdef RAW_FINAL_OUTPUT
                    for (j = 0; j < pair_count; j++) {
                        fprintf(fw_kmers_counted_raw, "%llu %llu\n", kcp[j].kmer, kcp[j].count);
                    }
                    fflush(fw_kmers_counted_raw);
#endif
                    pair_count = 0;
                    flushed = 1;
                }             
                // break;
            }
            kmer_count++;
        } 
        // printf("diff_kmers_count = %lu\n", diff_kmers_count);
        // printf("threshold_kmers_count = %lu\n", threshold_kmers_count);
#ifdef DEBUG
        // fprintf(fdebug, "TTT COUNT = %lu\n", ttt_count);
            fprintf(fdebug, "diff_kmers_count = %lu\n", diff_kmers_count);
            fprintf(fdebug, "threshold_kmers_count = %lu\n", threshold_kmers_count);
            fflush(fdebug);
#endif
             
        if (flushed == 0 && pair_count > 0) {
            fwrite(kcp, sizeof(KmerCountPair), pair_count, fw_kmers_counted);
#ifdef RAW_FINAL_OUTPUT
            for (j = 0; j < pair_count; j++) {
                fprintf(fw_kmers_counted_raw, "%llu %llu\n", kcp[j].kmer, kcp[j].count);
            }
						fflush(fw_kmers_counted_raw);
#endif
        }
        fclose(fw_kmers_counted);
        
#ifdef RAW_FINAL_OUTPUT
        fclose(fw_kmers_counted_raw);
#endif

#ifdef DEBUG
        printf("counting END (%d -> %d)\n", cp.start_nmer, cp.end_nmer);
        fprintf(fdebug, "counting END (%d -> %d)\n", cp.start_nmer, cp.end_nmer);
#endif


    } // [end] count_parts_count
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    // free each corePart in fileParts?
    
//     free(fileParts);
//     free(countParts);
//     free(kmers);
//     free(aux);

    // remove tmp files
    for (c = 0; c < nps; c++) {
        for (t = 0; t < nthreads * fileParts[c].num_of_core_parts; t++) {
            char filename_kmers[1024];
            char filename_nmers[1024];
            sprintf(filename_kmers, "%skmers_%d_%d", out_path, c, t);
            sprintf(filename_nmers, "%snmers_%d_%d", out_path, c, t);
            remove(filename_kmers);
            remove(filename_nmers);
        }
    }
    // save meta data ?
    // nmer size, hranice (subory)

    // final time
    gettimeofday(&t_end, NULL);
    double time_total = (t_end.tv_sec - t_start.tv_sec) * 1000 + (t_end.tv_usec - t_start.tv_usec)/1000.0;  
    printf("[%d]: total time = %lf\n", myid, time_total);

#ifdef DEBUG
    fprintf(fdebug, "total time = %lf\n", myid, time_total);
    fflush(fdebug);
    fclose(fdebug);
#endif
    MPI_Finalize();
    return 0;
}
