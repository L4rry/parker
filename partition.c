#include "partition.h"

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
    }
    free(buf);
    fseek(f, 0, SEEK_SET);
}
/*
void setCountPart(CountPart *cp, int id, int start_nmer, int end_nmer, size_t kmers_count, int files_count, size_t *start_offsets, size_t *end_offsets) {
typedef struct CountPart {
    int id;
    size_t start_nmer;
    size_t end_nmer;
    size_t kmers_count;
    int files_count;
    size_t *start_offsets; // v kazdom subore kde je zaciatok (kolko k-tic preskocit)
    size_t *end_offsets; // v kazdom subore kde je zaciatok (kolko k-tic preskocit - kde konci)
    // size_t *files_size;
} CountPart;    
}
*/

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

void printCountedParams(CountedParams countedParams) {
    printf("nmers_bits = %d\n", countedParams.nmers_bits);
    printf("counted_parts_count = %d\n", countedParams.counted_parts_count);
    int i;
    for (i = 0; i < countedParams.counted_parts_count; i++) {
        printf("%d: %d => %d\n", i, countedParams.start_nmers[i], countedParams.end_nmers[i]);
    }
}

void saveCountedParams(CountedParams countedParams, char out_path[1024]) {
    char filename_out[1024];
    sprintf(filename_out, "%scounted_params", out_path);
    FILE *f_counted_params = fopen(filename_out, "wb");
    fwrite(&countedParams.nmers_bits, sizeof(int), 1, f_counted_params);
    fwrite(&countedParams.counted_parts_count, sizeof(int), 1, f_counted_params);
    int i;
    for (i = 0; i < countedParams.counted_parts_count; i++) {
        fwrite(&countedParams.start_nmers[i], sizeof(int), 1, f_counted_params);
        fwrite(&countedParams.end_nmers[i], sizeof(int), 1, f_counted_params);
    }
    fclose(f_counted_params);
}

void loadCountedParams(char filename_counted_params[1024]) {
    FILE *f_counted_params = fopen(filename_counted_params, "rb");
    CountedParams countedParams;
    fread(&countedParams.nmers_bits, sizeof(int), 1, f_counted_params);
    fread(&countedParams.counted_parts_count, sizeof(int), 1, f_counted_params);
    countedParams.start_nmers = (int*)malloc(sizeof(int) * (countedParams.counted_parts_count));
    countedParams.end_nmers = (int*)malloc(sizeof(int) * (countedParams.counted_parts_count));
    int i;
    for (i = 0; i < countedParams.counted_parts_count; i++) {
        fread(&countedParams.start_nmers[i], sizeof(int), 1, f_counted_params);
        fread(&countedParams.end_nmers[i], sizeof(int), 1, f_counted_params);
    }
    fclose(f_counted_params);
    printCountedParams(countedParams);
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