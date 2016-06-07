#ifndef _PARTITION_H_
#define _PARTITION_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memset */

#include "defs.h"


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
 
typedef struct CorePart {
    int id; // inner id
    size_t start_index;
    size_t end_index;
    size_t size;
} CorePart;
 
 typedef struct ThreadPart {
	int id;
	size_t start_index; // local for buffer (not for file)
	size_t end_index;   // local for buffer (not for file)
	size_t size;
} ThreadPart;

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

typedef struct CountedParams {
    int nmers_bits;
    int counted_parts_count;
    int *start_nmers;
    int *end_nmers;
} CountedParams;

size_t getFileSize(FILE *f);
size_t getEndOfSequenceForward(char *buf, size_t length);
void getThreadPartsStartsAndEnds(char *buffer, CorePart corePart, size_t reserve, int nthreads, ThreadPart *threadParts);
void getFilePartsStartsAndEnds(FILE *f, size_t filesize, size_t reserve, FilePart *fileParts, int file_parts_count, size_t memory_for_kmers, size_t memory_for_raw_part);
/*
void setCountPart(CountPart *cp, int id, int start_nmer, int end_nmer, size_t kmers_count, int files_count, size_t *start_offsets, size_t *end_offsets);
*/

void printFilePartData(FilePart filePart);
void printCorePart(CorePart corePart);
void printCountPart(CountPart countPart);
void printThreadPart(ThreadPart threadPart);

void printCountedParams(CountedParams countedParams);
void saveCountedParams(CountedParams countedParams, char out_path[1024]);
void loadCountedParams(char filename_counted_params[1024]);

#ifdef DEBUG
void printFilePartDataToFile(FilePart filePart, FILE *f);
void printCorePartToFile(CorePart corePart, FILE *f);
void printCountPartToFile(CountPart countPart, FILE *f);
void printThreadPartToFile(ThreadPart threadPart, FILE *f);
#endif

#endif