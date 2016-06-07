#ifndef _KMERS_H_
#define _KMERS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defs.h"

typedef struct KmerCountPair { // [16B]
    uint64_t kmer;
    uint64_t count;
} KmerCountPair;

void createKmers(char *buf, size_t buf_size, int KMER_LENGTH, size_t kmers_limit, uint64_t *kmers, size_t *kmers_count, int *nmers, int nmers_count);
void saveKmers(uint64_t *kmers, int kmers_count, char out_path[1024], int a, int b);
void saveKmersRaw(uint64_t *kmers, int kmers_count, char out_path[1024], int a, int b);
void saveNmersCounts(int *nmers_counts, int nmers_count, char out_path[1024], int a, int b);
void saveNmersCountsRaw(int *nmers_counts, int nmers_count, char out_path[1024], int a, int b);

#endif