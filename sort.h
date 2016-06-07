#ifndef _SORT_H_
#define _SORT_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memset */
#include <math.h>

#include "defs.h"

void insertion_sort_uint64(uint64_t *a, int n);
void quick_sort_uint64 (uint64_t *data, int *pivInxs, int n);
void bucketSortSerial_uint64(uint64_t *data, uint64_t *aux, int dataCount, int baseBitShift, int lgNBUCKETS, short beNested);

#endif