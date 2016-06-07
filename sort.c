#include "sort.h"

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