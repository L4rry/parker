#include <stdio.h>
#include <kmers.h>

int main(int argc, char *argv[])
{
    FILE *fr_kmers_counted = fopen("kmers_counted_0_203", "r");
    
    if (fr_kmers_counted != NULL) {
        KmerCountPair *kcp = (KmerCountPair*)malloc(sizeof(KmerCountPair) * 1);
        // fread(kcp, sizeof(KmerCountPair), pair_count, fr_kmers_counted);
        fread(kcp, sizeof(KmerCountPair), 1, fr_kmers_counted);
        printf("%llu => %llu\n", kcp.kmer, kcp.count);
        free(kcp);
    }
    
    fclose(fr_kmers_counted);
}
    