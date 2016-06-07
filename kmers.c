#include "kmers.h"

// help function
void uint64_to_binary(uint64_t x) {
    char b[65];
    b[0] = '\0';

    uint64_t z = (uint64_t)4611686018427387904 << 1;
    for (z = (uint64_t)4611686018427387904 << 1; z > 0; z >>= 1) {
        strcat(b, ((x & z) == z) ? "1" : "0");
    }
    printf("%s\n", b);
}

void createKmers(char *buf, size_t buf_size, int KMER_LENGTH, size_t kmers_limit, uint64_t *kmers, size_t *kmers_count, int *nmers, int nmers_int) {
    char BASES[256];
	BASES[65]=0; BASES[67]=1; BASES[71]=2; BASES[84]=3;

	int i;
 	char BASES_PARSING[256]; // check if character is a base (equals 1)
 	for (i = 0; i < 256; i++) { BASES_PARSING[i]=0; }
 	BASES_PARSING[65]=1; BASES_PARSING[67]=1; BASES_PARSING[71]=1; BASES_PARSING[84]=1;

	uint64_t BASES_SHIFTED[64][85];
	for (i = 0; i < 64; i += 2) {
		BASES_SHIFTED[i][65] = (uint64_t)0;
		BASES_SHIFTED[i][67] = (uint64_t)1 << i;
		BASES_SHIFTED[i][71] = (uint64_t)2 << i;
		BASES_SHIFTED[i][84] = (uint64_t)3 << i;
	}

    for (i = 0; i < nmers_int; i++) {
        nmers[i] = 0; 
    }
    
	(*kmers_count) = 0;
	int kmer_length_counter = 0;
    int nmers_offset = (64 - log2(nmers_int));
	uint64_t kmer = 0;
	uint64_t nmers_mask = ((uint64_t)(nmers_int - 1) << nmers_offset);
	for (i = 0; i < buf_size; i++) {
			
			if ((BASES_PARSING[ buf[i] ]) && (*kmers_count) < kmers_limit) {

					if (kmer_length_counter < KMER_LENGTH) {
                        kmer |= BASES_SHIFTED[ 62 - kmer_length_counter * 2 ][buf[i]];
                        kmer_length_counter++;

                        if (kmer_length_counter == KMER_LENGTH) {
                            kmers[(*kmers_count)] = kmer;
                            nmers[(int)(kmer >> nmers_offset)]++;
                            (*kmers_count)++;
                        }
					} else {
                        kmer = (uint64_t)((kmer << 2) | BASES_SHIFTED[2][buf[i]]);
                        kmers[(*kmers_count)] = kmer;
                        nmers[(int)(kmer >> nmers_offset)]++;
                        (*kmers_count)++;
					}
			} else {
                kmer = 0;
                kmer_length_counter = 0;
			}
	
			if ((*kmers_count) >= kmers_limit) {
				break;
			}
	}
}

void saveKmers(uint64_t *kmers, int kmers_count, char out_path[1024], int a, int b) {
    char filename_out[1024];
    sprintf(filename_out, "%skmers_%d_%d", out_path, a, b);
    FILE *fw = fopen(filename_out, "wb");
    fwrite(kmers, sizeof(uint64_t), kmers_count, fw);
    fclose(fw);
}

void saveKmersRaw(uint64_t *kmers, int kmers_count, char out_path[1024], int a, int b) {
    int i;
    char filename_kmers[1024];
    sprintf(filename_kmers, "%skmers_raw_%d_%d.txt", out_path, a, b);
    FILE *fw_kmers = fopen(filename_kmers, "w");
    for (i = 0; i < kmers_count; i++) {
        fprintf(fw_kmers, "%llu\n", kmers[i]);
    }                        
    fclose(fw_kmers);
}

void saveNmersCounts(int *nmers_counts, int nmers_count, char out_path[1024], int a, int b) {
    char filename_nmers_bin[1024];
    sprintf(filename_nmers_bin, "%snmers_%d_%d", out_path, a, b);
    FILE *fw_nmers_bin = fopen(filename_nmers_bin, "wb");
    fwrite(nmers_counts, sizeof(int), nmers_count, fw_nmers_bin);
    fclose(fw_nmers_bin);
}

void saveNmersCountsRaw(int *nmers_counts, int nmers_count, char out_path[1024], int a, int b) {
    int i;
    char filename_nmers[1024];
    sprintf(filename_nmers, "%snmers_%d_%d.txt", out_path, a, b);
    FILE *fw_nmers = fopen(filename_nmers, "w");
    for (i = 0; i < nmers_count; i++) {
        fprintf(fw_nmers, "%d => %d\n", i, nmers_counts[i]);
    }                        
    fclose(fw_nmers);
}