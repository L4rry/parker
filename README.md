# parker

#compile
mpicc -o3 -fopenmp -lm main.c sort.c kmers.c partition.c defs.h -o output


#run
mpirun -n 12 /home/3xkuban/freqCount/final5/parker -t 2 -m 2048 -k 31 -l 5 path/to/input/file


-t:		pocet threadov, defaultne 2

-m:		"dostupna" pamat, davam menej ako max aby boli tie sobory mensie a lepsie isiel sort a tak...

-k:		velkost k-tice

-l:		threshold, minimalna hranica pri ukladani (aspon takuto pocetnost musi mat k-tica aby sa ulozila)

-o:		cesta k vystupu

- na konci je cesta k vstupnemu suboru




datasety:

Muska (Drosophila Melanogaster):

/work/projects/3DNA-2013/data/freqCountData/data/SRX040485/fa_part_aa_ab_ac_ad_ae.fa

- v adresari su aj tie mensie datasety, a aj tie raw (bez hlaviciek, len ready)

Vcela (Bombus Impatiens):
/work/projects/3DNA-2013/data/freqCountData/data/bombus/frag_1.fa