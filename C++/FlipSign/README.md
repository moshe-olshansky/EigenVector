To create the executable do:

**g++ -O --std=c++0x -o GWevIntra_special.exe GWevIntra_special.cpp getHiCInfo.cpp theBestEigen.c thdMul.c hgFlipSign.c straw.cpp -I. -lz -lcurl -lpthread**

Run **./GWevIntra_special.exe** to see the usage.

Unless **-o observed** is specified, the **o/e** matrices will be used.

If it run with **-v 2** flag it will output summary information for every chromosome.

Outputs the GW intrachromosomal eigenvector in wig format. Flips the sign so that the eigenvector is positive for the A compartment. Currently supports hg19 and hg38 only (other genomes are OK but no sign flipping will happen).
