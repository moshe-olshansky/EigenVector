**GW Intra**

To create the executable do:

**g++ -O -o newGW_Intra_Flip.exe GWevIntra_new.cpp theBestEigen.c thdMul.c hgFlipSign.c straw.cpp -I. -lz -lcurl -lpthread**

Run **./newGW_Intra_Flip.exe** to see the usage.

Unless **-o observed** is specified, the **o/e** matrices will be used.

If it is run with **-v 2** flag it will output summary information for every chromosome.

Outputs the GW intrachromosomal eigenvector in wig format. Flips the sign so that the eigenvector is positive for the A compartment. Currently supports hg19 and hg38 only (other genomes are OK but no sign flipping will happen).

**Per Chromosome**

To create the executable do:

**g++ -O -o newSingleFlip.exe theEigenVector_flip_new.cpp theBestEigen.c thdMul.c hgFlipSign.c straw.cpp -I. -lz -lcurl -lpthread**

Run **./newSingleFlip.exe** to see the usage.

Unless **-o observed** is specified, the **o/e** matrices will be used.

Unlike the GW case, only the eigenvector (with NaNs) is generated.


**IMPORTANT:** by default both GW Intral and single chromosome cakculations are done on the unnormalized matrix. Specify **-n desired_norm** to run on the desired_norm normalized matrix.
