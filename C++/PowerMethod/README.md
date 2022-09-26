Very fast and memory efficient functions for the computation of the principal eigenvector of a correlation matrix of a given matrix (usually originatimg from HiC).  

There are two ways to read the input file: either from hic file or from an ASCII file produced by juicer tools dump.

**Reading from hic file:**  
This package requires straw.cpp and straw.h from Aiden Lab ( https://github.com/aidenlab/straw/tree/master/C%2B%2B ).  

To create the executable, make sure that you are using version 7 of gcc and g++, go to the PowerMethod folder and do:    
**g++ -O2 -o ev.exe theEigenVector_flip_new.cpp theBestEigen.c thdMul.c hgFlipSign.c STRAW/C++/straw.cpp -I . -I STRAW/C++ -lz -lcurl -lpthread**    
where STRAW is the folder containing straw (May 2022 version).  
This will create an executable file **ev.exe**.  
Run  
./ev.exe  
or  
./ev.exe -h  
to see the usage and options. This way you can create an eigenvector for a particular chromosome.  

To produce eigenvector for every chromosome in the hic file in the wiggle (wig) format, do:  
**g++ -O2 -o GWevIntra.exe GWevIntra_new.cpp theBestEigen.c thdMul.c hgFlipSign.c STRAW/C++/straw.cpp -I . -I STRAW/C++ -lz -lcurl -lpthread**  
Run  
./GWevIntra.exe  
to see te usage  

**Note:** if run on hg19 or hg38 genome with a resolution which is a divisor of 100,000, these two functions will flip the sign of the eigenvector if necessary so that positive sign corresponds to A and negative sighn to B compartment.  

**Reading from file produced by juicer tools dump**  
First of all run juicer tools dump to produce the file, say **dump.out**. Here you can choose either observed or oe (observed over expected) and choose the desired normalization.

To create the ececutable do:  
**g++ -O2 -o mainEigen.exe mainEigen.c theBestEigen.c thdMul.c -lpthread**  
this will create an executable file **mainEigen.exe**.  
Run  
./mainEigen.exe  
or  
./mainEigen.exe -h  
to see the usage and options.  
Since we are not reading hic file there are several new options:  
**-m**  an upper bound on the number of records (can get the exact value by 'wc -l infile')  
**-s**  the size of the chromosome (in bp); if specified, the matrix is assumed to be of dimension n by n where n is (chomosome size)/binsize rounded up; if not specified, n will be inferred from the data and may be slightly smaller than the true one.
