Very fast and memory efficient functions for the computation of the principal eigenvector of a correlation matrix of a given matrix (usually originatimg from HiC).  

There are two ways to read the input file either from hic file or from an ASCII file produced by juicer tools dump.

**Reading from hic file:**  
This package requires straw.cpp and straw.h from Aiden Lab ( https://github.com/aidenlab/straw/tree/master/C%2B%2B ).  

To create the executable put all the file (including straw.cpp and straw.h) in one folder and while in that folder do:    
**g++ -O --std=c++0x -o theEigenVector.exe theEigenVector.cpp theBestEigen.c thMul.c getMatrix.cpp straw.cpp -I. -lz -lcurl -lpthread**  
This will create an executable file **theEigenVector.exe**.  
Run  
./theEigenVector.exe  
or  
./theEigenVector.exe -h  
to see the usage and options. This way you can create an eigenvector for a particular chromosome.  

To create **Genome Wide** eigenvector do:  
**g++ -O --std=c++0x -o createGWEigenVector.exe createGWEigenVector.cpp theBestEigen.c thMul.c getMatrix.cpp straw.cpp -I. -lz -lcurl -lpthread**  
This will create an executable file **createGWEigenVector.exe**.  
Run  
./createGWEigenVector.exe  
or  
./createGWEigenVector.exe -h  
to see the usage and options.  
It produces a Genome Wide eigenvector in Wiggle (**wig**) format.

**Reading from file produced by juicer tools dump**  
First of all run juicer tools dump to produce the file, say **dump.out**. Here you can chhoose either observed or oe (observed over expected) and choose the desired normalization.

To create the ececutable do:  
**g++ -O --std=c99 -o mainEigen.exe mainEigen.c theBestEigen.c thdMul.c -lpthread**  
this will create an executable file **mainEigen.exe**.  
Run  
./mainEigen.exe  
or  
./mainEigen.exe -h  
to see the usage and options.  
Since we are not reading hic file there are several new options:  
**-m**  an upper bount on the number of records (can get the exact value by 'wc -l infile')  
**-s**  the size of the chromosome (in bp); if specified the matrix is assumed to be of dimension n by n where n is (chomosome size)/binsize rounded up; if not specified, n will be inferred from the data and may be slightly smaller than the true one.
