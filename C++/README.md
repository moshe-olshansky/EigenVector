Very fast and memory efficient functions for the computation of the principal eigenvector of a correlation matrix of a given matrix (usually originatimg from HiC).  
This package requires stra.cpp and straw.h from Aiden Lab ( https://github.com/aidenlab/straw/tree/master/C%2B%2B ).  

To create the executable put all the file (including straw.cpp and straw.h) in one folder and while in that folder do:  
g++ -O --std=c++0x -o theEigenVector.exe theEigenVector.cpp theBestEigen.c thMul.c straw.cpp -I. -lz -lcurl -lpthread  
This will create an executable file theEigenVector.exe.  
Run  
./theEigenVector.exe  
or  
./theEigenVector.exe -h  
to see the usage and options.
