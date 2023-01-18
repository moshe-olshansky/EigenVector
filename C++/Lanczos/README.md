**Lanczos Method with Selective Ortohanalization**  

These functions need liblapacke (and hence libblas and liblapack). To install them do:  
sudo apt update  
sudo apt install libopenblas-base  
sudo apt install libopenblas-dev  
sudo apt install liblapack3  
sudo apt install liblapack-dev  
sudo apt install liblapacke-dev  

Then make sure that you are using version 7 of gcc and g++  

**Rmark:** In my case the folder containing **straw** is **~/HiC/straw_may_2022**. The user should replace it by their own folder.

To create an executable for computing few leading eigenvectors of the correlation matrix of contact matrix for a particular chromosome do:  
**g++ -O2 -o Lan.exe s_fLan.cpp s_fSOLan.c s_dthMul.c hgFlipSign.c ~/HiC/straw_may_2022/C++/straw.cpp -I. -I ~/HiC/straw_may_2022/C++ -lz -lcurl -lpthread -lblas -llapack -llapacke**  
Run  
./Lan.exe  
to see usage.  
By default it uses unnormalized observed over expected (o/e) matrix.  
Use -o flag to use observed matrix instead (usually not recommended)  
Use -n norm to use normalized matrix; norm can be NONE (no normalization - default), VC, VC_SQRT, KR, SCALE, SCALA, etc.  

To do the above for Genome Wide (GW) contact matrix do:  
**g++ -O2 -o GWev.exe s_fGW.cpp getGWMatrix.cpp s_fSOLan.c s_dthMul.c ~/HiC/straw_may_2022/C++/straw.cpp -I ~/HiC/straw_may_2022/C++ -lz -lcurl -lpthread -lblas -llapack -llapacke**  
Run  
./GWev.exe for usage.  
By default it uses **inter**chromosomal matrix. To use the full matrix specife the -f flag
