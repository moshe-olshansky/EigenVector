**A very fast and memory efficient functions for the computation of the principal eigenvectors of a correlation matrix of a given matrix (usually originatimg from HiC)**  

The **PowerMethod** folder contans a C implementation of the Power Method while the **Lanczos** methods contains the implementation of Lanczos method with Selective Ortohonalization.  
  
The C++ functions are faster than their R counterparts and can handle matrices with an arbitrary number of nonzeros - this number is only limited by the RAM size. We need slightly more than 12N bytes of RAM, where N is the number of nonzero elements is the (upper triangle of the)  matrix.
