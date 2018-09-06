The main function is bestEigen1.R. It takes uo to 5 arguments.  
x - the square contacts matrix in Matrix (sparse) format; REQUIRED  
a1,a2 - two vectors - approximations to the first and second eigenvectors; do not use unless you know what you are doing  
tol - relative error allowed in the eigenvector (default 1.0e-6)  
maxiter - maximum number of iteretions (default 100)  

The second function is bestEigRscript_h5_1.R - it is a wrapper for bestEigen built to run via Rscript so that it can be run in parallel on several chromosomes. 
A typical run is  
Rscript --vanilla bestEigRscript_h5_1.R fin fout binsize tol maxiter  
Only first 3 arguments are required. They must come in exactly this order.  
fin - input file as produced by juicer dump function (in sparse form); REQUIRED  
fout - file to store the eigenvector; Required  
binsize - binsize (in base-pairs) used with juicer dump command; REQUIRED  
tol - desired relative error in eigenvector (default 1.0e-6)  
maxiter - maximum numberof iterations (default 100)  

The run produces two files: fout and fout.out. The first one holds the eigenvector with the midpoints of the corresponding bins. The second one contains one line describing run statistics:  
lam1 - the approximation to the first eigenvalue  
lam2 - the approximation to the second eigenvalue  
er - L2 norm of C*v - lam1*v where v is the eigenvector (normalised to have L2 norm of 1) and C is the correlation matrix  
niter - number of iterations (comes in units of 10)  

Please note that a bound on the relative error in v is er/(lam1 - lam2) (if lam1 and lam2 were exact). We stop when this becomes < tol or we reach maxiter iterations. Since each time we do 10 iterations we may exceed maxiter to the nearest multiple of 10.

PLEASE NOTE:
You need to have Matrix R package installed.  
We do not compute the correlation matrix since we only need the result of it's multiplication with a vector. We manage to keep only one matrix in sparse Matrix format so the dimension of the matrix almost does not affect the performance. What matters is the number of non-zero elements in the contact matrix produced by juicer dump command.  
Do not attempt to use tolerance (tol) smaller than 1.0e-15 since it is comparable to the roundoff error and the process may never converge.

Author: Moshe Olshansky;  e-mail: moshe.olshansky@monash.edu
 
