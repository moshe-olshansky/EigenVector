library(optparse,quietly=TRUE)
library(Matrix,quietly=TRUE)
source("bestEigen3.R")

option_list <- list(
  make_option(c("-t","--tolerance"), help = "precision", type="double", default=1.0e-6),
  make_option(c("-m","--maxiter"), help = "maximum iterations", type="integer",default=NA),
  make_option(c("-s","--size"),help = "chromosome length (in basepairs)", type="integer",default=NA),
  make_option(c("-v","--verbose"),help = "verbose", type="logical",default=FALSE)
)
parser <- OptionParser(
   usage = paste("Rscript %prog [OPTIONS] inFfile outFile resolution (basepairs)",
                 "computes the principal eigenvector of the correlation matrix of HiC contacts matrix",
                 "",
                 "inFile = file containing the contacts for the chromosome (as produced by dump)", 
                 "outFile = file to output the eigenvector; outFile.report will also be produced to report run summary",
                 "resolution = resolution in basepairs", sep="\n"),
   epilogue = "inFile, outFile and resolution are required.\n",
   option_list=option_list)

arguments=NA
tryCatch(
{ arguments = parse_args(parser, positional_arguments=TRUE);},
   error = function(e) { })
if (all(is.na(arguments))) {
     stop (paste("Failed to parse command-line parameters",
                 "Use --help for help"))
}
opts = arguments$options
if (length(arguments$args) < 3) stop("Need inFile, outFile and resolution\nUse --help for help")
inFile = arguments$args[1]
outFile = arguments$args[2]
binsize = as.numeric(arguments$args[3])

verbose <- opts$verbose
tol <- as.numeric(opts$tolerance)
maxiter <- 100
if (!is.na(opts$maxiter)) maxiter <- as.numeric(opts$maxiter)

t1 <- system.time(y <- scan(inFile,what=list(s=integer(),e=integer(),v=double()),quiet=TRUE))
t1 <- t1["elapsed"]
k <- length(y$s)
if (verbose) print(paste("took",t1,"seconds to read",k,"records"),digits=6, quote=FALSE)
y$v[is.na(y$v)] <- 0
i <- y$s/binsize+1
j <- y$e/binsize+1
n <- max(j)
if (!is.na(opts$size)) n <- ceiling(as.numeric(opts$size)/binsize)
t2 <- system.time(x <- sparseMatrix(i=i,j=j,x=y$v,dims=c(n,n),symmetric=TRUE))
t2 <- t2["elapsed"]
if (verbose) print(paste("took",t2,"seconds to build sparse matrix"),digits=6, quote=FALSE)
rm(y)
rm(i)
rm(j)
for (i in 1:2) gc(verbose=FALSE)
t3 <- system.time(res <- bestEigen(x,maxiter=maxiter,tol=tol))
t3 <- t3["elapsed"]
if (verbose) {
	print(paste("took",t1,"seconds to read",k,"records"),digits=6, quote=FALSE)
	print(paste("took",t2,"seconds to build sparse matrix"),digits=6, quote=FALSE)
	print(paste("eigenvector computation took",t3,"seconds"),digits=6, quote=FALSE)
	print(paste("lam1 =",res$lam1,",er =",res$er,",lam2 =",res$lam2," niter=",res$niter),digits=6, quote=FALSE)
}	

write.table(res$v,outFile,quote=FALSE,row.names=FALSE,col.names=FALSE)
writeLines(
	   c(
	     paste("took",t1,"seconds to read",k,"records"),
	     paste("took",t2,"seconds to build sparse matrix"),
	     paste("eigenvector computation took",t3,"seconds"),
	     paste("lam1 =",res$lam1,",er =",res$er,",lam2 =",res$lam2," niter=",res$niter)
	     ),
	     paste0(outFile,".report")
	 )
