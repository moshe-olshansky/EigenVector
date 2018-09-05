require(Matrix)
source("bestEigen1.R")
args <- commandArgs(trailingOnly = TRUE)

fin <- args[1]
fout <- args[2]
binsize <- as.integer(args[3])
if (length(args) >= 4) tol  <- as.numeric(args[4])
if (length(args) >= 4) maxiter <- as.integer(args[5])

y <- scan(fin,what=list(s=integer(),e=integer(),v=double()))
y$v[is.na(y$v)] <- 0
i <- y$s/binsize+1
j <- y$e/binsize+1
x <- sparseMatrix(i=i,j=j,x=y$v,symmetric=TRUE)
if (length(args) == 3) res <- bestEigen1(x)
if (length(args) == 4) res <- bestEigen1(x,tol=tol)
if (length(args) == 5) res <- bestEigen1(x,maxiter=maxiter,tol=tol)
n <- length(res$v)
z <- cbind(as.integer(((1:n)-0.5)*binsize),res$v)
write.table(z,fout,quote=FALSE,row.names=FALSE,col.names=FALSE)
writeLines(paste("lam1 =",res$lam1,",er =",res$er,",lam2 =",res$lam2," niter=",res$niter),paste0(fout,".out"))

