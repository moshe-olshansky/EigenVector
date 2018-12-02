require(Matrix)
source("bestEigen1.R")
args <- commandArgs(trailingOnly = TRUE)

fin <- args[1]
fout <- args[2]
binsize <- as.integer(args[3])
if (length(args) >= 4) tol  <- as.numeric(args[4])
if (length(args) >= 5) maxiter <- as.integer(args[5])

t0 <- Sys.time()
y <- scan(fin,what=list(s=integer(),e=integer(),v=double()))
#y <- scan(fin,what=list(s=double(),e=double(),v=double()))
y$v[is.na(y$v)] <- 0
i <- y$s/binsize+1
j <- y$e/binsize+1
x <- sparseMatrix(i=i,j=j,x=y$v,symmetric=TRUE)
rm(y)
for (i in 1:5) gc()
t1 <- Sys.time()
if (length(args) == 3) res <- bestEigen(x)
if (length(args) == 4) res <- bestEigen(x,tol=tol)
if (length(args) == 5) res <- bestEigen(x,maxiter=maxiter,tol=tol)
n <- length(res$v)
z <- cbind(as.integer(((1:n)-0.5)*binsize),res$v)
write.table(z,fout,quote=FALSE,row.names=FALSE,col.names=FALSE)
t2 <- Sys.time()
writeLines(paste("lam1 =",res$lam1,"\ner =",res$er,"\nlam2 =",res$lam2,"\nniter=",res$niter,"\ntime to read and create matrix =",format(t1-t0,format="%H:%M:%S"), "\ntime to compute eigenvector = ",format(t2-t1,format="%H:%M:%S")),paste0(fout,".out"))

