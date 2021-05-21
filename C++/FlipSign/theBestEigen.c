#include <stdlib.h>
#include <math.h>

void norm(double *x,int n) {
	int i;
	double s;
	s = 0;
	for (i=0; i<n;i++) s += x[i]*x[i];
	s = 1.0/sqrt(s);
	for (i=0;i<n;i++) x[i] *= s;
}
	
double dot(double *x, double *y, int n) {
	int i;
	double s;
	s = 0;
	for (i=0;i<n;i++) s += x[i]*y[i];
	return s;
}

void utmvMul(int *i,int *j,double *x,long m,double *v,int k,double *res, int threads, double **space);

void fullMul(int *i,int *j,double *x,long m,double *a,int n,double *res, double *good,double *r,double *c,double *d, int threads, double **space) {
	double *temp1 = (double *) malloc(n*sizeof(double));
	double *temp2 = (double *) malloc(n*sizeof(double));
	int p;
	for (p=0;p<n;p++) temp1[p] = d[p]*a[p];
	utmvMul(i,j,x,m,temp1,n,temp2,threads,space);
	double z;
	z = dot(a,c,n);
	for (p=0;p<n;p++) a[p] = temp2[p] - good[p]*z; 
	utmvMul(i,j,x,m,a,n,temp1,threads,space);
	z = dot(good,a,n);
	for (p=0;p<n;p++) res[p] = d[p]*temp1[p] - c[p]*z;
	free(temp1);
	free(temp2);
}

int bestEigen(long m,int *i,int *j,double *x,int *N,double *l1,double *l2,double *a1,double *a2,double *ev,double *er,double tol,int maxiter,int threads) {
	double *s1,*s2,*good, *sd, *d, *c;
	int n,k;
	long p;
	if ((*N) < 0) {
		n = 0;
		for (p=0;p<m;p++) if (j[p] > n) n = j[p];
	}
	else n = *N;
	s1 = (double *) malloc(n*sizeof(double));
	s2 = (double *) malloc(n*sizeof(double));
	good = (double *) malloc(n*sizeof(double));
	sd = (double *) malloc(n*sizeof(double));
	d = (double *) malloc(n*sizeof(double));
	c = (double *) malloc(n*sizeof(double));
	double **space = (double **) malloc(threads*sizeof(double *));
	for (p=0;p<threads;p++) space[p] = (double *) malloc(n*sizeof(double));

	for (p=0;p<n;p++) {
		s1[p] = 0;
		s2[p] = 0;
	}
	for (p=0; p<m; p++) {
		s1[i[p]] += x[p];
		s1[j[p]] += x[p];
		s2[i[p]] += x[p]*x[p];
		s2[j[p]] += x[p]*x[p];
		if (i[p] == j[p]) s2[i[p]] += 2.0*x[p]*x[p];
	}

	k = 0;
	for (p=0;p<n;p++) {
		if (s1[p] == 0) good[p] = 0;
		else {
			good[p] = 1;
			k++;
		}
	}

	for (p=0;p<n;p++) {
		s1[p] /= k;
		s2[p] /= k;
		sd[p] = sqrt(s2[p] - s1[p]*s1[p]);
		if (good[p] == 1) d[p] = 1.0/sd[p];
		else d[p] = 0;
		c[p] = s1[p]*d[p];
	}

	double *r = (double *) malloc(n*sizeof(double));
	double *v = (double *) malloc(n*sizeof(double));
	double *ev1 = (double *) malloc(n*sizeof(double));
	double *ev2 = (double *) malloc(n*sizeof(double));
	double lam1, lam2;

	utmvMul(i,j,x,m,good,n,r,threads,space);

	int niter=0;
	double er1;

	while (niter < maxiter) {
		for (int junk=0;junk<10;junk++) {
			fullMul(i,j,x,m,a1,n,a1,good,r,c,d,threads,space);
			norm(a1,n);
		}
		for (p=0;p<n;p++) ev1[p] = a1[p];
		fullMul(i,j,x,m,a1,n,a1,good,r,c,d,threads,space);
		lam1 = sqrt(dot(a1,a1,n));
		for (p=0;p<n;p++) v[p] = a1[p] - lam1*ev1[p];
		er1 = sqrt(dot(v,v,n));
		norm(a1,n);

		double trash;
		trash = dot(a1,a2,n);
		for (p=0;p<n;p++) a2[p] = a2[p] - trash*a1[p];

		 for (int junk=0;junk<10;junk++) {
			fullMul(i,j,x,m,a2,n,a2,good,r,c,d,threads,space);
			trash = dot(a1,a2,n);
			for (p=0;p<n;p++) a2[p] = a2[p] - trash*a1[p];
                        norm(a2,n);
                }

		for (p=0;p<n;p++) ev2[p] = a2[p];
		fullMul(i,j,x,m,a2,n,a2,good,r,c,d,threads,space);
                lam2 = sqrt(dot(a2,a2,n));
                norm(a2,n);

		if (lam1 < lam2) {
			double *temp = a1;
                        a1 = a2;
                        a2 = temp;
                        double tmp = lam1;
                        lam1 = lam2;
                        lam2 = tmp;
		}
		
		niter += 10;
                if (er1/(lam1-lam2) < tol) break;
        }

	for (p=0;p<n;p++) ev[p] = (good[p] == 1)? ev1[p]: NAN;
	*er = er1/n;
	*l1 = lam1/n;
	*l2 = lam2/n;
	*N = n;
	return(niter);
}













				


