#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

void normalize(double *x,unsigned int n) {
	unsigned int i;
	double s;
	s = 0;
	for (i=0; i<n;i++) s += x[i]*x[i];
	s = 1.0/sqrt(s);
	for (i=0;i<n;i++) x[i] *= s;
}
	
double norm(double *x,unsigned int n) {
	unsigned int i;
	double s;
	s = 0;
	for (i=0; i<n;i++) s += x[i]*x[i];
	s = sqrt(s);
	return(s);
}
	
double dot(double *x, double *y, unsigned int n) {
	unsigned int i;
	double s;
	s = 0;
	for (i=0;i<n;i++) s += x[i]*y[i];
	return s;
}

void utmvMul(unsigned int *i,unsigned int *j,float *x,long m,double *v,unsigned int k,double *res, int threads, double **space);

void fullMul(unsigned int *i,unsigned int *j,float *x,long m,double *a,unsigned int n,double *res, double *good,double *c,double *d, int threads, double **space) {
	double *temp1 = (double *) malloc(n*sizeof(double));
	double *temp2 = (double *) malloc(n*sizeof(double));
	unsigned int p;
	for (p=0;p<n;p++) temp1[p] = d[p]*a[p];
	utmvMul(i,j,x,m,temp1,n,temp2,threads,space);
	double z;
	z = dot(temp1,c,n);
	for (p=0;p<n;p++) temp2[p] -= good[p]*z; 
	utmvMul(i,j,x,m,temp2,n,temp1,threads,space);
	z = dot(good,temp2,n);
	for (p=0;p<n;p++) res[p] = d[p]*(temp1[p] - c[p]*z);
	free(temp1);
	free(temp2);
}

int SOLan(long m,unsigned int *i,unsigned int *j,float *x,unsigned int *N,double *r,int nv,double *lam,double **ev,double *er,double tol,double eps,int maxiter,int threads) {
	double *s1,*s2,*good, *sd, *d, *c;
	unsigned int n,k;
	long p;
	if ((*N) <= 0) {
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
	int kk = 0;
	for (p=0;p<n;p++) {
		if (s2[p] == 0) good[p] = 0;
		else {
			good[p] = 1;
			kk++;
		}
	}
	for (p=0;p<n;p++) {
		s1[p] /= kk;
		s2[p] /= kk;
		sd[p] = s2[p] - s1[p]*s1[p];
		if (good[p] == 1 && sd[p] > 1.0e-8*s2[p]) k++;
		else good[p] = 0;
	}
	for (p=0;p<n;p++) {
		sd[p] = sqrt(kk*sd[p]);
		if (good[p] == 1) d[p] = 1.0/sd[p];
		else d[p] = 0;
		c[p] = s1[p]*good[p];
		r[p] *= good[p];
	}

	double *al = (double *) malloc((maxiter+2)*sizeof(double));
	double *bet = (double *) malloc((maxiter+2)*sizeof(double));
	double *diag = (double *) malloc((maxiter+2)*sizeof(double));
	double *diag2 = (double *) malloc((maxiter+2)*sizeof(double));
	double **q = (double **) malloc((maxiter+2)*sizeof(double *));
	double *u = (double *) malloc(n*sizeof(double));
	double *work = (double *) malloc((maxiter+2)*(maxiter+2)*sizeof(double));
	double *z = (double *) malloc((maxiter+2)*(maxiter+2)*sizeof(double));

	q[0] = (double *) malloc(n*sizeof(double));
	for (p=0;p<n;p++) q[0][p] = 0;
	bet[0] = norm(r,n);
	lapack_int jj,info,ldz;
	char compz = 'I';
//	char compz = 'V';
	int j0;
	double dev;
	for (j0=1;j0<maxiter;j0++) {
		q[j0] = (double *) malloc(n*sizeof(double));
		for (p=0;p<n;p++) q[j0][p] = r[p]/bet[j0-1];
		fullMul(i,j,x,m,q[j0],n,u,good,c,d,threads,space);
		for (p=0;p<n;p++) r[p] = u[p] - bet[j0-1]*q[j0-1][p];
		al[j0] = dot(q[j0],r,n);
		for (p=0;p<n;p++) r[p] = r[p] - al[j0]*q[j0][p];
		bet[j0] = norm(r,n);
		if (j0 > 2) {
			double del=0;
			for (p=0;p<j0;p++) {
				diag[p] = al[p+1];
				if (diag[p] > del) del = diag[p];
				diag2[p] = bet[p+1];
			}
			for (p=0;p<j0;p++) diag[p] += 0.1*del;
			jj = j0;
			ldz = j0;
			info = 0;
			LAPACK_dpteqr(&compz,&jj,diag,diag2,z,&ldz,work,&info);
			if (info != 0) return(100000+info);
			if (j0 >= nv+2) {
			double gap = diag[nv-1] - diag[nv];
				dev = 0;
//				for (p=0;p<nv+1;p++) if (dev < fabs(z[p*j0+j0-1])) dev = fabs(z[p*j0+j0-1]);
				for (p=0;p<nv;p++) if (dev < fabs(z[p*j0+j0-1])) dev = fabs(z[p*j0+j0-1]);
				if (dev*diag[0]/gap < tol) break;
			}
			double *y = (double *) malloc(n*sizeof(double));
			for (p=0;p<j0 && p < nv+2;p++) {
				if (fabs(z[p*j0+j0-1]) < diag[0]*eps) {
					for (int ii=0; ii<n; ii++) {
						y[ii] = 0;
						for (int junk=0;junk <j0;junk++) y[ii] = y[ii] + q[junk+1][ii]*z[p*j0+junk];
					}
					normalize(y,n);
					double s = dot(r,y,n);
					for (int ii=0;ii<n;ii++) r[ii] = r[ii] - s*y[ii];
				}
			}
			bet[j0] = norm(r,n);
		}
		for (p=0;p<n;p++) r[p] = r[p]*good[p];
	}
	int over = 0;
	if (j0 == maxiter) over = 1;
	j0 -= over;
	double del=0;
        for (p=0;p<j0;p++) {
        	diag[p] = al[p+1];
                if (diag[p] > del) del = diag[p];
                diag2[p] = bet[p+1];
	}
	for (p=0;p<j0;p++) diag[p] += 0.1*del;
	jj = j0;
	ldz = j0;
	info = 0;
	LAPACK_dpteqr(&compz,&jj,diag,diag2,z,&ldz,work,&info);
	if (info != 0) return(200000+info);
	for (int ii=0;ii<nv+2;ii++) lam[ii] = diag[ii]-0.1*del;
	for (int pos=0;pos<nv+2;pos++) {
		for (p=0;p<n;p++) {
			if (good[p] == 0) {
				ev[pos][p]=NAN;
				continue;
			}
			ev[pos][p]=0;
			for (int junk=0;junk <j0;junk++) ev[pos][p] += q[junk+1][p]*z[pos*j0+junk];
		}
	}
	*er = dev;
	*N = n;
	return(j0+over);
}













				


