#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>

struct th2 {
	int *i;
	int *j;
	double *x;
	double *v;
	double *res;
	long m;
};

void *Mul(void *threadid) {
	struct th2 *a = (struct th2 *) threadid;
	int *i = a->i;
	int *j = a->j;
	double *x = a->x;
	double *v = a->v;
	double *res = a->res;
	long m = a->m;
	long p;
	for (p=0;p<m;p++) {
               res[i[p]] += x[p]*v[j[p]];
               res[j[p]] += x[p]*v[i[p]];
	}
}

void utmvMul(int *i,int *j,double *x,long m,double *v,int k,double *res, int nth, double **rs) {
	long *mm = (long *) malloc(nth*sizeof(long));
	int t;
	long n = (long) (m/nth);
	for (t=0;t<nth;t++) mm[t] = n;
	mm[nth-1] = m - n*(nth-1);
	long p;
	struct th2 a[nth];
	for (t=0;t<nth;t++) {
		for (p=0;p<k;p++) rs[t][p] = 0;
		a[t].m = mm[t];
		a[t].i = i + t*n;
		a[t].j = j + t*n;
		a[t].x = x + t*n;
		a[t].v = v;
		a[t].res = rs[t];
	}

	pthread_t threads[nth];
	int rc[nth];
	for (t=0;t<nth;t++) rc[t] = pthread_create(&threads[t], NULL, Mul, (void *) &a[t]);
	for (t=0;t<nth;t++) pthread_join(threads[t], NULL);
	for (p=0;p<k;p++) {
		res[p] = 0;
		for (t=0;t<nth;t++) res[p] += rs[t][p];
	}
}

