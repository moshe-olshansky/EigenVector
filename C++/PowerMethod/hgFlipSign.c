#include <string.h>
#include <math.h>
#include <human_100K_regions.h>

int flipSign(const char *genome, double *x, int n, char *chr, int binsize) {
	int *starts;
	int i,j,l;
	for (i=0;i<24;i++) if (strcmp(chr,hchroms[i]) == 0) break;	
	if (i == 24) return(0);
	if (strcmp(genome,"hg19") == 0) starts = hg19[i];
	else if (strcmp(genome,"hg38") == 0) starts = hg38[i];
	else return(0);
	int p=0, m=0;
	int k=100000/binsize;
	for (i=0;starts[i] >= 0; i++) {
		j = starts[i]/binsize;
		for (l=j;l<(j+k);l++) {
			if (isnan(x[j])) continue;
			if (x[j] > 0) p++;
			else m++;
		}
	}
	if (p>m) return(1);
	else {
		for (i=0;i<n;i++) x[i] = -x[i];
		return(-1);
	}
}



