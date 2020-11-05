#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <stdlib.h>
#include <sys/timeb.h>
#include <time.h>
#include <cmath>
#include <unistd.h>
using namespace std;

int bestEigen(long m,int *i,int *j,double *x,int *N,double *l1,double *l2,double *a1,double *a2,double *ev,double *er,double tol,int maxiter,int threads);

int getMatrix(string fname, int binsize, string norm, int **i, int **j, double **x, long *m, vector<std::string> &CH,vector<int> &CHL);

static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [-t tol][-I max_iterations][-n normalization][-T threads][-v verbose] <hicfile> <outfile> <resolution> \n", argv0);
  fprintf(stderr, "  <hicfile>: hic file\n");
  fprintf(stderr, "  <outfile>: eigenvector output  file\n");
  fprintf(stderr, "  <resolution>: resolution in bp\n");
}

int main(int argc, char *argv[]) {
        int version;
        string norm("NONE");
        string unit("BP");
        ifstream fin;
        struct timeb t0,t1;

        double tol=1.0e-6;
        int maxiter=100;
        int threads = 1;
        int verb = 1;
        int opt;

        while ((opt = getopt(argc, argv, "t:I:T:n:v:h")) != -1) {
                switch (opt) {
                                case 't':
                                        tol = atof(optarg);
                                        break;
                                case 'I':
                                        maxiter=atoi(optarg);
                                        break;
                                case 'n':
                                        norm=optarg;
                                        break;
                                case 'T':
                                        threads=atoi(optarg);
                                        break;
                                case 'v':
                                        verb=atoi(optarg);
                                        break;
                                case 'h':
                                        usage(argv[0]);
                                        exit(EXIT_SUCCESS);
                                default:
                                        usage(argv[0]);
                                        exit(EXIT_FAILURE);
                }
        }

        if (argc - optind < 3) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

	string fname = argv[optind++];
	if (verb) printf("Reading hic file from %s\n", argv[optind-1]);
	char *out_name = argv[optind++];

	FILE *fout = fopen(out_name,"w");
	if (fout==NULL) {
		fprintf(stderr, "Error! File %s cannot be opened for writing\n", argv[optind-1]);
		exit(EXIT_FAILURE);
	  }
	else if (verb) printf("Writing norm vector to %s\n", argv[optind-1]);

	int binsize = atoi(argv[optind++]);

	int *ii;
	int *jj;
	double *xx;
	long p;
	long m;
	int N;
	vector<std::string> chroms;
	vector<int> chrLen;

	ftime(&t0);
	N = getMatrix(fname, binsize, norm, &ii, &jj, &xx, &m,chroms,chrLen);
	ftime(&t1);
	if (N < 0) {
		cerr << "Error! File " << fname << " cannot be opened for reading" << endl;
		exit(EXIT_FAILURE);
	}
	if (verb) printf("took %ld seconds for %ld records\n",(long) (t1.time - t0.time),m);

	for (p=0; p<m;p++) if (jj[p] == ii[p]) xx[p] *= 0.5;

        double l1,l2,er;
        double *a1 = (double *) malloc(N*sizeof(double));
        double *a2 = (double *) malloc(N*sizeof(double));
        double *ev = (double *) malloc(N*sizeof(double));

        long init_seed = 1234563;
        srand48(init_seed);
        for (p=0;p<N;p++) {
                a1[p] = drand48();
                a2[p] = drand48();
        }

        ftime(&t0);
        int iter;
        iter = bestEigen(m,ii,jj,xx,&N,&l1,&l2,a1,a2,ev,&er,tol,maxiter,threads);
        ftime(&t1);

        if (verb) {
                printf("total %d iterations\n",iter);
                printf("iterations took %15.10lf seconds\n",((double) (t1.time - t0.time)) + 0.001*(t1.millitm - t0.millitm));
                printf("lam1 = %20.6f; lam2 = %20.6f; er = %g; error in EV = %g\n",l1,l2,er,er/(l1-l2));
        }
	fprintf(fout,"track type=wiggle_0\n");
        int n = 0;
        for (int i=0; i<chroms.size(); i++) {
       		fprintf(fout,"fixedStep ");
		char *chr = const_cast<char*> (chroms.at(i).c_str());
		char *chr1 = (char *) malloc((strlen(chr)+4)*sizeof(char));
		if (!strstr(chr,"chr")) strcpy(chr1,"chr");
		else strcpy(chr1,"");
		strcat(chr1,chr);
		if (strcmp(chr1,"chrMT") == 0)  strcpy(chr1,"chrM");
		fprintf(fout,"chrom=%s ",chr1);
       		fprintf(fout,"start=1 step=%d span=%d\n",binsize,binsize);
		for (int j=0;j< ((int) ceil(chrLen.at(i)/((double) binsize)));j++) {
			if (j == ((int) ceil(chrLen.at(i)/((double) binsize))) - 1) fprintf(fout,"fixedStep chrom=%s start=%d step=%d span=%d\n",chr1,j*binsize+1,chrLen.at(i) % binsize, chrLen.at(i) % binsize );
	        	if (!isnan(ev[n])) fprintf(fout,"%20.10f\n",ev[n++]);
                        else {
                                fprintf(fout,"%20s\n","NaN");
                                n++;
                        }
                }
        }
        fclose(fout);
        return(iter);
}

