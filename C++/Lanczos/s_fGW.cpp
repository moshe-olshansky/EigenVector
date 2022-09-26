#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <sys/timeb.h>
#include <ctime>
#include <cmath>
#include <unistd.h>
using namespace std;

int SOLan(long m,unsigned int *i,unsigned int *j,float *x,unsigned int *N,double *r,int nv,double *lam,double **ev,double *er,double tol,double eps,int maxiter,int threads);

unsigned int getMatrix(string fname, int binsize, string norm, string ob, bool interOnly, unsigned int **i, unsigned int **j, float **x, long *m, vector<std::string> &CH,vector<int> &CHL);

static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [-f (for full matrix)][-t tol][-e eps][-I max_iterations][-T threads][-v verbose] <hicfile> <outbase> <resolution> <[nv]>\n", argv0);
  fprintf(stderr, "  <hicfile>: hic file\n");
  fprintf(stderr, "  <outbase>: Eigenvector base name\n");
  fprintf(stderr, "  <resolution>: resolution in bp\n");
  fprintf(stderr, "  <number of eigenvectors>: optional; default is 2\n");
}

int main(int argc, char *argv[]) {
	string norm("NONE");
	string unit("BP");
	string ob("observed");
	ifstream fin;
	time_t t0,t1;

	double tol=1.0e-7;
	double eps=1.0e-8;
	int maxiter=200;
	int threads = 1;
	int verb = 1;
	bool interOnly = true;
	unsigned int N;
	int opt;
	int nv=2;

	while ((opt = getopt(argc, argv, "ft:e:I:T:v:h")) != -1) {
		switch (opt) {
				case 'f':
					interOnly = false;
					break;
				case 't':
					tol = atof(optarg);
					break;
				case 'e':
					eps = atof(optarg);
					break;
				case 'I':
					maxiter=atoi(optarg);
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
	if (verb) printf("Eigenvectors base name is: %s\n", out_name);
	int binsize = atoi(argv[optind++]);
	if (argc > optind) nv=atoi(argv[optind++]);

	printf("\n");

	unsigned int *i;
        unsigned int *j;
        float *x;
        long p;
        long m;
        vector<std::string> chroms;
        vector<int> chrLen;
        time(&t0);
        N = getMatrix(fname, binsize, norm, ob, interOnly, &i, &j, &x, &m,chroms,chrLen);
        time(&t1);
        if (N <= 0) {
                cerr << "Error! File " << fname << " cannot be opened for reading" << endl;
                exit(EXIT_FAILURE);
        }
        if (verb) printf("took %ld seconds for %ld records\n",(long) (t1 - t0),m);

        for (p=0; p<m;p++) if (j[p] == i[p]) x[p] *= 0.5;
        int iter;
	
        double er;
        double *r = (double *) malloc(N*sizeof(double));
        for (p=0;p<N;p++) r[p] = 1.0;
        double *lam = (double *) malloc((nv+2)*sizeof(double));
        double **ev = (double **) malloc((nv+2)*sizeof(double *));
        for (p=0;p<nv+2;p++) ev[p] = (double *) malloc(N*sizeof(double));

	time(&t0);
	iter = SOLan(m,i,j,x,&N,r,nv,lam,ev,&er,tol,eps,maxiter,threads);
	time(&t1);
	if (iter > 90000) {
		printf("return code is %d\n",iter);
		exit(EXIT_FAILURE);
	}

	if (verb) {
		printf("total %d iterations\n",iter);
		printf("iterations took %ld seconds\n",t1 - t0);
	}

	for (int j0=0;j0<nv+2;j0++) printf("%lg ",lam[j0]);
        printf("\n");
        for (int j0=0;j0<nv+2;j0++) {
                char *curout = (char *) malloc((10+strlen(out_name))*sizeof(char));
                char *temp = (char *) malloc(50);
                sprintf(temp,"_Ev%d",j0+1);
                strcpy(curout,out_name);
                strcat(curout,temp);
                strcat(curout,".wig");
                FILE *fout = fopen(curout,"w");
                if (fout==NULL) {
                        fprintf(stderr, "Error! File %s cannot be opened for writing\n", curout);
                        exit(EXIT_FAILURE);
                }

        	fprintf(fout,"track type=wiggle_0\n");
        	int a = 0;
        	for (int b=0; b<chroms.size(); b++) {
                	fprintf(fout,"fixedStep ");
                	char *chr = const_cast<char*> (chroms.at(b).c_str());
                	char *chr1 = (char *) malloc((strlen(chr)+4)*sizeof(char));
                	if (!strstr(chr,"chr")) strcpy(chr1,"chr");
                	else strcpy(chr1,"");
                	strcat(chr1,chr);
                	if (strcmp(chr1,"chrMT") == 0)  strcpy(chr1,"chrM");
                	fprintf(fout,"chrom=%s ",chr1);
                	fprintf(fout,"start=1 step=%d span=%d\n",binsize,binsize);
                	for (int c=0;c< ((int) ceil(chrLen.at(b)/((double) binsize)));c++) {
                        	if (c == ((int) ceil(chrLen.at(b)/((double) binsize))) - 1) fprintf(fout,"fixedStep chrom=%s start=%d step=%d span=%d\n",chr1,c*binsize+1,chrLen.at(b) % binsize, chrLen.at(b) % binsize );
                        	if (!std::isnan(ev[j0][a])) fprintf(fout,"%lg\n",ev[j0][a++]);
                        	else {
                                	fprintf(fout,"%s\n","0");
                                	a++;
                        	}
                	}
        	}

                fclose(fout);
        }
        return(iter);
}

