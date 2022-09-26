#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <sys/timeb.h>
#include <time.h>
#include <cmath>
#include <unistd.h>
using namespace std;

int SOLan(long m,unsigned int *i,unsigned int *j,float *x,unsigned int *N,double *r,int nv,double *lam,double **ev,double *er,double tol,double eps,int maxiter,int threads);

map<string, chromosome> readHeader(istream &fin, long &masterIndexPosition, string &genomeID, int &numChromosomes, int &version, long &nviPosition, long &nviLength);

int flipSign(const char *genome, double *x, int n, char *chr, int binsize);

static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [-o (observed)][-t tol][-e eps][-I max_iterations][-n normalization][-T threads][-v verbose] <hicfile> <chromosome> <outbase> <resolution> <[nv]>\n", argv0);
  fprintf(stderr, "  <hicfile>: hic file\n");
  fprintf(stderr, "  <chromosome>: chromosomee\n");
  fprintf(stderr, "  <outbase>: Eigenvectors output base name\n");
  fprintf(stderr, "  <resolution>: resolution in bp\n");
  fprintf(stderr, "  <number of eigenvectors>: optional; default is 2\n");
}

int main(int argc, char *argv[]) {
	string norm("NONE");
	string unit("BP");
	string ob("oe");
	ifstream fin;
	time_t t0,t1;

	double tol=1.0e-7;
	double eps=1.0e-8;
	int maxiter=200;
	int threads = 1;
	int verb = 1;
	unsigned int N;
	int opt;

	while ((opt = getopt(argc, argv, "ot:e:I:T:n:v:h")) != -1) {
		switch (opt) {
				case 'o':
					ob = "observed";
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

	if (argc - optind < 4) {
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	string fname = argv[optind++];
	fin.open(fname, fstream::in);
	if (fin.is_open()) {
		if (verb) printf("Reading hic file from %s\n", argv[optind-1]);
	}
	else {
		fprintf(stderr, "Error! File %s cannot be opened for reading\n", argv[optind-1]);
		exit(EXIT_FAILURE);
	  }
	int nv=2;
	string chrom(argv[optind++]);
	char *out_name = argv[optind++];
	if (verb) printf("Eigenvectors base name is: %s\n", out_name);
	int binsize = atoi(argv[optind++]);
	if (argc > optind) nv=atoi(argv[optind++]);

	printf("\n");

// chromosome map for finding matrix
    long master = 0L;
    map<string, chromosome> chromosomeMap;
    string genomeID;
    int numChromosomes = 0;
    int version = 0;
    long nviPosition = 0;
    long nviLength = 0;
    long totalFileSize;

	chromosomeMap = readHeader(fin, master, genomeID, numChromosomes, version, nviPosition, nviLength);
	map<string,chromosome>::iterator itr0 = chromosomeMap.find(chrom);
	if (itr0 != chromosomeMap.end()) N = (int) ceil(itr0->second.length/((double) binsize));
	else {
		cout << "chromosome " << chrom << " is not found" << endl;
		exit(1);
	}
	fin.close();

        string hg19("hg19");
        string hg38("hg38");
        if (genomeID != hg19 && genomeID != hg38) if (verb) {
                cout << genomeID;
                cout << " is not currently supported; no sign flip will be attempted!" << endl << endl;
        }

	time(&t0);
	vector<contactRecord> records = straw(ob, norm, fname, chrom, chrom, unit, binsize);
	long nonZer = records.size();
	unsigned int *i = (unsigned int *) malloc(nonZer*sizeof(int));
	unsigned int *j = (unsigned int *) malloc(nonZer*sizeof(int));
	float *x = (float *) malloc(nonZer*sizeof(float));
	for (long k=0; k<nonZer; k++) {
				i[k] = records[k].binX/binsize; 
				j[k] = records[k].binY/binsize; 
				x[k] = (float) records[k].counts;
				if (isnan(x[k])) x[k] = 0;
	}
	time(&t1);
	if (verb) printf("took %ld seconds for %ld records\n",t1 - t0,nonZer);
	records.clear();
	records.shrink_to_fit();

	long p;
	long m = nonZer;
	for (p=0; p<m;p++) if (j[p] == i[p]) x[p] *= 0.5;
	double er;
	double *r = (double *) malloc(N*sizeof(double));
	for (p=0;p<N;p++) r[p] = 1.0;
	double *lam = (double *) malloc((nv+2)*sizeof(double));
	double **ev = (double **) malloc((nv+2)*sizeof(double *));
	for (p=0;p<nv+2;p++) ev[p] = (double *) malloc(N*sizeof(double));

	time(&t0);
	int iter;
	iter = SOLan( m,i,j,x,&N,r,nv,lam,ev,&er,tol,eps,maxiter,threads);
	time(&t1);
	if (iter > 90000) {
		printf("return code is %d\n",iter);
		exit(EXIT_FAILURE);
	}

	if (verb) {
		printf("total %d iterations\n",iter);
		printf("iterations took %ld seconds\n",t1-t0);
	}

        char *chr = const_cast<char*> (chrom.c_str());
        char *chr1 = (char *) malloc((strlen(chr)+4)*sizeof(char));
        if (!strstr(chr,"chr")) strcpy(chr1,"chr");
        else strcpy(chr1,"");
        strcat(chr1,chr);
        if (strcmp(chr1,"chrMT") == 0)  strcpy(chr1,"chrM");
       	char *genome1 = const_cast<char*> (genomeID.c_str());
        if (strcmp(chr1,"chrY")!=0 && strcmp(chr1,"chrM")!=0 && (100000 % binsize == 0)) {
		int junk;
		for (p=0;p<nv+2;p++) junk = flipSign(genome1,ev[p],N,chr1,binsize);
	}

	for (int j0=0;j0<nv+2;j0++) printf("%lg ",lam[j0]);
	printf("\n");
	for (int j0=0;j0<nv+2;j0++) {
		char *curout = (char *) malloc((10+strlen(out_name))*sizeof(char));
		char *temp = (char *) malloc(50);
		sprintf(temp,".Ev%d",j0+1);
		strcpy(curout,out_name);
		strcat(curout,temp);
		FILE *fout = fopen(curout,"w");
		if (fout==NULL) {
			fprintf(stderr, "Error! File %s cannot be opened for writing\n", curout);
			exit(EXIT_FAILURE);
		}

		for (p=0;p<N;p++) {
			if (!isnan(ev[j0][p])) fprintf(fout,"%lg\n",ev[j0][p]);
			else fprintf(fout,"%s\n","NaN");
		}
		fprintf(fout,"\n");
		fclose(fout);
	}
	return(iter);
}

