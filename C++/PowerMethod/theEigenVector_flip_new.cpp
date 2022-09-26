#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <time.h>
#include <cmath>
#include <unistd.h>
using namespace std;

int bestEigen(long m,int *i,int *j,double *x,int *N,double *l1,double *l2,double *a1,double *a2,double *ev,double *er,double tol,int maxiter,int threads);

int flipSign(const char *genome, double *x, int n, char *chr, int binsize);

map<string, chromosome> readHeader(istream &fin, long &masterIndexPosition, string &genomeID, int &numChromosomes,
                                   int &version, long &nviPosition, long &nviLength);

static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [-o (observed)][-t tol][-I max_iterations][-n normalization][-T threads][-v verbose] <hicfile> <chromosome> <outfile> <resolution> \n", argv0);
  fprintf(stderr, "  <hicfile>: hic file\n");
  fprintf(stderr, "  <chromosome>: chromosomee\n");
  fprintf(stderr, "  <outfile>: eigenvector output  file\n");
  fprintf(stderr, "  <resolution>: resolution in bp\n");
}

int main(int argc, char *argv[]) {
	string norm("NONE");
	string unit("BP");
	string ob("oe");
	ifstream fin;
	struct timespec t0,t1;

	double tol=1.0e-6;
	int maxiter=100;
	int threads = 1;
	int verb = 1;
	int N;
	int opt;

	while ((opt = getopt(argc, argv, "ot:I:T:n:v:h")) != -1) {
		switch (opt) {
				case 'o':
					ob = "observed";
					break;
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
	string chrom(argv[optind++]);
	char *out_name = argv[optind++];
	int binsize = atoi(argv[optind++]);
	FILE *fout = fopen(out_name,"w");
	if (fout==NULL) {
		fprintf(stderr, "Error! File %s cannot be opened for writing\n", argv[optind-1]);
		exit(EXIT_FAILURE);
	}
	else {
		if (verb) printf("Writing eigenvector to %s\n", out_name);
	}


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

	clock_gettime(CLOCK_REALTIME,&t0);
	vector<contactRecord> records = straw(ob, norm, fname, chrom, chrom, unit, binsize);
	long nonZer = records.size();
	int *ii = (int *) malloc(nonZer*sizeof(int));
	int *jj = (int *) malloc(nonZer*sizeof(int));
	double *xx = (double *) malloc(nonZer*sizeof(double));
	for (long k=0; k<nonZer; k++) {
				ii[k] = records[k].binX/binsize; 
				jj[k] = records[k].binY/binsize; 
				xx[k] = (double) records[k].counts;
				if (isnan(xx[k])) xx[k] = 0;
	}
	clock_gettime(CLOCK_REALTIME,&t1);
	if (verb) printf("took %10.3f seconds for %ld records\n",((double) (t1.tv_sec - t0.tv_sec)) + ((double) (t1.tv_nsec - t0.tv_nsec))/1e9,nonZer);
	records.clear();
	records.shrink_to_fit();

	long p;
	long m = nonZer;
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

	clock_gettime(CLOCK_REALTIME,&t0);
	int iter;
	iter = bestEigen(m,ii,jj,xx,&N,&l1,&l2,a1,a2,ev,&er,tol,maxiter,threads);
	clock_gettime(CLOCK_REALTIME,&t1);

	if (verb) {
		printf("total %d iterations\n",iter);
		printf("iterations took %10.3f seconds\n",((double) (t1.tv_sec - t0.tv_sec)) + ((double) (t1.tv_nsec - t0.tv_nsec))/1e9);
		printf("lam1 = %g; lam2 = %g; lam1/lam2 = %g; er = %g; error in EV = %g\n",l1,l2,l1/l2,er,er/(l1-l2));
	}

	        char *chr = const_cast<char*> (chrom.c_str());

              char *chr1 = (char *) malloc((strlen(chr)+4)*sizeof(char));
              if (!strstr(chr,"chr")) strcpy(chr1,"chr");
              else strcpy(chr1,"");
              strcat(chr1,chr);
              if (strcmp(chr1,"chrMT") == 0)  strcpy(chr1,"chrM");

	char *genome1 = const_cast<char*> (genomeID.c_str());
        if (strcmp(chr1,"chrY")!=0 && strcmp(chr1,"chrM")!=0 && (100000 % binsize == 0)) int junk = flipSign(genome1,ev,N,chr1,binsize);

	for (int j=0;j<N;j++) {
		if (!isnan(ev[j])) fprintf(fout,"%20.10f\n",ev[j]);
		else fprintf(fout,"%20s\n","NaN");
	}
	fclose(fout);
	return(iter);
}

