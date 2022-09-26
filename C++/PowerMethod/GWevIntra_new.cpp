#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <time.h>
#include <cmath>
#include <unistd.h>
using namespace std;

int bestEigen(long m,int *i,int *j,double *x,int *N,double *l1,double *l2,double *a1,double *a2,double *ev,double *er,double tol,int maxiter,int threads);

map<string, chromosome> readHeader(istream &fin, long &masterIndexPosition, string &genomeID, int &numChromosomes,
                                   int &version, long &nviPosition, long &nviLength);

int flipSign(const char *genome, double *x, int n, char *chr, int binsize);

static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [-o observed][-t tol][-I max_iterations][-n normalization][-T threads][-v verbose] <hicfile> <outfile> <resolution> \n", argv0);
  fprintf(stderr, "  <hicfile>: hic file\n");
  fprintf(stderr, "  <outfile>: eigenvectors output file\n");
  fprintf(stderr, "  <resolution>: resolution in bp\n");
}

int main(int argc, char *argv[]) {
        string norm("NONE");
        string unit("BP");
        ifstream fin;
        struct timespec t0,t1,st,en;

        double tol=1.0e-6;
        int maxiter=500;
        int threads = 1;
        int verb = 1;
	string ob("oe");
        int opt;

	clock_gettime(CLOCK_REALTIME,&st);

	while ((opt = getopt(argc, argv, "o:t:I:T:n:v:h")) != -1) {
        	switch (opt) {
                                case 'o':
                                        if (strcmp(optarg,"observed") == 0) ob = "observed";
					else {
						usage(argv[0]);
                                        	exit(EXIT_FAILURE);
					}
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

	if (argc - optind < 3) {
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

	char *out_name = argv[optind++];
	FILE *fout = fopen(out_name,"w");
	if (fout==NULL) {
		fprintf(stderr, "Error! File %s cannot be opened for writing\n", argv[optind-1]);
		exit(EXIT_FAILURE);
	  }
	else if (verb) printf("Writing eigenvectors to %s\n\n", argv[optind-1]);
	int binsize = atoi(argv[optind++]);
	fprintf(fout,"track type=wiggle_0\n");

	int *ii;
	int *jj;
	double *xx;
	long p,q;
	long m;
	int k;
	int i;

	string hg19("hg19");
	string hg38("hg38");

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
    fin.close();

	if (genomeID != hg19 && genomeID != hg38) if (verb) {
		cout << genomeID;
		cout << " is not currently supported; no sign flip will be attempted!" << endl << endl;
	}

	vector<contactRecord> records;
	int iter0 = 0;

//	loop accross all chromosomes
	for (map<string,chromosome>::iterator itr =  chromosomeMap.begin(); itr !=  chromosomeMap.end(); ++itr) {
		string chrom = itr->second.name;
		char *chr = const_cast<char*> (chrom.c_str());
		if (strcmp(chr,"Y")==0 || strcmp(chr,"M")==0 || strcmp(chr,"MT")==0) continue;
		if (strcmp(chr,"chrY")==0 || strcmp(chr,"chrM")==0 || strcmp(chr,"chrMT")==0) continue;
		if (strcmp(chr,"ALL")==0 || strcmp(chr,"All")==0 || strcmp(chr,"all")==0) continue;
		clock_gettime(CLOCK_REALTIME,&t0);
        	records = straw(ob, norm, fname, chrom, chrom, unit, binsize);
		if (records.size() == 0) exit(EXIT_FAILURE);
        	m = records.size();
		ii = (int *) malloc(m*sizeof(int));
		jj = (int *) malloc(m*sizeof(int));
		xx = (double *) malloc(m*sizeof(double));
		p=0;
                for (q=0; q<m; q++) {
                        if (!isnan((double) records[q].counts)) {
                		ii[p] = records[q].binX/binsize;
                        	jj[p] = records[q].binY/binsize;
                        	xx[p] = (double) records[q].counts;
				p++;
			}
                }
        	records.clear();
        	records.shrink_to_fit();
		clock_gettime(CLOCK_REALTIME,&t1);
		if (verb > 1) printf("chromosome %s:  took %10.3f seconds for %ld records\n",chr,((double) (t1.tv_sec - t0.tv_sec)) + ((double) (t1.tv_nsec - t0.tv_nsec))/1e9,m);
		m = p;
		for (p=0; p<m;p++) if (jj[p] == ii[p]) xx[p] *= 0.5;

		int N = (int) ceil(itr->second.length/((double) binsize));

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
		if (verb > 1) {
			printf("total %d iterations\n",iter);
                	printf("iterations took %10.3f seconds\n",((double) (t1.tv_sec - t0.tv_sec)) + ((double) (t1.tv_nsec - t0.tv_nsec))/1e9);
			printf("lam1 = %g; lam2 = %g; lam1/lam2 = %g; er = %g; error in EV = %g\n",l1,l2,l1/l2,er,er/(l1-l2));
			printf("                           -------------------- \n");
		}

		fprintf(fout,"fixedStep ");
                char *chr1 = (char *) malloc((strlen(chr)+4)*sizeof(char));
                if (!strstr(chr,"chr")) strcpy(chr1,"chr");
                else strcpy(chr1,"");
                strcat(chr1,chr);
                if (strcmp(chr1,"chrMT") == 0)  strcpy(chr1,"chrM");
                fprintf(fout,"chrom=%s ",chr1);
                fprintf(fout,"start=1 step=%d span=%d\n",binsize,binsize);
		char *genome1 = const_cast<char*> (genomeID.c_str());
		if (100000 % binsize == 0) int junk = flipSign(genome1,ev,N,chr1,binsize);
                for (int j=0;j< N;j++) {
			if (j == N-1) fprintf(fout,"fixedStep chrom=%s start=%d step=%ld span=%ld\n",chr1,j*binsize+1,itr->second.length % binsize, itr->second.length % binsize );
			if (!isnan(ev[j])) fprintf(fout,"%f\n",ev[j]);
			else fprintf(fout,"%s\n","0");
//			else fprintf(fout,"%s\n","NaN");
		}

		free(ii);
		free(jj);
		free(xx);
	}
	fclose(fout);
	clock_gettime(CLOCK_REALTIME,&en);
	if (verb) printf("\n**************    all together took %10.3f seconds\n",((double) (en.tv_sec - st.tv_sec)) + ((double) (en.tv_nsec - st.tv_nsec))/1e9);
	return(0);
}

