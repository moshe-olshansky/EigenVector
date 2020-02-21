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

int bestEigen(long m,int *i,int *j,double *x,int *N,double *l1,double *l2,double *a1,double *a2,double *ev,double *er,double tol,int maxiter,int threads);

static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [-t tol][-I max_iterations][-n normalization][-T threads][-v verbose] <hicfile> <chromosome> <outfile> <resolution> \n", argv0);
  fprintf(stderr, "  <hicfile>: hic file\n");
  fprintf(stderr, "  <chromosome>: chromosomee\n");
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
	int N;
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

	string str;
	getline(fin, str, '\0' );
//	cout << str << "\n";
	fin.read((char*)&version, sizeof(int));
	if (version < 6) {
		cerr << "Version " << version << " no longer supported" << endl;
		 exit(1);
	}
//	cout << "version = " << version << "\n";
	long master;
	fin.read((char*)&master, sizeof(long));
	string genome;
	getline(fin, genome, '\0' );
	int nattributes;
	fin.read((char*)&nattributes, sizeof(int));
//	cout << master << " " << genome << " " << nattributes << "\n";

// reading and ignoring attribute-value dictionary
	for (int i=0; i<nattributes; i++) {
		string key, value;
		getline(fin, key, '\0');
		getline(fin, value, '\0');
//		cout << key << " " << value << "\n";
//		cout << key << "\n";
	}
	int nChrs;
	fin.read((char*)&nChrs, sizeof(int));
	bool found=false;
// chromosome map for finding matrix
	for (int i=0; i<nChrs; i++) {
		string name;
		int length;
		getline(fin, name, '\0');
		fin.read((char*)&length, sizeof(int));
		if (name == chrom) {
			N = (int) ceil(length/((double) binsize));
			found = true;
			break;
		}
    	}
	if (!found) {
		cout << "chromosome " << chrom << " is not found" << endl;
		exit(1);
	}
	fin.close();
	ftime(&t0);
	vector<contactRecord> records = straw(norm, fname, chrom, chrom, unit, binsize);
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
	ftime(&t1);
	if (verb) printf("took %ld seconds for %ld records\n",(long) (t1.time - t0.time),nonZer);
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

	ftime(&t0);
	int iter;
	iter = bestEigen(m,ii,jj,xx,&N,&l1,&l2,a1,a2,ev,&er,tol,maxiter,threads);
	ftime(&t1);

	if (verb) {
		printf("total %d iterations\n",iter);
		printf("iterations took %15.10lf seconds\n",((double) (t1.time - t0.time)) + 0.001*(t1.millitm - t0.millitm));
		printf("lam1 = %20.6f; lam2 = %20.6f; er = %g; error in EV = %g\n",l1,l2,er,er/(l1-l2));
	}
	for (int j=0;j<N;j++) {
		if (!isnan(ev[j])) fprintf(fout,"%20.10f\n",ev[j]);
		else fprintf(fout,"%20s\n","NaN");
	}
	fclose(fout);
	return(iter);
}

