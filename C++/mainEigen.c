#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/timeb.h>
#include <time.h>

int bestEigen(long m,int *i,int *j,double *x,int *N,double *l1,double *l2,double *a1,double *a2,double *ev,double *er,double tol,int maxiter,int threads);

static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [-m records][-t tol][-I max_iterations][-s size][-T threads][-v verbose] <infile> <outfile> <resolution> \n", argv0);
  fprintf(stderr, "  <infile>: the \'triplets\' file\n");
  fprintf(stderr, "  <outfile>: eigenvector output  file\n");
  fprintf(stderr, "  <resolution>: resolution in bp\n");
}


int main(int argc, char *argv[]) {
	struct timeb t0,t1;

        double tol=1.0e-6;
        int maxiter=500;
        int threads = 1;
        int verb = 1;
        int s=-1;
        long m1=1e9;
        int opt;

        while ((opt = getopt(argc, argv, "m:t:I:T:s:v:h")) != -1) {
                switch (opt) {
                                case 't':
                                        tol = atof(optarg);
                                        break;
                                case 'I':
                                        maxiter=atoi(optarg);
                                        break;
                                case 'm':
                                        m1=atoi(optarg);
                                        break;
                                case 's':
                                        s=atoi(optarg);
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


  FILE *fin = fopen(argv[optind++],"r");
  if (fin==NULL) {
    fprintf(stderr, "Error! File %s cannot be opened for reading\n", argv[optind-1]);
    exit(EXIT_FAILURE);
  }
  else if (verb) printf("Reading from %s", argv[optind-1]);

  FILE *fout = fopen(argv[optind++],"w");
  if (fout==NULL) {
    fprintf(stderr, "Error! File %s cannot be opened for writing\n", argv[optind-1]);
    exit(EXIT_FAILURE);
  }
  else if (verb) printf(" and writing to %s\n", argv[optind-1]);
  
  int binsize = atoi(argv[optind++]);

  long m,m0,p;
  int *i,*j;
  double *x;

  i = (int *) malloc(m1*sizeof(int));
  if (i == NULL) {
    fprintf(stderr,"Failed to allocate %ld bytes. Exiting\n",m1*sizeof(int));
    exit(EXIT_FAILURE);
  }
  j = (int *) malloc(m1*sizeof(int));
  if (j == NULL) {
    fprintf(stderr,"Failed to allocate %ld bytes. Exiting\n",m1*sizeof(int));
    exit(EXIT_FAILURE);
  }
  x = (double *) malloc(m1*sizeof(double));
  if (x == NULL) {
    fprintf(stderr,"Failed to allocate %ld bytes. Exiting\n",m1*sizeof(double));
    exit(EXIT_FAILURE);
  }
  m = 0; 
  m0 = m1;
  ftime(&t0);
  while(fscanf(fin,"%d %d %lg",&i[m],&j[m],&x[m]) == 3) {
    i[m] /= binsize;	
    j[m] /= binsize;	
    if (isnan(x[m])) x[m] = 0;
    m++;
    if (m == m1) {
      m1 += m0;
      i = (int *) realloc(i,m1*sizeof(int));
      if (i == NULL) {
        fprintf(stderr,"Failed to allocate additional %ld bytes. Exiting\n",m0*sizeof(int));
        exit(EXIT_FAILURE);
      }
      j = (int *) realloc(j,m1*sizeof(int));
      if (j == NULL) {
        fprintf(stderr,"Failed to allocate additional %ld bytes. Exiting\n",m0*sizeof(int));
        exit(EXIT_FAILURE);
      }
      x = (double *) realloc(x,m1*sizeof(double));
      if (x == NULL) {
        fprintf(stderr,"Failed to allocate additional %ld bytes. Exiting\n",m0*sizeof(double));
        exit(EXIT_FAILURE);
      }
    }
  }
  ftime(&t1);
  fclose(fin);
  if (verb) {
	  printf("finished reading\n");
  	  printf("took %ld seconds to read %ld records\n",(long) (t1.time - t0.time),m);
  }

  int k;
  if (s > 0) k = (int) ceil(((double) s)/binsize);
  else {
  	k = 0;
  	for (p=0; p<m;p++) if (j[p] > k) k=j[p];
  	k++;
  }

  for (p=0; p<m;p++) if (j[p] == i[p]) x[p] *= 0.5;

  double l1,l2,er;
  double *a1 = (double *) malloc(k*sizeof(double));
  double *a2 = (double *) malloc(k*sizeof(double));
  double *ev = (double *) malloc(k*sizeof(double));

  long init_seed = 1234563;
  srand48(init_seed);
  for (p=0;p<k;p++) {
  	a1[p] = drand48();
  	a2[p] = drand48();
  }

  ftime(&t0);
  int iter;
  iter = bestEigen(m,i,j,x,&k,&l1,&l2,a1,a2,ev,&er,tol,maxiter,threads);
  ftime(&t1);

   if (verb) {
          printf("total %d iterations\n",iter);
          printf("iterations took %15.10lf seconds\n",((double) (t1.time - t0.time)) + 0.001*(t1.millitm - t0.millitm));
          printf("lam1 = %20.6f; lam2 = %20.6f; er = %g; error in EV = %g\n",l1,l2,er,er/(l1-l2));
  }
  for (p=0;p<k;p++) {
          if (!isnan(ev[p])) fprintf(fout,"%20.10f\n",ev[p]);
          else fprintf(fout,"%20s\n","NaN");
  }
  fclose(fout);
  return(iter);
}

  
