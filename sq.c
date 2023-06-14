#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define ATOMS 1
#define NDIM 3
#define anint(x) ((x >= .5) ? (1.0) : (x <= -.5) ? (-1.0) : (0.0))
#define QMAX 200
#define VECMAX 300
#define NA 602.205
#define MASS 18.016

typedef struct {double r, i;} complex;

int N;

void read_qvectors(int qvector[QMAX][VECMAX][NDIM], int nvec[VECMAX]);

void main(int argc, char **argv)
{
  double **x, rho, L, norm, sf[QMAX];
  complex rhoq;
  float *xx;
  int i, j, n, k, t, vec, nconf, nlin, nlog, start, end, repeat, idum, vdum;
  int q, qi, qj, qk;
  double qx, qy, qz, qmin, qdotr, ddum, l1, l2;
  int qvector[QMAX][VECMAX][NDIM], nvec[VECMAX];
  char fname[100], outfile[100], command[100], ignore[1024];;
  FILE *fp;

  if (argc != 4) {
    printf("\nUsage: msd <temp> <start> <end> \n\n");
    exit(0);
  }

  start = atoi(argv[2]);
  end = atoi(argv[3]);

 /* read number of molecules and box size from coordinate file*/
  sprintf(fname, "%s.%02d.xyz", argv[1], start);
  if ((fp = fopen(fname, "r"))==NULL){
    printf("\nError opening file %s;  Program aborted!\n\n", fname);
    exit(1);
  }
  fgets(ignore, sizeof(ignore), fp);
  fgets(ignore, sizeof(ignore), fp);
  fgets(ignore, sizeof(ignore), fp);
  fscanf(fp, "%d", &N);
  fgets(ignore, sizeof(ignore), fp);
  fgets(ignore, sizeof(ignore), fp);
  fscanf(fp, "%lf %lf", &l1, &l2);
  fclose(fp);
  L = l2 - l1;
  printf("Box Size: %lf\nNumber of Atoms: %d\n", L, N);

  /* get starting and ending coordinate file */
  sprintf(command, "ls %s.*.xyz | cut -d. -f2 | sort -n | head -1", argv[1]);
  if ((fp = popen(command, "r")) == NULL) {
    printf("Could not get repeat frequency from log file\n");
    exit(0);
  }
  fscanf(fp, "%d", &start);
  pclose(fp);

  sprintf(command, "ls %s.*.xyz | cut -d. -f2 | sort -n | tail -1", argv[1]);
  if ((fp = popen(command, "r")) == NULL) {
    printf("Could not get repeat frequency from log file\n");
    exit(0);
  }
  fscanf(fp, "%d", &end);
  pclose(fp);

 /* read number of molecules from xyz file*/
  // sprintf(fname, "%s.%03d.xyz", argv[1], start);
  // if ((fp = fopen(fname, "r"))==NULL){
  //   printf("\nError opening file %s;  Program aborted!\n\n", fname);
  //   exit(1);
  // }
  // fscanf(fp, "%d", &N);
  N /= ATOMS;
  // fclose(fp);

  /* number of configurations per datafile */
  // sprintf(command, "grep ^run %s | head -1 | awk '{print $2}'", argv[2]);
  // if ((fp = popen(command, "r")) == NULL) {
  //   printf("Could not get repeat frequency from log file\n");
  //   exit(0);
  // }
  // fscanf(fp, "%d", &repeat);
  // pclose(fp);
  // nlog = ceil(log(repeat)/log(2.0))+1;

  // /* double-check nlog against logdata file */
  // sprintf(command, "wc %s.001.xyz  | awk '{print $1/%d}'", argv[1], N*ATOMS+2);
  // if ((fp = popen(command, "r")) == NULL) {
  //   printf("Could not get repeat frequency from log file\n");
  //   exit(0);
  // }
  // fscanf(fp, "%d", &idum);
  // pclose(fp);
  // if ((idum-1) != nlog) {
  //   printf("Inconsistency in value of nlog from xyz and lammps-log files\n");
  //   exit(0);
  // }

  // /* number of data files */
  // sprintf(command, "ls %s.*.xyz | wc | awk '{print $1}'", argv[1]);
  // if ((fp = popen(command, "r")) == NULL) {
  //   printf("Could not get repeat frequency from log file\n");
  //   exit(0);
  // }
  // fscanf(fp, "%d", &nlin);
  // pclose(fp);

  // nconf = nlin*nlog;

  nlin = end - start;
  nlog = ceil(log(32768)/log(2.0))+1;
  nconf = nlin*nlog;

  printf("N: %d\tnconf: %d\tL: %lf\tRepeat: %d\n", N, nconf, L, repeat);

  /* Calculate qmin*/
  qmin = 2.0*M_PI/L;
  printf("qmin = %lf\n", qmin);
  
  /* read qvectors */
  printf("Reading q-vectors to calculate...\n");
  fflush(NULL);
  read_qvectors(qvector, nvec);

  /* allocate needed memory */
  x = (double**) calloc(N, sizeof(double*));
  for (i=0; i<N; i++) 
    x[i] = (double*) calloc(NDIM, sizeof(double));
  
  xx = (float*) calloc(N*ATOMS*NDIM, sizeof(float));
  
  for (q=0; q<QMAX; q++) 
    sf[q]=0;
  
  printf("\nCalculating...    ");
  fflush(NULL);
  /* loop over all files */
  for(n = 0; n < nlin; n++) {
    
    printf("\b\b\b\b%3.0f%%", 100.0*(double)n/(double)nlin);
    fflush(NULL);

    sprintf(fname, "%s.%02d.xyz.starr", argv[1], n+1);
    if ((fp = fopen(fname, "r"))==NULL){
      printf("\nError opening file %s;  Program aborted!\n\n", fname);
      exit(1);
    }
    
    fscanf(fp, "%d", &idum);
    //printf("%d", idum);
    fgets(command, 100, fp);
    fgets(command, 100, fp);

    for (i=0; i<N; i++) {
      fscanf(fp, "%d %lf %lf %lf", &vdum, &x[i][0], &x[i][1], &x[i][2]);
      // fscanf(fp, "%d %lf %lf %lf", &idum, &ddum, &ddum, &ddum);
      // fscanf(fp, "%d %lf %lf %lf", &idum, &ddum, &ddum, &ddum);
    }
  
    fclose(fp);
    
    for (q=1; q<QMAX; q++) {
      for(vec=0; vec<nvec[q]; vec++) {
	rhoq.r = rhoq.i = 0.0;
	
	qx = qmin*qvector[q][vec][0];
	qy = qmin*qvector[q][vec][1];
	qz = qmin*qvector[q][vec][2];
	
	for (i=0; i<N; i++) {
	  qdotr = x[i][0]*qx + x[i][1]*qy + x[i][2]*qz; 
	  rhoq.r += cos(qdotr);
	  rhoq.i += sin(qdotr);
	}
	sf[q] += rhoq.r*rhoq.r + rhoq.i*rhoq.i;
      }
      
    }
  }
  
  norm = 1.0/((double)nlin);
  norm /= (double)N;
  
  sprintf(outfile, "sq.dat");
  fp = fopen(outfile, "w");
  
  for (q=1; q<QMAX; q++) {
    fprintf(fp, "%lf\t%lf\n", (double)(q+1)*0.5*qmin, 
	    sf[q]*norm/(double)nvec[q]);
  }
  fclose(fp);
  
  printf("\nStructure factor in file sq.dat\n");
  
}

void read_qvectors(int qvector[QMAX][VECMAX][NDIM], int nvec[QMAX])
{
  FILE *fp;
  char fname[100];
  int q;
  
  for (q=1; q<QMAX; q++) {

    sprintf(fname, "/zfshomes/fstarr/qvector/qvector.%03d", q+1);
    if ((fp = fopen(fname, "r"))==NULL){
      printf("\nError opening file %s;  Program aborted!\n\n", fname);
      exit(1);
    }
    
    nvec[q] = 0;
    fscanf(fp, "%d %d %d", qvector[q][nvec[q]], qvector[q][nvec[q]]+1, 
	   qvector[q][nvec[q]]+2);
    while(!feof(fp)) {
      nvec[q]++;
      fscanf(fp, "%d %d %d", qvector[q][nvec[q]], qvector[q][nvec[q]]+1, 
	     qvector[q][nvec[q]]+2);
    }
    
    fclose(fp);
  }
}
