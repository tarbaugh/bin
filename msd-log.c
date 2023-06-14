#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NDIM 3
#define ATOMS 1
#define magic 6755399441055744.0
#define anint(x) ((x >= 0.5) ? (1.0) : (x <= -0.5) ? (-1.0) : (0.0))
#define DT 0.001

// function to swap elements
void swapID(int *a, int *b) {
  int t = *a;
  *a = *b;
  *b = t;
}

void swapPos(double **a, double **b) {
  double *t = *a;
  *a = *b;
  *b = t;
}

// function to find the partition position
int partition(int *array, double **pos, int low, int high) {
  
  // select the rightmost element as pivot
  int pivot = array[high];
  
  // pointer for greater element
  int i = (low - 1);

  // traverse each element of the array
  // compare them with the pivot
  for (int j = low; j < high; j++) {
    if (array[j] <= pivot) {
        
      // if element smaller than pivot is found
      // swap it with the greater element pointed by i
      i++;
      
      // swap element at i with element at j
      swapID(&array[i], &array[j]);
      swapPos(&pos[i], &pos[j]);
    }
  }

  // swap the pivot element with the greater element at i
  swapID(&array[i + 1], &array[high]);
  swapPos(&pos[i + 1], &pos[high]);
  
  // return the partition point
  return (i + 1);
}

void quickSort(int *array, double **pos, int low, int high) {
  if (low < high) {
    
    // find the pivot element such that
    // elements smaller than pivot are on left of pivot
    // elements greater than pivot are on right of pivot
    int pi = partition(array, pos, low, high);
    
    // recursive call on the left of pivot
    quickSort(array, pos, low, pi - 1);
    
    // recursive call on the right of pivot
    quickSort(array, pos, pi + 1, high);
  }
}

int main(int argc, char *argv[])
{
  double ***x, L, l1, l2, ***dr, rtot[NDIM], ddum;
  int i, k, n, nconf, nskip, nlog, ncalc, nlin, start, end, vdum, idum, **types, **ids;
  int N;
  int t, tt, *time;
  int *countlin, **countlog;
  double **msd_log, *msd_lin;
  char fname[100], command[100], ignore[2048];
  FILE *fp;

  if (argc != 4) {
    printf("\nUsage: msd <temp> <start[inclusive]> <end[inclusive]> \n\n");
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
  fscanf(fp,"%*s %*s");
  fscanf(fp,"%*d");
  fscanf(fp,"%*s %*s %*s %*s");
  fscanf(fp, "%d", &N);
  fscanf(fp,"%*s %*s %*s %*s %*s %*s");
  fscanf(fp, "%lf %lf", &l1, &l2);
  fclose(fp);
  L = l2 - l1;
  printf("Box Size: %lf\nNumber of Atoms: %d\n", L, N);

  nlin = end - start + 1;
  nlog = ceil(log(32768)/log(2.0))+1;
  nconf = nlin*nlog;
  
  printf("nlog: %d", nlog);
  printf("\nnconf: %d", nconf);


  x = (double***) calloc(2*nconf, sizeof(double**));  
  for (n=0; n<nconf; n++) {
    x[n] = (double**) calloc(N, sizeof(double*));
    for (i=0; i<N; i++) 
      x[n][i] = (double*) calloc(NDIM, sizeof(double));
  }
  
  types = (int**) calloc(2*nconf, sizeof(int*));  
  for (n=0; n<nconf; n++) {
    types[n] = (int*) calloc(N, sizeof(int));
  }

  ids = (int**) calloc(2*nconf, sizeof(int*));  
  for (n=0; n<nconf; n++) {
    ids[n] = (int*) calloc(N, sizeof(int));
  }

  time = (int*) calloc(2*nconf, sizeof(int));
  countlin = (int*) calloc(2*nlin, sizeof(int));
  countlog = (int**) calloc(2*nlog, sizeof(int*));
  for (i=0; i<nlog; i++)
    countlog[i] = (int*) calloc(2*nlog, sizeof(int));  
  msd_log = (double**) calloc(2*nlog, sizeof(double*));
  msd_lin = (double*) calloc(2*nlin, sizeof(double));
  for (i=0; i<nlog; i++) {
    msd_log[i] = (double*) calloc(2*nlog, sizeof(double));
  }
  dr = (double***) calloc(2*nconf, sizeof(double**));
  for (n=0; n<nconf; n++) {
    dr[n] = (double**) calloc(N, sizeof(double*));
    for (i=0; i<N; i++) 
      dr[n][i] = (double*) calloc(NDIM, sizeof(double));
    
  }
  
  printf("\nReading data files...     ");
  
  /* loop over all files */
  for(n = 0; n < nlin; n++){

    printf("\b\b\b\b%3.0f%%", 100.0*(double)n/(double)nlin);
    fflush(NULL);

    sprintf(fname, "%s.%02d.xyz", argv[1], n+start);
    if ((fp = fopen(fname, "r"))==NULL){
      printf("\nError opening file %s;  Program aborted!\n\n", fname);
      exit(1);
    }
        
    /* skip 0th config */
    if (n == 0) {
      fscanf(fp,"%*s %*s");
      fscanf(fp,"%*d");
      fscanf(fp,"%*s %*s %*s %*s");
      fscanf(fp,"%*d");
      fscanf(fp,"%*s %*s %*s %*s %*s %*s");
      fscanf(fp,"%*lf %*lf");
      fscanf(fp,"%*lf %*lf");
      fscanf(fp,"%*lf %*lf");
      fscanf(fp,"%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");

      for (i=0; i<N; i++) {
        fscanf(fp, "%*d %*d %*lf %*lf %*lf %*lf %*lf %*lf");
      }
    }
  
    for (k=0; k<nlog; k++ ) {
      fscanf(fp,"%*s %*s");
      fscanf(fp, "%d", &time[n*nlog+k]);
      fscanf(fp,"%*s %*s %*s %*s");
      fscanf(fp,"%*d");
      fscanf(fp,"%*s %*s %*s %*s %*s %*s");
      fscanf(fp,"%*lf %*lf");
      fscanf(fp,"%*lf %*lf");
      fscanf(fp,"%*lf %*lf");
      fscanf(fp,"%s %*s %*s %*s %*s %*s %*s %*s %*s %*s", ignore);

      time[n*nlog+k] += n*pow(2,nlog-1);
      for (i=0; i<N; i++) {
        fscanf(fp, "%d %d %lf %lf %lf %*lf %*lf %*lf", &ids[n*nlog+k][i], &types[n*nlog+k][i], &x[n*nlog+k][i][0], &x[n*nlog+k][i][1], &x[n*nlog+k][i][2]);
      }
    }
    fclose(fp);
  }

  for(n = 0; n < nlin; n++) {
    for (k=0; k<nlog; k++ ) {
      quickSort(ids[n*nlog+k], x[n*nlog+k], 0, N-1);
      // for (i=0; i<N; i++) {
      //   x[n*nlog+k][i][0] += l1;
      //   x[n*nlog+k][i][1] += l1;
      //   x[n*nlog+k][i][2] += l1;
      // }
    }
  }
  
  for (t=0; t<nlin; t++) {
      countlin[t] = 0;
      msd_lin[t]  = 0;
  }
  for (t=0; t<nlog; t++) 
    for (tt=0; tt<nlog; tt++){
      countlog[t][tt] = 0;  
      msd_log[t][tt]  = 0;
    } 
  

  for (t=1; t<nconf; t++) {
    printf("%lf\n",x[t][0][2]);
  }
  printf("\nCaculating incremental diplacements....     ");
  for (i=0; i<N; i++) {
    printf("\b\b\b\b%3.0f%%", 100.0*(double)i/(double)N);
    fflush(NULL);
    for (t=1; t<nconf; t++) 
        for (k=0; k<NDIM; k++) {
	  dr[t][i][k] = x[t][i][k] - x[t-1][i][k];
        /*dr[t][i][k] -= L*(dr[t][i][k]/L+magic-magic);*/
	  dr[t][i][k] -= L*anint(dr[t][i][k]/L);
    }
  }
  

  printf("\nCaculating MSD for type 0...     ");
  for (i=0; i<N; i++) {
    if (types[0][i] == 0) continue;
    printf("\b\b\b\b%3.0f%%", 100.0*(double)i/(double)N);
    fflush(NULL);
    
  /*printf("\nCalculating log-spaced msd...    ");*/
  /*for (t=0; t<nconf; t++) {  USE THIS LINE TO CALCULATE ALL CORR */
    for (t=0; t<nconf; t+=nlog) {
      rtot[0] = rtot[1] = rtot[2] = 0.0;
      start = (t/nlog)*nlog;
      
      for (tt=t+1; tt<start+nlog; tt++) {
        for (k=0; k<NDIM; k++) {
          rtot[k] += dr[tt][i][k];
          msd_log[t-start][tt-start] += rtot[k]*rtot[k];
        }
        countlog[t-start][tt-start]++;
      }
    }    
    
    
  /*printf("\nCalculating linear-spaced msd...    ");*/
    for (t=0; t<nconf-nlog; t++) {
      rtot[0] = rtot[1] = rtot[2] = 0.0;
	  
      for (tt=t+1; tt<nconf; tt++) {
        if (((tt-t)%nlog)==0) {
          countlin[(tt-t)/nlog]++;
          for (k=0; k<NDIM; k++) {
            rtot[k] += dr[tt][i][k]; 
            msd_lin[(tt-t)/nlog] += rtot[k]*rtot[k];
          }
        }
        else
          for (k=0; k<NDIM; k++) {
	    rtot[k] += dr[tt][i][k]; 
          }
	    
      }
    }
  }
  
  fp = fopen("msd0.dat", "w");
  fprintf(fp, "#t\tr^2\n"); 
      
  /* first log correlations */
  for (tt=0; tt<nlog; tt++) 
    for (t=tt-1; t>=0; t--)
      if (msd_log[t][tt] != 0.0)
        fprintf(fp, "%lf\t%le\n", 
	        DT*(time[tt]-time[t]),
	        msd_log[t][tt]/(double)(countlog[t][tt]));
  
  /*  now linear correlations */
  for (t=0; t<nlin; t++) 
    if (msd_lin[t] != 0.0)
      fprintf(fp, "%lf\t%le\n", 
              DT*(time[t*nlog]-time[0]),	
              msd_lin[t]/(double)(countlin[t]));
  fclose(fp);
    
  printf("\nData in msd0.dat\n");

  for (t=0; t<nlin; t++) {
      countlin[t]  = 0;
      msd_lin[t]   = 0;
  }
  for (t=0; t<nlog; t++) 
    for (tt=0; tt<nlog; tt++){
      countlog[t][tt] = 0;  
      msd_log[t][tt]  = 0;
    } 
  

  printf("\nCaculating MSD for type 1...     ");
  for (i=0; i<N; i++) {
    if (types[0][i] == 1) continue;
    printf("\b\b\b\b%3.0f%%", 100.0*(double)i/(double)N);
    fflush(NULL);
    
  /*printf("\nCalculating log-spaced msd...    ");*/
  /*for (t=0; t<nconf; t++) {  USE THIS LINE TO CALCULATE ALL CORR */
    for (t=0; t<nconf; t+=nlog) {
      rtot[0] = rtot[1] = rtot[2] = 0.0;
      start = (t/nlog)*nlog;
      
      for (tt=t+1; tt<start+nlog; tt++) {
        for (k=0; k<NDIM; k++) {
          rtot[k] += dr[tt][i][k];
          msd_log[t-start][tt-start] += rtot[k]*rtot[k];
        }
        countlog[t-start][tt-start]++;
      }
    }    
    
    
  /*printf("\nCalculating linear-spaced msd...    ");*/
    for (t=0; t<nconf-nlog; t++) {
      rtot[0] = rtot[1] = rtot[2] = 0.0;
	  
      for (tt=t+1; tt<nconf; tt++) {
        if (((tt-t)%nlog)==0) {
          countlin[(tt-t)/nlog]++;
          for (k=0; k<NDIM; k++) {
            rtot[k] += dr[tt][i][k]; 
            msd_lin[(tt-t)/nlog] += rtot[k]*rtot[k];
          }
        }
        else
          for (k=0; k<NDIM; k++) {
	    rtot[k] += dr[tt][i][k]; 
          }
	    
      }
    }
  }
  
  fp = fopen("msd1.dat", "w");
  fprintf(fp, "#t\tr^2\n"); 
      
  /* first log correlations */
  for (tt=0; tt<nlog; tt++) 
    for (t=tt-1; t>=0; t--)
      if (msd_log[t][tt] != 0.0)
        fprintf(fp, "%lf\t%le\n", 
	        DT*(time[tt]-time[t]),
	        msd_log[t][tt]/(double)(countlog[t][tt]));
  
  /*  now linear correlations */
  for (t=0; t<nlin; t++) 
    if (msd_lin[t] != 0.0)
      fprintf(fp, "%lf\t%le\n", 
              DT*(time[t*nlog]-time[0]),	
              msd_lin[t]/(double)(countlin[t]));
  fclose(fp);
    
  printf("\nData in msd1.dat\n");

}