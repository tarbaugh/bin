#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

// function to swap elements
void swapID(int *a, int *b) {
  int t = *a;
  *a = *b;
  *b = t;
}

void swapPos(double **a, double **b) {
  double *t;
  t = (double*) calloc(3, sizeof(double));
  t = *a;
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

int main(int argc, char *argv[]) {
    time_t t;
    srand((unsigned) time(&t));

    double ***x;
    int i, j, n, N, nconf, **types;
    nconf = 2;
    N = 5;

    x = (double***) calloc(2*nconf, sizeof(double**));  
    for (n=0; n<nconf; n++) {
        x[n] = (double**) calloc(N, sizeof(double*));
        for (i=0; i<N; i++) {
            x[n][i] = (double*) calloc(3, sizeof(double));
            for (j=0; j<3; j++) {
                x[n][i][j] = (float)rand()/(float)(2147483647/20);
            }
        }
    }
  
  types = (int**) calloc(2*nconf, sizeof(int*));  
  for (n=0; n<nconf; n++) {
    types[n] = (int*) calloc(N, sizeof(int));
    for (j=0; j<N; j++) {
        types[n][j] = rand() % (N + 1 - 1) + 1;
    }
  }
for (i=0; i<N; i++) {
            printf("%d\t%lf\t%lf\t%lf\n", types[0][i], x[0][i][0], x[0][i][1], x[0][i][2]);
        }
printf("\n\n");
  quickSort(types[0], x[0], 0, N-1);
    // for (n=0; n<nconf; n++) {
for (i=0; i<N; i++) {
    printf("%d\t%lf\t%lf\t%lf\n", types[0][i], x[0][i][0], x[0][i][1], x[0][i][2]);
}
    // }
}