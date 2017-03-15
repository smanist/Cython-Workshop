#include <math.h>
#include <stdio.h>
#include "pc_utils.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

void cAddMat(double *M,
             double *coef,
             double *m0,
             double *m1,
             unsigned int block,
             unsigned int band,
             unsigned int size) {
  unsigned int i, idx, j, k;
  unsigned int N = (size+1)*2;

/* #pragma omp parallel num_threads(4) */
  {
/* #pragma omp for */
    for (i = 0; i < size; i++) {
      idx = i*block;
      for (j = 0; j < band; j++) {
        for (k = 0; k < band; k++)
          M[(idx+j)*N+idx+k] += coef[i]*m0[band*j+k] + coef[i+1]*m1[band*j+k];
      }
    }
  }
}  // cAddMat

void cAddVec(double *V,
             double *coef,
             double *v0,
             double *v1,
             unsigned int block,
             unsigned int band,
             unsigned int size) {
  unsigned int i, j;

  for (i = 0; i < size; i++) {
    for (j = 0; j < band; j++) {
      V[i*block+j] += coef[i]*v0[j] + coef[i+1]*v1[j];
    }
  }
}  // cAddVec

void cDecomp(double *M,
             unsigned int band,
             unsigned int size) {
  unsigned int i, idx, j, jdx, k;

  for (i = 0; i < size; i++) {
    idx = MIN(i+band, size);
    M[i*size+i] = sqrt(M[i*size+i]);
    for (j = i+1; j < idx; j++)
      M[j*size+i] /= M[i*size+i];
    for (j = i+1; j < idx; j++) {
      jdx = MIN(j+band, size);
      for (k = j; k < jdx; k++)
        M[k*size+j] -= M[j*size+i]*M[k*size+i];
    }
  }
}  // cDecomp

void cSolve(double *M,
            double *b,
            double *x,
            double *y,
            unsigned int band,
            unsigned int size) {
  unsigned int i, idx, j, jdx;
  double val;

  y[0] = b[0] / M[0];
  for (i = 1; i < size; i++) {
    idx = MAX(i-band, 0);
    val = 0;
    for (j = idx; j < i; j++)
      val += y[j]*M[i*size+j];
    y[i] = (b[i] - val) / M[i*size+i];
  }
  for (i = 0; i < size; i++) {
    jdx = size-1-i;
    idx = MIN(jdx+band, size);
    val = 0;
    for (j = jdx+1; j < idx; j++)
      val += x[j]*M[j*size+jdx];
    x[jdx] = (y[jdx] - val) / M[jdx*size+jdx];
  }
}  // cSolve
