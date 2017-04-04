void cAddMat(double *M,
             double *coef,
             double *m0,
             double *m1,
             unsigned int block,
             unsigned int band,
             unsigned int size);

void cAddVec(double *V,
             double *coef,
             double *v0,
             double *v1,
             unsigned int block,
             unsigned int band,
             unsigned int size);

void cDecomp(double *M,
             unsigned int band,
             unsigned int size);

void cSolve(double *M,
            double *b,
            double *x,
            double *y,
            unsigned int band,
            unsigned int size);
