import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

cdef extern from "pc_utils.h":
    void cAddMat(double *, double *, double *, double *,
                 unsigned int, unsigned int, unsigned int)
    void cAddVec(double *, double *, double *, double *,
                 unsigned int, unsigned int, unsigned int)
    void cDecomp(double *, unsigned int, unsigned int);
    void cSolve(double *, double *, double *, double *,
                unsigned int, unsigned int);

def AddMat(np.ndarray[DTYPE_t, ndim=2] M,
           np.ndarray[DTYPE_t, ndim=1] coef,
           np.ndarray[DTYPE_t, ndim=2] m0,
           np.ndarray[DTYPE_t, ndim=2] m1,
           unsigned int block,
           unsigned int band):
    cdef unsigned int N = coef.shape[0] - 1
    cAddMat(<double *> M.data,
            <double *> coef.data,
            <double *> m0.data,
            <double *> m1.data,
            block, band, N)

def AddVec(np.ndarray[DTYPE_t, ndim=1] V,
           np.ndarray[DTYPE_t, ndim=1] coef,
           np.ndarray[DTYPE_t, ndim=1] v0,
           np.ndarray[DTYPE_t, ndim=1] v1,
           unsigned int block,
           unsigned int band):
    cdef unsigned int N = coef.shape[0] - 1
    cAddVec(<double *> V.data,
            <double *> coef.data,
            <double *> v0.data,
            <double *> v1.data,
            block, band, N)

def Decomp(np.ndarray[DTYPE_t, ndim=2] M,
           unsigned int q):
    cdef unsigned int N = M.shape[0]
    cDecomp(<double *> M.data, q, N)

def Solve(np.ndarray[DTYPE_t, ndim=2] M,
          np.ndarray[DTYPE_t, ndim=1] b,
          np.ndarray[DTYPE_t, ndim=1] x,
          unsigned int q):
    cdef int N = M.shape[0]
    cdef np.ndarray[DTYPE_t, ndim=1] y = np.zeros((N,))
    
    cSolve(<double *> M.data,
           <double *> b.data,
           <double *> x.data,
           <double *> y.data,
           q, N)

    return x
