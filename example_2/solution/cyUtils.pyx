import numpy as np
cimport numpy as np
cimport cython

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)
def AddMat(np.ndarray[DTYPE_t, ndim=2] M,
           np.ndarray[DTYPE_t, ndim=1] coef,
           np.ndarray[DTYPE_t, ndim=2] m0,
           np.ndarray[DTYPE_t, ndim=2] m1,
           unsigned int block,
           unsigned int band):
    cdef unsigned int N = coef.shape[0] - 1
    cdef unsigned int i, idx, j, k

    for i in range(N):
        idx = i*block
        for j in range(band):
            for k in range(band):
                M[idx+j, idx+k] = M[idx+j, idx+k] + coef[i]*m0[j,k] + coef[i+1]*m1[j,k]

@cython.boundscheck(False)
@cython.wraparound(False)
def AddVec(np.ndarray[DTYPE_t, ndim=1] V,
           np.ndarray[DTYPE_t, ndim=1] coef,
           np.ndarray[DTYPE_t, ndim=1] v0,
           np.ndarray[DTYPE_t, ndim=1] v1,
           unsigned int block,
           unsigned int band):
    cdef unsigned int N = coef.shape[0] - 1
    cdef unsigned int i, idx, j

    for i in range(N):
        idx = i*block
        for j in range(band):
            V[idx+j] = V[idx+j] + coef[i]*v0[j] + coef[i+1]*v1[j]

@cython.boundscheck(False)
@cython.wraparound(False)
def Decomp(np.ndarray[DTYPE_t, ndim=2] M,
           unsigned int q):
    cdef unsigned int N = M.shape[0]
    cdef unsigned int i, idx, j, jdx, k
    
    for i in range(N):
        idx = min(i+q, N)
        M[i,i] = np.sqrt(M[i,i])
        for j in range(i+1, idx):
            M[j,i] = M[j,i] / M[i,i]
        for j in range(i+1, idx):
            jdx = min(j+q, N)
            for k in range(j, jdx):
                M[k,j] = M[k,j] - M[j,i]*M[k,i]

@cython.boundscheck(False)
@cython.wraparound(False)
def Solve(np.ndarray[DTYPE_t, ndim=2] M,
          np.ndarray[DTYPE_t, ndim=1] b,
          np.ndarray[DTYPE_t, ndim=1] x,
          unsigned int q):
    cdef unsigned int N = M.shape[0]
    cdef unsigned int i, idx, j
    cdef np.ndarray[DTYPE_t, ndim=1] y = np.zeros((N,))
    cdef DTYPE_t val
    
    y[0] = b[0] / M[0,0]
    for i in range(1, N):
        idx = max(i-q, 0)
        val = 0
        for j in range(idx, i):
            val = val + y[j]*M[i,j]
        y[i] = (b[i] - val) / M[i,i]
    for i in range(N-1, -1, -1):
        idx = min(i+q, N)
        val = 0
        for j in range(i+1, idx):
            val = val + x[j]*M[j,i]
        x[i] = (y[i] - val) / M[i,i]
    return x
