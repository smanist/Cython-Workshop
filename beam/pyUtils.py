import numpy as np

def AddMat(M, coef, m0, m1, block, band):
    N = len(coef)-1
    for i in xrange(N):
        idx = i*block
        for j in xrange(band):
            for k in xrange(band):
                M[idx+j, idx+k] += coef[i]*m0[j,k] + coef[i+1]*m1[j,k]

def AddVec(V, coef, v0, v1, block, band):
    N = len(coef)-1
    for i in xrange(N):
        idx = i*block
        for j in xrange(band):
            V[idx+j] += coef[i]*v0[j] + coef[i+1]*v1[j]

def Decomp(M, q):
    N = len(M)
    for i in xrange(N):
        idx = min(i+q, N)
        M[i,i] = np.sqrt(M[i,i])
        for j in xrange(i+1, idx):
            M[j,i] = M[j,i] / M[i,i]
        for j in xrange(i+1, idx):
            jdx = min(j+q, N)
            for k in xrange(j, jdx):
                M[k,j] = M[k,j] - M[j,i]*M[k,i]

def Solve(M, b, x, q):
    N = len(M)
    y = np.zeros((N,))
    y[0] = b[0] / M[0,0]
    for i in xrange(1, N):
        idx = max(i-q, 0)
        val = 0
        for j in range(idx, i):
            val = val + y[j]*M[i,j]
        y[i] = (b[i] - val) / M[i,i]
    for i in xrange(N-1, -1, -1):
        idx = min(i+q, N)
        val = 0
        for j in range(i+1, idx):
            val = val + x[j]*M[j,i]
        x[i] = (y[i] - val) / M[i,i]
    return x
