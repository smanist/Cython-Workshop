import numpy as np

def AddMat(M, coef, m0, m1, block, band):
    for i in xrange(len(coef)-1):
        idx = i*block
        jdx = idx+band
        M[idx:jdx, idx:jdx] += coef[i]*m0 + coef[i+1]*m1

def AddVec(V, coef, v0, v1, block, band):
    for i in xrange(len(coef)-1):
        idx = i*block
        V[idx:idx+band] += coef[i]*v0 + coef[i+1]*v1

def Decomp(M, q):
    N = len(M)
    for i in xrange(N):
        idx = min(i+q, N)
        M[i,i] = np.sqrt(M[i,i])
        M[i+1:idx,i] = M[i+1:idx,i] / M[i,i]
        for j in xrange(i+1, idx):
            jdx = min(j+q, N)
            M[j:jdx,j] = M[j:jdx,j] - M[j,i]*M[j:jdx,i]

def Solve(M, b, q):
    N = len(M)
    x = np.zeros((N,))
    y = np.zeros((N,))
    y[0] = b[0] / M[0,0]
    for i in xrange(1, N):
        idx = max(i-q, 0)
        y[i] = (b[i] - np.dot(y[idx:i], M[i,idx:i])) / M[i,i]
    for i in xrange(N-1, -1, -1):
        idx = min(i+q, N)
        x[i] = (y[i] - np.dot(x[i+1:idx], M[i+1:idx,i])) / M[i,i]
    return x
