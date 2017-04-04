import numpy as np
from scipy import linalg

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
    pass

def Solve(M, b, x, q):
    np.copyto(x, linalg.solve(M, b, sym_pos=True))
