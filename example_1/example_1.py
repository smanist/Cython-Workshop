import time
from prime import prime_python
from prime import prime_cython

def runTime(f, arg=1000, N=10):
    t0 = time.time()
    for i in xrange(N):
        f(arg)
    t = (time.time()-t0)/N
    print("Time elapsed: {0:5.4}ms".format(t*1000))

runTime(prime_python.prime1, 1000)
runTime(prime_cython.prime2, 1000)
runTime(prime_cython.prime3, 1000)
