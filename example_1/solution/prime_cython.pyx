import cython
cimport cython

def prime2(Range):
    """
    Find the first few primes - Cythonized Python version
    """
    # Set upper limit to 1000
    Range  = min(Range, 1000)
    # Initialize working variables
    primes = [0]*Range
    nPrime = 0
    number = 2
    while nPrime < Range:
        # Check if *number* is a prime
        i = 0
        while i < nPrime and number % primes[i] != 0:
            i = i + 1
        # Collect a prime number
        if i == nPrime:
            primes[nPrime] = number
            nPrime = nPrime + 1
        number = number + 1
    return primes

def prime3(int Range):
    """
    Find the first few primes - Cython version
    """
    # Set upper limit to 1000
    Range  = min(Range, 1000)
    # Initialize working variables
    primes = [0]*Range
    cdef int work[1000]
    cdef int nPrime = 0
    cdef int number = 2
    cdef int i
    while nPrime < Range:
        # Check if *number* is a prime
        i = 0
        while i < nPrime and number % work[i] != 0:
            i = i + 1
        # Collect a prime number
        if i == nPrime:
            primes[nPrime] = number
            work[nPrime] = number
            nPrime = nPrime + 1
        number = number + 1
    return primes
