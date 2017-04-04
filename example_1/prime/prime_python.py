def prime1(Range):
    """
    Find the first few primes - Pure Python version
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
