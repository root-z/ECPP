"""
All methods defined here return one of a factor of given integer.
When 1 is returned, the method has failed to factor,
but 1 is a factor anyway.

'verbose' boolean flag can be specified for verbose reports.
"""

import logging
import nzmath.arith1 as arith1
import nzmath.bigrandom as bigrandom
import nzmath.gcd as gcd
import nzmath.prime as prime

_log = logging.getLogger('nzmath.factor.find')

# Pollard's rho method
def rhomethod(n, **options):
    """
    Find a non-trivial factor of n using Pollard's rho algorithm.
    The implementation refers the explanation in C.Pomerance's book.
    """
    # verbosity
    verbose = options.get('verbose', False)
    if not verbose:
        _silence()

    if n <= 3:
        return 1

    g = n
    while g == n:
        # x^2 + a is iterated. Starting value x = u.
        a = bigrandom.randrange(1, n-2)
        u = v = bigrandom.randrange(0, n-1)
        _log.info("%d %d" % (a, u))
        g = gcd.gcd((v**2 + v + a) % n - u, n)
        while g == 1:
            u = (u**2 + a) % n
            v = ((pow(v, 2, n) + a)**2 + a) % n
            g = gcd.gcd(v - u, n)
    if not verbose:
        _verbose()
    return g

# p-1 method
def pmom(n, **options):
    """
    This function tries to find a non-trivial factor of n using
    Algorithm 8.8.2 (p-1 first stage) of Cohen's book.
    In case of N = pow(2,i), this program will not terminate.
    """
    # verbosity
    verbose = options.get('verbose', False)
    if not verbose:
        _silence()

    # initialize
    x = y = 2
    primes = []
    if 'B' in options:
        B = options['B']
    else:
        B = 10000

    for q in prime.generator():
        primes.append(q)
        if q > B:
            if gcd.gcd(x-1, n) == 1:
                if not verbose:
                    _verbose()
                return 1
            x = y
            break
        q1 = q
        l = B // q
        while q1 <= l:
            q1 *= q
        x = pow(x, q1, n)
        if len(primes) >= 20:
            if gcd.gcd(x-1, n) == 1:
                primes, y = [], x
            else:
                x = y
                break

    for q in primes:
        q1 = q
        while q1 <= B:
            x = pow(x, q, n)
            g = gcd.gcd(x-1, n)
            if g != 1:
                if not verbose:
                    _verbose()
                if g == n:
                    return 1
                return g
            q1 *= q

def trialDivision(n, **options):
    """
    Return a factor of given integer by trial division.

    options can be either:
    1) 'start' and 'stop' as range parameters.
    2) 'iterator' as an iterator of primes.
    If both options are not given, prime factor is searched from 2
    to the square root of the given integer.
    """
    # verbosity
    verbose = options.get('verbose', False)
    if not verbose:
        _silence()

    if 'start' in options and 'stop' in options:
        if 'step' in options:
            trials = range(options['start'], options['stop'], options['step'])
        else:
            trials = range(options['start'], options['stop'])
    elif 'iterator' in options:
        trials = options['iterator']
    elif n < 1000000:
        trials = prime.generator_eratosthenes(arith1.floorsqrt(n))
    else:
        trials = prime.generator()

    limit = arith1.floorsqrt(n)
    for p in trials:
        if limit < p:
            break
        if 0 == n % p:
            if not verbose:
                _verbose()
            return p
    if not verbose:
        _verbose()
    return 1

def _silence():
    """
    Stop verbose outputs.
    """
    _log.setLevel(logging.NOTSET)

def _verbose():
    """
    Stop silencing.
    """
    _log.setLevel(logging.DEBUG)
