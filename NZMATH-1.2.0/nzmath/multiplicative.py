"""
Multiplicative number theoretic functions.
"""

import nzmath.factor.misc as factor_misc
import nzmath.prime as prime

def euler(n):
    """
    Euler totient function.
    It returns the number of relatively prime numbers to n smaller than n.
    """
    if n == 1:
        return 1
    if prime.primeq(n):
        return n - 1
    t = 1
    for p, e in factor_misc.FactoredInteger(n):
        if e > 1:
            t *= pow(p, e - 1) * (p - 1)
        else:
            t *= p - 1
    return t

def moebius(n):
    """
    Moebius function.
    It returns:
      -1  if n has odd distinct prime factors,
       1  if n has even distinct prime factors, or
       0  if n has a squared prime factor.
    """
    if n == 1:
        return 1
    if prime.primeq(n):
        return -1
    m = 1
    for p, e in factor_misc.FactoredInteger(n):
        if e > 1:
            return 0
        m = -m
    return m

def sigma(m, n):
    """
    Return the sum of m-th powers of the factors of n.
    """
    if n == 1:
        return 1
    if prime.primeq(n):
        return 1 + n**m
    s = 1
    for p, e in factor_misc.FactoredInteger(n):
        t = 1
        for i in range(1, e + 1):
            t += (p**i)**m
        s *= t
    return s

