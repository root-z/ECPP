'''
author: Z.Lin
'''

import random
from jacobi import jacobi
from cornacchia_smith import cornacchia_smith
from hilbert import hilbert
from nzmath import equation
from nzmath.arith1 import inverse
from nzmath import factor
from nzmath import prime
import mpmath
from nzmath import ecpp

small_primes = factor.mpqs.eratosthenes(10**6)  # for small_primes


def atkin_morain(n):
    d, ms = choose_discriminant(n)
    m_found = False
    while not m_found:
        for m in ms:
            factors = factor_orders(m, n)
            if factors is not None:
                k, q = factors
                m_found = True
                break
        # if no proper m can be found. Go back to choose_discriminant()
        d, ms = choose_discriminant(n, d)
    params = curve_parameters(d, n)

    # Test to see if the order of the curve is really m


def factor_orders(m, n):
    '''
    m = kq
    Args:
        m: for factorizations
        n: the integer for prime proving. Used to compute bound

    Returns:

    '''
    k = 1
    q = m
    bound = (mpmath.power(n, 1/4.0) + 1) ** 2
    for p in small_primes:
        '''
        Check again.
        '''
        # if q becomes too small, fails anyway
        if q <= bound:
            break
        # if p is over the bound, and p is a factor of q. A solution is found.
        if p > bound and q % p == 0:
            k *= q/p
            q = p
            break
        # When q is > bound. Try to factorize through small primes.
        while q % p == 0:
            q = q/p
            k *= p

    if q <= bound or k == 1:
        return None
    if not prime.millerRabin(q):
        return None
    return k, q


def choose_discriminant(n, start=0):
    '''
    First step of Algorithm
    Args:
        n:
        start:

    Returns:
        d: The discriminant
        ms: The possible orders

    '''
    d = gen_discriminant(start)
    uv = cornacchia_smith(n, d)
    while jacobi(d, n) != 1 or uv is None:
        d = gen_discriminant(d)
        uv = cornacchia_smith(n, d)
    u, v = uv

    default = [n+1+u, n+1-u]
    if d == -4:
        ms = default + [n+1+2*v, n+1-2*v]
    elif d == -3:
        ms = default + [n+1+(u+3*v)/2, n+1-(u+3*v)/2, n+1+(u-3*v)/2, n+1-(u-3*v)/2]
    else:
        ms = default

    return d, ms

def curve_parameters(d, p):
    '''
    Modified Algorithm 7.5.9 for the use of ecpp
    Args:
        d: discriminant
        p: number for prime proving

    Returns:
        a list of (a, b) parameters for
    '''
    g = gen_QNR(p)
    u, v = cornacchia_smith(p, d)
    # go without the check for result of cornacchia because it's done by previous methods.
    if jacobi(d, p) != 1:
        raise ValueError('jacobi(d, p) not equal to 1.')

        # check for -d = 3 or 4
    # Choose one possible output
    # Look at param_gen for comparison.
    answer = []

    if d == -3:
        x = -1
        for i in range(0, 6):
            answer.append((0, x))
            x = (x * g) % p
        return answer

    if d == -4:
        x = -1
        for i in range(0, 4):
            answer.append((x, 0))
            x = (x * g) % p
        return answer

    # otherwise compute the hilbert polynomial
    _, t, _ = hilbert(d)
    s = [i % p for i in t]
    j = equation.root_Fp(s, p) # Find a root for s in Fp. Algorithm 2.3.10
    c = j * inverse(j - 1728, p) % p
    r = -3 * c % p
    s = 2 * c % p

    return [(r, s), (r * g * g % p, s * (g**3) % p)]

def generate_curve(p):
    '''
    Essentially Algorithm 7.5.9
    Args:
        p:

    Returns:
        parameters a, b for the curve
    '''
    # calculate quadratic nonresidue
    g = gen_QNR(p)
    # find discriminant
    new_d = gen_discriminant(0)
    uv = cornacchia_smith(p, new_d)
    while jacobi(new_d, p) != 1 or uv is None:
        new_d = gen_discriminant(new_d)
        uv = cornacchia_smith(p, new_d)
    u, v = uv   # storing the result of cornacchia. u^2 + v^2 * |D| = 4*p

    # check for -d = 3 or 4
    # Choose one possible output
    # Look at param_gen for comparison.
    answer = []

    if new_d == -3:
        x = -1
        for i in range(0, 6):
            answer.append((0, x))
            x = (x * g) % p
        return answer

    if new_d == -4:
        x = -1
        for i in range(0, 4):
            answer.append((x, 0))
            x = (x * g) % p
        return answer

    # otherwise compute the hilbert polynomial
    _, t, _ = hilbert(new_d)
    s = [i % p for i in t]
    j = equation.root_Fp(s, p) # Find a root for s in Fp. Algorithm 2.3.10
    c = j * inverse(j - 1728, p) % p
    r = -3 * c % p
    s = 2 * c % p

    return [(r, s), (r * g * g % p, s * (g**3) % p)]


def gen_discriminant(start=0):
    # really should consider generating list of discriminants
    # find a fundamental discriminant. Use start as a starting point.
    d = start
    # compute the odd part
    d_found = False
    while not d_found:
        '''find the next fundamental discriminant'''
        d -= 1
        odd = odd_part(-d)
        # check if the odd part is square free.
        if not square_free_odd(odd):
            continue
        if not ((-d) % 16 in {3, 4, 7, 8, 11, 15}):
            continue
        return d
        # Have the next discriminant. See if it is good


def odd_part(n):
    #compute the odd part of a number
    oddPart = n
    while (oddPart % 2 == 0):
        oddPart = oddPart//2
    return oddPart


def square_free_odd(n):
    #check if the odd part is square free. No efficient algorithm is known.
    x = 1
    x_square = x * x
    while x_square <= n:
        if x_square == n:
            return False
        else:
            x += 2
            x_square = x * x
    return True
    

def gen_QNR(p):
    '''
    generating quardratic residue
    p -- prime
    '''
    #generate random number g
    g = random.randrange(2, p)
    #1. g has to be quadratic nonresidue
    #2. repeat if p = 1 (mod 3) and g^((p-1)/3) = 1 (mod p)
    while(jacobi(g, p) != -1 or (p % 3 == 1 and pow(g, (p-1)/3, p) == 1)):
        g = random.randrange(2, p)
    return g
    
if __name__=='__main__':
    print(gen_QNR(17))

    #primes = factor.mpqs.eratosthenes(10 ** 6)
    #print primes[0:100]
    print factor_orders(49, 1)

