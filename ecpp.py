'''
author: Z.Lin
'''

import random
from jacobi import jacobi
from cornacchia_smith import cornacchia_smith
from hilbert import hilbert
from nzmath import equation
from nzmath.arith1 import inverse

def atkin_morain():
    pass


def generate_curve(p):
    '''

    Args:
        p:

    Returns:
        parameters a, b for the curve
    '''
    # calculate quadratic nonresidue
    g = gen_QNR(p)
    # find discriminant
    new_d = gen_discriminant(p)
    uv = cornacchia_smith(p, new_d)
    while jacobi(new_d, p) != 1 or uv is None:
        new_d = gen_discriminant(p, new_d)
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
    _, t, _ = hilbert(d)
    s = [i % p for i in t]
    j = equation.root_Fp(s, p) # Find a root for s in Fp. Algorithm 2.3.10
    c = j * inverse(j - 1728, p) % p
    r = -3 * c % p
    s = 2 * c % p

    return [(r, s), (r * g * g % p, s * (g**3) % p)]

    
    
def gen_discriminant(p, start = 0):
    # really should consider generating list of discriminants
    # find a fundamental discriminant. Use start as a starting point.
    d = start
    # compute the odd part
    dfound = False
    while not dfound:
        '''find the next fundamental discriminant'''
        d -= 1
        odd = odd_part(-d)
        #check if the odd part is square free.
        if not square_free_odd(odd):
            continue
        if not ((-d) % 16 in {3, 4, 7, 8, 11, 15}):
            continue
        return d
        #Have the next discriminant. See if it is good
        '''
        if (jacobi(d % p , p) != 1):
            continue
        '''
        #connarchia

def odd_part(n):
    #compute the odd part of a number
    oddPart = n
    while (oddPart % 2 == 0):
        oddPart = oddPart//2
    return oddPart
 
def square_free_odd(n):
    #check if the odd part is square free. No efficient algorithm is known.
    x = 1
    xSquare = x * x
    while (xSquare <= n):
        if xSquare == n:
            return False
        else:
            x += 2
            xSquare = x * x
    return True
    
 
'''
generating quardratic residue
p -- prime
'''
def gen_QNR(p):
    #generate random number g
    g = random.randrange(2, p)
    #1. g has to be quadratic nonresidue
    #2. repeat if p = 1 (mod 3) and g^((p-1)/3) = 1 (mod p)
    while(jacobi(g, p) != -1 or (p % 3 == 1 and pow(g, (p-1)/3, p) == 1)):
        g = random.randrange(2, p)
    return g
    
if __name__=='__main__':
    print(gen_QNR(17))
