'''
The modified Connacchia-Smith Algorithm. Implemented through Algorithm 2.3.13.
'''
from jacobi import jacobi
from nzmath.arith1 import modsqrt, floorsqrt, issquare
import mpmath

def cornacchia_smith(p, d):
    '''
    modified Cornacchia's Algorithm to solve a^2 + b^2 |D| = 4p for a and b
    Args:
        p:
        d:

    Returns:
        a, b such that a^2 + b^2 |D| = 4p

    '''
    # check input
    if not -4 * p < d < 0:
        raise ValueError(" -4p < D < 0 not true.")
    elif not (d % 4 in {0, 1}):
        raise ValueError(" D = 0, 1 (mod 4) not true.")
        
    # case where p=2
    if p == 2:
        r = sqrt(d + 8)
        if r != -1:
            return r, 1
        else:
            return None
    # test for solvability
    if jacobi(d % p, p) < 1:
        return None
    
    x = modsqrt(d, p)
    if (x % 2) != (d % 2):
        x = p - x
    # euclid chain
    a, b = (2*p, x)
    c = floorsqrt(4 * p)
    while b > c:
        a, b = b, a % b 
    
    t = 4 * p - b*b
    if t % (-d) != 0:
        return None
    if not issquare(t/(-d)):
        return None
    return b, int(mpmath.sqrt(t / -d))
    

def sqrt(x):
    """trial method for calculating square root. """
    i = 1
    while (i * i < x):
        i += 1
    if i * i == x:
        return i
    return -1

'''
def modSqrt(a, p):
    a = a % p
    if p % 8 in {3, 7}:
        x = pow(a, (p+1)/4, p)
        return x
    if p % 8 == 5:
        x = pow(a, (p+3)/8, p)
        c = pow(x, 2, p)
        if c != a:
            x = x * pow(2, (p-1)/4, p) 
            x = x % p
        return x
    if p % 8 == 1:
        d = random.randrange(2, p)
        while (jacobi(d, p) != -1):
            d = random.randrange(2, p)
'''


if __name__ == '__main__':
    print (cornacchia_smith(7, -3))
