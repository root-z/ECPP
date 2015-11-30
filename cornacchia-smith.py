'''
The modified Connacchia-Smith Algorithm. Implemented through Algorithm 2.3.13.
'''
from jacobi import jacobi
import math
from nzmath.arith1 import modsqrt, floorsqrt, issquare

def solve(p, D):
    #check input
    if not -4 * p < D < 0:
        raise ValueError(" -4p < D < 0 not true.")
    elif not (D % 4 in {0,1}):
        raise ValueError(" D = 0, 1 (mod 4) not true.")
        
    #case where p=2
    if p == 2:
        r = sqrt(D+8)
        if r != -1:
            return (r, 1)
        else:
            return None
    #test for solvability
    if jacobi(D % p, p) < 1:
        return None
    
    x = modsqrt(D, p)
    if ((x % 2) != (D % 2)):
        x = p - x
    #euclid chain
    a, b = (2*p, x)
    c = floorsqrt(4 * p)
    while b > c:
        a, b = b, a % b 
    
    t = 4 * p - D
    if t % (-D) ==  0:
        return None
    if not issquare(t/(-D)):
        return None
    return (b, int(math.sqrt(t/(-D))))
    

def sqrt(x):
    '''trial method for calculating square root. Good enough for purpose'''
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
unfinished
'''


if __name__=='__main__':
    print(solve(103 , -15))
