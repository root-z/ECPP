'''
Compute Hilbert Class Polynomials
'''
from mpmath import floor, sqrt, power
from dedekind import dedekind
import mpmath
import nzmath.ecpp


def solve(d):
    # initialize
    t = [1]
    b = d % 2
    r = floor(sqrt((-d)/3))
    h = 0
    red = {}

    # outer loop
    while (b <= r) :
        m = (b*b - d) / 4
        a = 0
        while a * a <= m:
            if m % a != 0:
                continue
            c = m/a
            if b > a:
                continue
            # optional polynomial setup
            tau = (-b + 1j * sqrt(-d)) / (2*a)
            f = power(dedekind(2 * tau) / dedekind(tau), 24)
            j = power((256 * f + 1), 3) / f

            if b==a or c==a or b==0:
                t = polynomialMul(t, [-1j, 1])

def delta(q):
    return q

def polynomialMul(p1, p2):
    if len(p1)==0 or len(p2)==0:
        raise ValueError('Polynomial Array empty.')
    m = [0]* (len(p1) + len(p2) -1)
    for i in range(0, len(p1)):
        for j in range(0, len(p2)):
            m[i+j] += p1[i] * p2[j]
    return m

if __name__ == '__main__':
    '''
    print(mpmath.sqrt(1 + 1j))
    print(mpmath.sqrt(211))
    print power(1j, 2)
    print nzmath.ecpp.hilbert(-15)
    '''
    print polynomialMul([0,1], [-1j, 1])
