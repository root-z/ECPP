'''
Compute Hilbert Class Polynomials
'''
from mpmath import floor, sqrt, power
from dedekind import dedekind
import mpmath
import nzmath.quad as quad


def solve(d):
    # initialize
    t = [1]
    b = d % 2
    r = floor(sqrt((-d)/3))-1
    h = 0
    red = set()

    h1, reduced_forms = quad.class_group(d)
    print h1
    harmonic_sum = sum(1/mpmath.mpf(form[0]) for form in reduced_forms)
    floatpre = floor(mpmath.pi*mpmath.sqrt(-d)*harmonic_sum / mpmath.log(10)+0.5) + 10
    mpmath.mp.dps = floatpre



    # outer loop
    while b <= r:
        m = (b*b - d) / 4
        m_sqrt = int(floor(sqrt(m)))
        for a in range(1, m_sqrt+1):
            if m % a != 0:
                continue
            c = m/a
            if b > a:
                continue
            # optional polynomial setup
            tau = (-b + 1j * sqrt(-d)) / (2*a)
            f = power(dedekind(2 * tau, floatpre) / dedekind(tau, floatpre), 24)
            j = power((256 * f + 1), 3) / f

            if b==a or c==a or b==0:
                # T = T * (X-j)
                t = polynomial_mul(t, [-j, 1])
                h += 1
                red.add((a, b, c))
            else:
                poly = [j.real * j.real + j.imag * j.imag, -2 * j.real, 1]
                t = polynomial_mul(t, poly)
                h += 2
                red.add((a, b, c))
                red.add((a, -b, c))
        b += 1

    return h, [floor(mpmath.re(p) + 0.5) for p in t], red


def delta(q):
    return q


def polynomial_mul(p1, p2):
    '''
    Used to Compute T = T * (X-j)
    '''
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
    print solve(-4)