'''
Compute Hilbert Class Polynomials
'''
from mpmath import floor, sqrt, power
import mpmath
import nzmath.quad as quad


round = lambda x : mpmath.floor(x + 0.5)

def hilbert(d):
    '''

    Args:
        d:

    Returns:

    '''
    # initialize
    t = [1]
    b = d % 2
    r = floor(sqrt((-d)/3))
    h = 0
    red = set()

    h1, reduced_forms = quad.class_group(d)
    # print h1
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
        b += 2

    return h, [int(floor(mpmath.re(p) + 0.5)) for p in t], red


def delta(q):
    return q


def dedekind(tau, floatpre):
    """
    Algorithm 22 (Dedekind eta)
    Input : tau in the upper half-plane, k in N
    Output : eta(tau)
    """
    a = 2 * mpmath.pi / mpmath.mpf(24)
    b = mpmath.exp(mpmath.mpc(0, a))
    #print 'b=', b
    # b = e^(2pi*i/24)
    p = 1   # What is this for?
    m = 0

    # Now I know which two functions are used in the for loops
    #But why keep the absolute value greater than 1?
    while m <= 0.999:
        n = round(tau.real)
        if n != 0:
            tau -= n
            # Removing the integer part of the tau.real
            p *= b**n
        m = tau.real*tau.real + tau.imag*tau.imag
        # Keep m >= 1? Why?
        if m <= 0.999:
            ro = mpmath.sqrt(mpmath.power(tau, -1)*1j)
            # ro = sqrt((tau^-1)*i)
            if ro.real < 0:
                ro = -ro
            p = p*ro
            tau = (-p.real + p.imag*1j) / m
    #print 'tau=', tau, '\n p =', p

    q1 = mpmath.exp(a*tau*1j)
    #print 'q1=', q1
    # q1 = e^(2pi*tau*i/24)
    q = q1**24
    # q = e^(2pi*tau*i)
    s = 1
    qs = mpmath.mpc(1, 0)
    qn = 1
    des = mpmath.mpf(10)**(-floatpre)

    '''I believe this is the sum function listed in the book and website'''
    # http://mathworld.wolfram.com/DedekindEtaFunction.html
    while abs(qs) > des:
        t = -q*qn*qn*qs
        # t = -q, q^5
        qn = qn*q
        # qn = q, q^2
        qs = qn*t
        # qs = -q^2, q^7
        s += t + qs
        # s = 1 - q - q^2, 1-q-q^2+q^5+q^7

    return p*q1*s
    # Compare to wolfram alpha the result is correct.


def polynomial_mul(p1, p2):
    '''
    Used to Compute T = T * (X-j)
    '''
    if len(p1) == 0 or len(p2) == 0:
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
    print hilbert(-23)