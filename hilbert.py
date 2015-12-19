'''
Compute Hilbert Class Polynomials
'''
from mpmath import *
import mpmath
import ecpp

round = lambda x: mpmath.floor(x + 0.5)


def reduced_form(d):
    '''

    Args:
        d:

    Returns:

    '''
    # initialize
    b = d % 2
    r = floor(sqrt((-d)/3))
    h = 0
    red = set()

    # outer loop
    while b <= r:
        m = (b*b - d) / 4
        m_sqrt = int(floor(sqrt(m)))
        for a in range(1, m_sqrt+1):
            if m % a != 0:
                continue
            c = m / a
            if b > a:
                continue
            # optional polynomial setup
            if b==a or c==a or b==0:
                # T = T * (X-j)
                h += 1
                red.add((a, b, c))
            else:
                h += 2
                red.add((a, b, c))
                red.add((a, -b, c))
        b += 2

    return red


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

    reduced_forms = reduced_form(d)    # print h1
    a_inverse_sum = sum(1/mpf(form[0]) for form in reduced_forms)
    precision = round(pi*sqrt(-d)*a_inverse_sum / log(10)) + 10
    mpmath.mp.dps = precision

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
            f = power(dedekind(2 * tau, precision) / dedekind(tau, precision), 24)
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

    if red != reduced_forms:
        raise ValueError('Reduced form inconsistent.')

    return h, [int(floor(mpmath.re(p) + 0.5)) for p in t], red


def delta(q):
    return q


def dedekind(tau, precision):
    """

    Args:
        tau:
        precision:

    Returns:

    """
    # a = 2 * mpmath.pi / mpmath.mpf(24)
    # b = mpmath.exp(mpmath.mpc(0, a))

    x = exp(mpc(0, 2 * pi / mpf(24)))
    # b = e^(2pi*i/24)
    outer = 1   # What is this for?
    absolute = 0

    # functional equations
    while absolute <= 1 - 0.1**5:
        real_tau = round(tau.real)
        if real_tau != 0:
            tau -= real_tau
            outer *= x ** real_tau
        absolute = fabs(tau)
        if absolute > 1 - 0.1**5:
            break
        ro = mpmath.sqrt(mpmath.power(tau, -1)*1j)
        # ro = sqrt((tau^-1)*i)
        if ro.real < 0:
            ro = -ro
        outer = outer*ro
        tau = (-outer.real + outer.imag*1j) / absolute
    #print 'tau=', tau, '\n p =', p

    q1 = mpmath.exp((pi/12) * tau * 1j)
    q = q1**24
    # q = e^(2pi*tau*i)
    sum = 1
    qs = mpmath.mpc(1, 0)
    qn = 1
    bound = mpmath.mpf(10)**(-precision-2)

    while fabs(qs) > bound:
        t = -q*qn*qn*qs
        qn *= q
        qs = qn*t
        sum += t + qs

    return outer*q1*sum
    # Compare to wolfram alpha the result is correct.


def polynomial_mul(p1, p2):
    '''
    Used to Compute T = T * (X-j)
    '''
    if len(p1) == 0 or len(p2) == 0:
        raise ValueError('Polynomial Array empty.')
    m = [0] * (len(p1) + len(p2) - 1)
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
    #print hilbert(-23)
    a, b = reduced_form(ecpp.gen_discriminant(-10000))
    print len(a), len(b)