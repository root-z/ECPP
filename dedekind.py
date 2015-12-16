import mpmath
from mpmath import power
from hilbert import dedekind as dd
from mpmath import *

def nearest_integer(x):
    """
    Round x to nearest integer.
    """
    return mpmath.floor(x + 0.5)

def eta(tau, precision):
    """
    Algorithm 22 (Dedekind eta)
    Input : tau in the upper half-plane, k in N
    Output : eta(tau)
    """
    a = 2 * mpmath.pi
    b = mpmath.exp(mpmath.mpc(0, a))
    #print 'b=', b
    # b = e^(2pi*i/24)
    p = 1   # What is this for?
    m = 0

    # Now I know which two functions are used in the for loops
    #But why keep the absolute value greater than 1?
    while m <= 0.999:
        n = nearest_integer(tau.real)
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
        n = nearest_integer(tau.real)
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

if __name__=='__main__':
    floatpre = 8
    tau = 5 + 2j
    print dedekind(tau, floatpre)
    print dd(tau, floatpre)
    # print eta(1 + 2j, 8)
    '''
    f = power(dedekind(2 * tau, floatpre) / dedekind(tau, floatpre), 24)
    f2 = power(dedekind(2 * tau, floatpre), 24) / power(dedekind(tau, floatpre), 24)

    #print f
    #print f2
    '''
    q1 = mpmath.exp((pi/12) * tau * 1j)
    q12 = mpmath.exp((2*pi/24) * tau * 1j)
    print q1, q12
    #print 'q1=', q1
    # q1 = e^(2pi*tau*i/24)
    q = q1**24
    q2 = q12 ** 24
    print q, q2
