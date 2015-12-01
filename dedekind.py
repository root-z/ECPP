import mpmath

def nearest_integer(x):
    """
    Round x to nearest integer.
    """
    return mpmath.floor(x + 0.5)


def dedekind(tau, floatpre):
    """
    Algorithm 22 (Dedekind eta)
    Input : tau in the upper half-plane, k in N
    Output : eta(tau)
    """
    a = 2 * mpmath.pi / mpmath.mpf(24)
    b = mpmath.exp(mpmath.mpc(0, a))
    # b = e^(2pi*i/24)
    p = 1
    m = 0
    # What the hell is happening in this for loop
    '''Removing the integer part of the tau.real '''
    while m <= 0.999:
        n = nearest_integer(tau.real)
        if n != 0:
            tau -= n
            p *= b**n
        m = tau.real*tau.real + tau.imag*tau.imag
        if m <= 0.999:
            ro = mpmath.sqrt(mpmath.power(tau, -1)*1j)
            if ro.real < 0:
                ro = -ro
            p = p*ro
            tau = (-p.real + p.imag*1j) / m
    print('tau=', tau, '\n p =', p)
    q1 = mpmath.exp(a*tau*1j)
    # q1 = e^(2pi*tau*i/24)
    q = q1**24
    # q = e^(2pi*tau*i)
    s = 1
    qs = mpmath.mpc(1, 0)
    qn = 1
    des = mpmath.mpf(10)**(-floatpre)
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

if __name__=='__main__':
    print(dedekind(0.2 + 0.1j, 8))
