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
    p = 1
    m = 0
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
    q1 = mpmath.exp(a*tau*1j)
    q = q1**24
    s = 1
    qs = mpmath.mpc(1, 0)
    qn = 1
    des = mpmath.mpf(10)**(-floatpre)
    while abs(qs) > des:
        t = -q*qn*qn*qs
        qn = qn*q
        qs = qn*t
        s += t + qs
    return p*q1*s

if __name__=='__main__':
    print(dedekind(13+4j, 8))
