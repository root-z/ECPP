from __future__ import division
import os
import csv
import logging
import nzmath.arith1 as arith1
import nzmath.bigrandom as bigrandom
import nzmath.bigrange as bigrange
import nzmath.equation as equation
import nzmath.factor as factor
import nzmath.intresidue as intresidue
import nzmath.prime as prime
import nzmath.quad as quad
import nzmath.squarefree as squarefree
import nzmath.compatibility
from nzmath.config import DATADIR, HAVE_MPMATH, HAVE_NET

if HAVE_MPMATH:
    import mpmath

    log = mpmath.log

    def nearest_integer(x):
        """
        Round x to nearest integer.
        """
        return mpmath.floor(x + 0.5)

else:
    import math

    log = math.log

    def nearest_integer(x):
        """
        Round x to nearest integer.
        """
        return math.floor(x + 0.5)


_log = logging.getLogger('nzmath.ecpp')

# mapping from disc to list of order functions.
default_orders = [lambda n, u, v: n + 1 + u,
                  lambda n, u, v: n + 1 - u]
orders = {-3: default_orders + [lambda n, u, v: n + 1 - ((u + 3*v) >> 1),
                                lambda n, u, v: n + 1 + ((u + 3*v) >> 1),
                                lambda n, u, v: n + 1 - ((u - 3*v) >> 1),
                                lambda n, u, v: n + 1 + ((u - 3*v) >> 1)],
          -4: default_orders + [lambda n, u, v: n + 1 + (v << 1),
                                lambda n, u, v: n + 1 - (v << 1)]}
# default absolute value bound of discriminant
DEFAULT_ABS_BOUND = 200000

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

def hilbert_mpmath(D):
    """
    Algorithm 23&24 (Floating-point precision & Hilbert class polynomial)
    Input : a negative fundamental discriminant D
    Output : the class number h(D) and the Hilbert class polynomial T in Z[X]
    """
    if D >= 0:
        raise ValueError("Input number must be a negative integer.")
    T = [complex(1), 0]
    mpmath.mp.default()
    h, reduced_forms = quad.class_group(D)
    harmonic_sum = sum(1/mpmath.mpf(form[0]) for form in reduced_forms)
    floatpre = nearest_integer(mpmath.pi*mpmath.sqrt(-D)*harmonic_sum / mpmath.log(10)) + 10
    mpmath.mp.dps = floatpre
    sqrtD = mpmath.sqrt(D)
    for a, b, c in reduced_forms:
        if b < 0:continue
        tau = (-b + sqrtD)/(2*a)
        f = mpmath.power(dedekind(2*tau, floatpre) / dedekind(tau, floatpre), 24)
        j = ((256*f + 1)**3) / f
        TC = tuple(T)
        if b == a or c == a or b == 0:
            T[0] = TC[0]*(-j)
            if len(T) > 2:
                for i in range(1, len(T)):
                    T[i] = TC[i - 1] + (TC[i]*(-j))
                T.append(1)
            elif T[1] != 0:
                T[1] = TC[0] + (TC[1]*(-j))
                T.append(1)
            else:
                T[1] = 1
        else:
            rej = -2*j.real
            abj = j.real**2 + j.imag**2
            T[0] = TC[0]*abj
            T[1] = TC[0]*rej + TC[1]*abj
            if len(T) > 2:
                for i in range(2, len(T)):
                    T[i] = TC[i-2] + TC[i-1]*rej + TC[i]*abj
            T.extend([TC[-2] + TC[-1]*rej, 1])
    T = [(int(nearest_integer(c.real)) if hasattr(c, 'real') else c) for c in T]
    return h, T

def hilbert(D):
    """
    wrapper to split mpmath part and others.
    """
    if HAVE_NET:
        try:
            import urllib2
            url = 'http://hilbert-class-polynomial.appspot.com/%d/' % D
            result = urllib2.urlopen(url).read()
            if result == 'null':
                pass
            else:
                hcp = eval(result)
                return len(hcp) - 1, hcp
        except Exception:
            pass

    if HAVE_MPMATH:
        return hilbert_mpmath(D)

    raise prime.ImplementLimit("impliment limit")


def outputHP(absbound=DEFAULT_ABS_BOUND):
    """
    Output Hilbert class polynomials to file 'Hilbert_Poly.txt'.
    """
    f = open('Hilbert_Poly.txt', 'a')
    for disc in disc_gen(absbound):
        h, T = hilbert(disc)
        f.write("%s %s %s\n" % (str(disc), str(h), str(T)))
    f.close()


class Elliptic(object):
    """
    The class is for operations of elliptic curve points.
    """
    def __init__(self, coefficient, modulus):
        self.modulus = modulus
        self.a, self.b = coefficient
        self.f = lambda x: ((x*x + self.a)*x + self.b) % self.modulus

    def add(self, P, Q):
        """
        Return P + Q
        """
        infinity = [0]
        if P == infinity:
            return Q
        elif Q == infinity:
            return P
        else:
            if P[0] == Q[0]:
                if not (P[1] + Q[1]):
                    return infinity
                s = 3*P[0]*P[0] + self.a
                t = 2*P[1]
            else:
                s = Q[1] - P[1]
                t = Q[0] - P[0]
            m = s/t
            x_3 = m*m - P[0] - Q[0]
            return [x_3, m*(P[0] - x_3) - P[1]]

    def sub(self, P, Q):
        """
        return P - Q
        """
        infinity = [0]
        if Q == infinity:
            return P
        R = [Q[0], -Q[1]]
        if P == infinity:
            return R
        else:
            return self.add(P, R)

    def mul(self, k, P):
        """
        this returns [k]*P
        """
        if k >= 0:
            l = arith1.expand(k, 2)
            Q = [0]
            for j in range(len(l) - 1, -1, -1):
                Q = self.add(Q, Q)
                if l[j] == 1:
                    Q = self.add(Q, P)
            return Q
        else:
            return self.sub([0], self.mul(-k, P))

    def choose_point(self):
        """
        Choose point on E_{a,b}(Z_n)
        Algorithm 27 (Atkin-morain ECPP) Step5
        """
        n, f = self.modulus, self.f
        x = bigrandom.randrange(n)
        Q = f(x)
        while arith1.legendre(Q, n) == -1:
            x = bigrandom.randrange(n)
            Q = f(x)
        y = arith1.modsqrt(Q, n)
        return [intresidue.IntegerResidueClass(t, n) for t in (x, y)]


def cmm(p):
    """
    Algorithm 25 (CM methods)
    Input : a prime number p (>=3)
    Output : curve parameters for CM curves
    """
    disc, g = _cmm_find_disc(p)[:2]
    return [param for param in param_gen(disc, p, g)]

def cmm_order(p):
    """
    Algorithm 25 (CM methods)
    Input : a prime number p (>=3)
    Output : curve parameters for CM curves and orders
    """
    disc, g, u, v = _cmm_find_disc(p)
    m = [f(p, u, v) for f in orders.get(disc, default_orders)]
    params = []
    for param in param_gen(disc, p, g):
        E = Elliptic(param, p)
        base_point = E.choose_point()
        for mi in m:
            Q = E.mul(mi, base_point)
            if Q == [0]:
                params.append(param + (mi,))
                m.remove(mi)
                break
    return params

def _cmm_find_disc(p, absbound=DEFAULT_ABS_BOUND):
    """
    Input : a prime number p (>=3)
    Output : (disc, g, u, v) where:
            - disc: discriminant appropriate for CM method
            - g: quasi primitive element of p
            - u, v: a solution of u**2 + disc * v**2 = 4*p
    """
    for disc in disc_gen(absbound):
        if arith1.legendre(disc, p) == 1:
            try:
                u, v = cornacchiamodify(disc, p)
                g = quasi_primitive(p, disc == -3)
                return disc, g, u, v
            except ValueError:
                continue
    raise ValueError("no discriminant found")

def cornacchiamodify(d, p):
    """
    Algorithm 26 (Modified cornacchia)
    Input : p be a prime and d be an integer such that d < 0 and d > -4p with
            d = 0, 1 (mod 4)
    Output : the solution of u^2 -d * v^2 = 4p.
    """
    q = 4 * p
    if (d >= 0) or (d <= -q):
        raise ValueError("invalid input")
    if p == 2:
        b = arith1.issquare(d + 8)
        if b:
            return (b, 1)
        else:
            raise ValueError("no solution")
    if arith1.legendre(d, p) == -1:
        raise ValueError("no solution")
    x0 = arith1.modsqrt(d, p)
    if (x0 - d) & 1:
        x0 = p - x0
    a = 2 * p
    b = x0
    l = arith1.floorsqrt(q)
    while b > l:
        a, b = b, a % b
    c, r = divmod(q - b * b, -d)
    if r:
        raise ValueError("no solution")
    t = arith1.issquare(c)
    if t:
        return (b, t)
    else:
        raise ValueError("no solution")

def quasi_primitive(p, disc_is_minus3):
    """
    Return g s.t.
    0) 1 < g < p
    1) quadratic nonresidue modulo p
    and
    2) cubic nonresidue modulo p if p == 1 mod 3
    with p >= 3, a prime.

    For p, not a prime, condition 1 means the Jacobi symbol
    (g/p) != -1, and condition 2 means g doesn't have order
    dividing (p - 1)//3.

    if disc_is_minus3 is True, cubic nonresidue check becomes
    stricter so that g**((p-1)//3) must match with one of the
    primitive third roots of unity.
    """
    third_order = (p - 1) // 3
    for g in bigrange.range(2, p):
        # if g is not quadratic nonresidue, then skipped.
        legendre_jacobi = arith1.legendre(g, p)
        if legendre_jacobi != -1:
            continue
        # if p != 1 mod 3 or g is cubic nonresidue, it's what we need.
        if p % 3 != 1:
            return g
        cubic_symbol = pow(g, third_order, p)
        if cubic_symbol != 1:
            # g is cubic nonresidue.
            if disc_is_minus3 and pow(cubic_symbol, 3, p) != 1:
                # stricter check
                continue
            return g
    else:
        raise ValueError("p is not prime.")

def param_gen(disc, n, g):
    """
    Return generator to generate curve params.
    """
    if disc == -3:
        return _generate_params_for_minus3(n, g)
    elif disc == -4:
        return _generate_params_for_minus4(n, g)
    else:
        return _generate_params_for_general_disc(disc, n, g)

def _generate_params_for_minus3(n, g):
    """
    Generate curve params for discriminant -3.
    """
    yield (0, n - 1)
    gk = -1
    for i in range(5):
        gk = (gk*g)%n
        yield (0, gk)

def _generate_params_for_minus4(n, g):
    """
    Generate curve params for discriminant -4.
    """
    yield (n - 1, 0)
    gk = -1
    for i in range(3):
        gk = (gk*g)%n
        yield (gk, 0)

def _generate_params_for_general_disc(disc, n, g):
    """
    Generate curve params for discriminant other than -3 or -4.
    """
    jr = equation.root_Fp([c % n for c in hilbert(disc)[-1]], n)
    c = (jr * arith1.inverse(jr - 1728, n)) % n
    r = (-3 * c) % n
    s = (2 * c) % n
    yield (r, s)
    g2 = (g * g) % n
    yield ((r * g2) % n, (s * g2 * g) % n)


def ecpp(n, era=None):
    """
    Algorithm 27 (Atkin-Morain ECPP)
    Input : a natural number n
    Output : If n is prime, retrun True. Else, return False
    """
    _log.info('n = %d' % n)
    disc_and_order_info = _ecpp_find_disc(n, era)
    if not disc_and_order_info:
        return False
    # disc, m = k*q, era
    disc, m, k, q, era = disc_and_order_info
    # non(quadratic|cubic) residue
    g = quasi_primitive(n, disc == -3)

    try:
        # find param corresponding to order m
        for param in param_gen(disc, n, g):
            E = Elliptic(param, n)
            P = E.choose_point()
            if E.mul(m, P) == [0]:
                break
        # find a point [k]P is not the unit.
        U = E.mul(k, P)
        while U == [0]:
            P = E.choose_point()
            U = E.mul(k, P)
        # compute [q]([k]P)
        V = E.mul(q, U)
    except ZeroDivisionError:
        return False
    if V == [0]:
        if q > 1000000000:
            return ecpp(q, era)
        elif prime.trialDivision(q):
            return True
        else:
            return False
    else:
        return False

def _ecpp_find_disc(n, era=None):
    """
    Return (disc, m, k, q, era) where:
    - disc: appropriate discriminant for ecpp
    - m: one of possible orders for disc
    - k: smooth part of m
    - q: rough part of m which is greater than (n^(1/4) + 1)^2
    - era: list of small primes
    """
    up = (arith1.floorpowerroot(n, 4) + 1)**2
    for disc in disc_gen(n):
        legendre = arith1.legendre(disc, n)
        if legendre == 0:
            return False
        elif legendre == -1:
            continue
        # legendre == 1:
        try:
            u, v = cornacchiamodify(disc, n)
        except ValueError:
            continue
        for m in (f(n, u, v) for f in orders.get(disc, default_orders)):
            k, q, era = _factor_order(m, up, era)
            if k:
                return disc, m, k, q, era
    else:
        # no usable discriminat found
        return False

def _factor_order(m, u, era=None):
    """
    Return triple (k, q, primes) if m is factored as m = kq,
    where k > 1 and q is a probable prime > u.
    u is expected to be (n^(1/4)+1)^2 for Atkin-Morain ECPP.
    If m is not successfully factored into desired form,
    return (False, False, primes).
    The third component of both cases is a list of primes
    to do trial division.

    Algorithm 27 (Atkin-morain ECPP) Step3
    """
    if era is None:
        v = min(500000, int(log(m)))
        era = factor.mpqs.eratosthenes(v)
    k = 1
    q = m
    for p in era:
        if p > u:
            break
        e, q = arith1.vp(q, p)
        k *= p**e
    if q < u or k == 1:
        return False, False, era
    if not prime.millerRabin(q):
        return False, False, era
    return k, q, era

def next_disc(d, absbound):
    """
    Return fundamental discriminant D, such that D < d.
    If there is no D with |d| < |D| < absbound, return False
    """
    # -disc % 16
    negdisc_mod16 = (3, 4, 7, 8, 11, 15)
    for negdisc in bigrange.range(-d + 1, absbound):
        if negdisc & 15 not in negdisc_mod16:
            continue
        if negdisc & 1 and not squarefree.trial_division(negdisc):
            continue
        if arith1.issquare(negdisc):
            continue
        return -negdisc
    return False

def disc_gen(absbound=DEFAULT_ABS_BOUND):
    """
    Generate negative discriminants.

    The order of discriminants depends on how to produce them.
    By default, discriminants are in ascending order of class
    numbers while reading from precomputed data file, then
    in descending order of discriminant itself.

    If discriminant reach the absbound, then StopIteration will be
    raised.  The default value of absbound is DEFAULT_ABS_BOUND.
    """
    csvfile = open(os.path.join(DATADIR, "discriminant.csv"))
    for disc_str in csv.reader(csvfile):
        disc = int(disc_str[0])
        if -disc >= absbound:
            raise StopIteration("absbound reached")
        yield disc
    disc = next_disc(disc, absbound)
    while disc:
        yield disc
        disc = next_disc(disc, absbound)
    raise StopIteration("absbound reached")
