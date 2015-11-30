"""
equation -- methods to solve algebraic equations.

If you would like to solve an equation in the complex field and need
more precision, then you can import and use plugins.SETPRECISION.

e.g:
from nzmath.plugins import SETPRECISION
SETPRECISION(200)

Then, the following computations will be done in such precision, if
plugins.PRECISION_CHANGEABLE is true.
"""

from __future__ import division
import logging

import nzmath.arith1 as arith1
import nzmath.bigrandom as bigrandom
import nzmath.bigrange as bigrange
import nzmath.finitefield as finitefield
import nzmath.imaginary as imaginary
import nzmath.poly.uniutil as uniutil

from nzmath.plugins import MATHMODULE as math, CMATHMODULE as cmath, \
     FLOATTYPE as Float, COMPLEXTYPE as Complex

_log = logging.getLogger('nzmath.equation')
_log.setLevel(logging.DEBUG)

# x is (list,tuple)
# t is variable
def e1(x):
    """
    0 = x[0] + x[1]*t
    """
    if x[1] == 0:
        raise ZeroDivisionError("No Solution")
    else:
        return -x[0]/x[1]

def e1_ZnZ(x, n):
    """
    Return the solution of x[0] + x[1]*t = 0 (mod n).
    x[0], x[1] and n must be positive integers.
    """
    try:
        return (-x[0] * arith1.inverse(x[1], n)) % n
    except ZeroDivisionError:
        raise ValueError("No Solution")

def e2(x):
    """
    0 = x[0] + x[1]*t + x[2]*t**2
    """
    c, b, a = x
    d = b**2 - 4*a*c
    if d >= 0:
        sqrtd = math.sqrt(d)
    else:
        sqrtd = cmath.sqrt(d)
    return ((-b + sqrtd)/(2*a), (-b - sqrtd)/(2*a))

def e2_Fp(x, p):
    """
    p is prime
    f = x[0] + x[1]*t + x[2]*t**2
    """
    c, b, a = [_x % p for _x in x]
    if a == 0:
        return [e1_ZnZ([c, b], p)]
    if p == 2:
        solutions = []
        if x[0] & 1 == 0:
            solutions.append(0)
        if (x[0] + x[1] + x[2]) & 1 == 0:
            solutions.append(1)
        if len(solutions) == 1:
            return solutions * 2
        return solutions
    d = b**2 - 4*a*c
    if arith1.legendre(d, p) == -1:
        return []
    sqrtd = arith1.modsqrt(d, p)
    a = arith1.inverse(2*a, p)
    return [((-b+sqrtd)*a)%p, ((-b-sqrtd)*a)%p]

def e3(x):
    """
    0 = x[0] + x[1]*t + x[2]*t**2 + x[3]*t**3
    """
    x3 = Float(x[3])
    a = x[2] / x3
    b = x[1] / x3
    c = x[0] / x3
    p = b - (a**2)/3
    q = 2*(a**3)/27 - a*b/3 + c
    w = (-1 + cmath.sqrt(-3)) / 2
    W = (1, w, w.conjugate())
    k = -q/2 + cmath.sqrt((q**2)/4 + (p**3)/27)
    l = -q/2 - cmath.sqrt((q**2)/4 + (p**3)/27)
    m = k ** (1/3)
    n = l ** (1/3)

    # choose n*W[i] by which m*n*W[i] be the nearest to -p/3
    checker = [abs(3*m*n*z + p) for z in W]
    n = n * W[checker.index(min(checker))]

    sol = []
    for i in range(3):
        sol.append(W[i]*m + W[-i]*n - a/3)
    return sol

def e3_Fp(x, p):
    """
    p is prime
    0 = x[0] + x[1]*t + x[2]*t**2 + x[3]*t**3
    """
    x.reverse()
    lc_inv = finitefield.FinitePrimeFieldElement(x[0], p).inverse()
    coeff = []
    for c in x[1:]:
        coeff.append((c * lc_inv).n)
    sol = []
    for i in bigrange.range(p):
        if (i**3 + coeff[0]*i**2 + coeff[1]*i + coeff[2]) % p == 0:
            sol.append(i)
            break
    if len(sol) == 0:
        return sol
    X = e2_Fp([coeff[1] + (coeff[0] + sol[0])*sol[0], coeff[0] + sol[0], 1], p)
    if len(X) != 0:
        sol.extend(X)
    return sol

def Newton(f, initial=1, repeat=250):
    """
    Compute x s.t. 0 = f[0] + f[1] * x + ... + f[n] * x ** n
    """
    length = len(f)
    df = []
    for i in range(1, length):
        df.append(i * f[i])
    l = initial
    for k in range(repeat):
        coeff = Float(0)
        dfcoeff = Float(0)
        for i in range(1, length):
            coeff = coeff * l + f[-i]
            dfcoeff += dfcoeff * l + df[-i]
        coeff = coeff * l + f[0]
        if coeff == 0:
            return l
        elif dfcoeff == 0: # Note coeff != 0
            raise ValueError("There is not solution or Choose different initial")
        else:
            l = l - coeff / dfcoeff
    return l


def SimMethod(f, NewtonInitial=1, repeat=250):
    """
    Return zeros of a polynomial given as a list.
    """
    if NewtonInitial != 1:
        ni = NewtonInitial
    else:
        ni = None
    return _SimMethod(f, newtoninitial=ni, repeat=repeat)


def _SimMethod(g, initials=None, newtoninitial=None, repeat=250):
    """
    Return zeros of a polynomial given as a list.

    - g is the list of the polynomial coefficient in ascending order.
    - initial (optional) is a list of initial approximations of zeros.
    - newtoninitial (optional) is an initial value for Newton method to
      obtain an initial approximations of zeros if 'initial' is not given.
    - repeat (optional) is the number of iteration. The default is 250.
    """
    if initials is None:
        if newtoninitial:
            z = _initialize(g, newtoninitial)
        else:
            z = _initialize(g)
    else:
        z = initials

    f = uniutil.polynomial(enumerate(g), imaginary.theComplexField)
    deg = f.degree()
    df = f.differentiate()

    value_list = [f(z[i]) for i in range(deg)]
    for loop in range(repeat):
        sigma_list = [0] * deg
        for i in range(deg):
            if not value_list[i]:
                continue
            sigma = 0
            for j in range(i):
                sigma += 1 / (z[i] - z[j])
            for j in range(i+1, deg):
                sigma += 1 / (z[i] - z[j])
            sigma_list[i] = sigma

        for i in range(deg):
            if not value_list[i]:
                continue
            ratio = value_list[i] / df(z[i])
            z[i] -= ratio / (1 - ratio*sigma_list[i])
            value_list[i] = f(z[i])

    return z

def _initialize(g, newtoninitial=None):
    """
    create initial values of equation given as a list g.
    """
    g = [float(c) for c in g[:]]
    q = [-abs(c) for c in g[:-1]]
    q.append(abs(g[-1]))
    if newtoninitial is None:
        r = Newton(q, _upper_bound_of_roots(q))
    else:
        r = Newton(q, newtoninitial)

    deg = len(g) - 1
    center = -g[-2]/(deg*g[-1])
    about_two_pi = 6
    angular_step = cmath.exp(Complex(0, 1) * about_two_pi / deg)
    angular_move = r
    z = []
    for i in range(deg):
        z.append(center + angular_move)
        angular_move *= angular_step

    return z

def _upper_bound_of_roots(g):
    """
    Return an upper bound of roots.
    """
    weight = len(filter(None, g))
    assert g[-1]
    return max([pow(weight*abs(c/g[-1]), 1/len(g)) for c in g])

def root_Fp(g, p, flag=True):
    """
    Return a root over F_p of nonzero polynomial g.
    p must be prime.
    If flag = False, return a root randomly
    """
    if isinstance(g, list):
        if not isinstance(g[0], tuple):
            g = zip(range(len(g)), g)
    Fp = finitefield.FinitePrimeField(p)
    g = uniutil.FinitePrimeFieldPolynomial(g, Fp)
    h = uniutil.FinitePrimeFieldPolynomial({1:-1, p:1}, Fp)
    g = g.gcd(h)
    deg_g = g.degree()
    if g[0] == 0:
        deg_g = deg_g - 1
        g = g.shift_degree_to(deg_g)
    while True:
        if deg_g == 0:
            return None
        if deg_g == 1:
            return (-g[0]/g[1]).toInteger()
        elif deg_g == 2:
            d = g[1]*g[1] - 4*g[0]
            e = arith1.modsqrt(d.toInteger(), p)
            return ((-g[1]-e)/(2*g[2])).toInteger()
        deg_h = 0
        x = uniutil.FinitePrimeFieldPolynomial({0:-1, (p-1)>>1:1}, Fp)
        if flag:
            a = 0
            while deg_h == 0 or deg_h == deg_g:
                b = uniutil.FinitePrimeFieldPolynomial({0:-a, 1:1}, Fp)
                v = g(b)
                h = x.gcd(v)
                a = a + 1
                deg_h = h.degree()
                b = uniutil.FinitePrimeFieldPolynomial({0:a-1, 1:1}, Fp)
        else:
            while deg_h == 0 or deg_h == deg_g:
                a = bigrandom.randrange(p)
                b = uniutil.FinitePrimeFieldPolynomial({0:-a, 1:1}, Fp)
                v = g(b)
                h = x.gcd(v)
                deg_h = h.degree()
                b = uniutil.FinitePrimeFieldPolynomial({0:a, 1:1}, Fp)
        g = h(b)
        deg_g = deg_h

def allroots_Fp(g, p):
    """
    Return roots over F_p of nonzero polynomial g.
    p must be prime.
    """
    if isinstance(g, list):
        if not isinstance(g[0], tuple):
            g = zip(range(len(g)), g)
    Fp = finitefield.FinitePrimeField(p)
    g = uniutil.FinitePrimeFieldPolynomial(g, Fp)
    h = uniutil.FinitePrimeFieldPolynomial({1:-1, p:1}, Fp)
    g = g.gcd(h)
    deg_g = g.degree()
    roots = []
    if g[0] == 0:
        roots.append(0)
        deg_g = deg_g - 1
        g = g.shift_degree_to(deg_g)
    return roots + roots_loop(g, deg_g, p, Fp)

def roots_loop(g, deg_g, p, Fp):
    if deg_g == 0:
        return []
    if deg_g == 1:
        return [(-g[0]/g[1]).toInteger()]
    elif deg_g == 2:
        d = g[1]*g[1]-4*g[0]
        e = arith1.modsqrt(d.toInteger(), p)
        g1 = -g[1]
        g2 = 2*g[2]
        return [((g1 - e) / g2).toInteger(), ((g1 + e) / g2).toInteger()]
    deg_h = 0
    x = uniutil.FinitePrimeFieldPolynomial({0:-1, (p-1)>>1:1}, Fp)
    a = 0
    while deg_h == 0 or deg_h == deg_g:
        b = uniutil.FinitePrimeFieldPolynomial({0:-a, 1:1}, Fp)
        v = g(b)
        h = x.gcd(v)
        a = a + 1
        deg_h = h.degree()
    b = uniutil.FinitePrimeFieldPolynomial({0:a-1, 1:1}, Fp)
    s = h(b)
    deg_s = deg_h
    t = g.exact_division(s)
    deg_t = t.degree()
    return roots_loop(s, deg_s, p, Fp) + roots_loop(t, deg_t, p, Fp)
