"""
hensel -- Hensel Lift classes for factorization
          of integer coefficient polynomials


If you need lifted factors only once, you may want to use lift_upto
function.  On the other hand, if you might happen to need another
lift, you would directly use one of the lifters and its lift method
for consecutive lifts.
"""

from __future__ import division
import sys
import nzmath.arith1 as arith1
import nzmath.rational as rational
import nzmath.poly.ring as polyring


# module globals
the_ring = polyring.PolynomialRing(rational.theIntegerRing)
the_one = the_ring.one
the_zero = the_ring.zero


def _extgcdp(f, g, p):
    """
    _extgcdp(f,g,p) -> u,v,w

    Find u,v,w such that f*u + g*v = w = gcd(f,g) mod p.
    p should be a prime number.

    This is a private function.
    """
    modp = lambda c: c % p
    u, v, w, x, y, z = the_one, the_zero, f, the_zero, the_one, g
    while z:
        zlc = z.leading_coefficient()
        if zlc != 1:
            normalizer = arith1.inverse(zlc, p)
            x = (x * normalizer).coefficients_map(modp)
            y = (y * normalizer).coefficients_map(modp)
            z = (z * normalizer).coefficients_map(modp)
        q = w.pseudo_floordiv(z)
        u, v, w, x, y, z = x, y, z, u - q*x, v - q*y, w - q*z
        x = x.coefficients_map(modp)
        y = y.coefficients_map(modp)
        z = z.coefficients_map(modp)
    if w.degree() == 0 and w[0] != 1:
        u = u.scalar_exact_division(w[0]) # u / w
        v = v.scalar_exact_division(w[0]) # v / w
        w = w.getRing().one # w / w
    return u, v, w

def _init_cofactors(factors):
    """
    Return list of bi's; bi = a{i+1}*...*ar where ai's are of factors.

    This is a private function.
    """
    bis = [the_one]
    for ai in factors[-1:0:-1]:
        bis.append(bis[-1] * ai)
    bis.reverse()
    return bis


class HenselLiftPair(object):
    """
    A class represents integer polynomial pair which will be lifted by
    Hensel's method.
    """
    def __init__(self, target, factor1, factor2, ladder1, ladder2, p, q=None):
        """
        HenselLiftPair(f, a1, a2, u1, u2, p, q)

        The parameters satisfy that:
          a1 and a2 are monic,
          f == a1*a2 (mod q) and
          a1*u1 + a2*u2 == 1 (mod p),
        with positive integers p dividing q.
        
        If p==q, q can be omit from the argument.
        """
        self.f = target
        self.a1 = factor1
        self.a2 = factor2
        self.u1 = ladder1
        self.u2 = ladder2
        self.p = p
        if q is None:
            self.q = p
        else:
            self.q = q

    @classmethod
    def from_factors(cls, target, factor1, factor2, p):
        """
        Create and return an instance of HenselLiftPair.

        HenselLiftPair.from_factors(f, a1, a2, p)

        The parameters satisfy that:
          f == a1*a2 (mod p)
        with a prime number p.
        """
        u1, u2 = _extgcdp(factor1, factor2, p)[:2]
        return cls(target, factor1, factor2, u1, u2, p)

    def lift_factors(self):
        """
        Update factors by lifted integer coefficient polynomials Ai's:
          f == A1*A2 (mod p*q) and
          Ai == ai (mod q) (i = 1, 2).
        Moreover, q is updated by p*q

        PRECONDITIONS (automatically satisfied):
          f == a1*a1 (mod q)
          a1*u1 + a2*u2 == 1 (mod p)
          p | q
        """
        assert not (self.a1*self.u1 + self.a2*self.u2 - the_one).coefficients_map(lambda c: c % self.p)
        higher = (self.f - self.a1 * self.a2).scalar_exact_division(self.q)
        quot, rem = (self.u2 * higher).pseudo_divmod(self.a1)
        filterq = lambda c: (c % self.p) * self.q
        self.a1 += rem.coefficients_map(filterq)
        self.a2 += (self.u1 * higher + self.a2 * quot).coefficients_map(filterq)
        self.q *= self.p
        assert not (self.a1*self.u1 + self.a2*self.u2 - the_one).coefficients_map(lambda c: c % self.p)

    def lift_ladder(self):
        """
        Update u1 and u2 with U1 and U2:
          a1*U1 + a2*U2 == 1 (mod p**2)
          Ui == ui (mod p) (i = 1, 2)
        then, update p with p**2.

        PRECONDITIONS (automatically satisfied):
          a1*u1 + a2*u2 == 1 (mod p)
        """
        higher = (the_one - self.u1 * self.a1 - self.u2 * self.a2).scalar_exact_division(self.p)
        quot, rem = (self.u2 * higher).pseudo_divmod(self.a1)
        filterp = lambda c: (c % self.p) * self.p
        self.u1 += (self.u1 * higher + self.a2 * quot).coefficients_map(filterp)
        self.u2 += rem.coefficients_map(filterp)
        self.p *= self.p

    def lift(self):
        """
        The lift.
        """
        self.lift_factors()
        self.lift_ladder()

    def _get_factors(self):
        """
        getter for factors of target mod q
        """
        return [self.a1, self.a2]

    factors = property(_get_factors)


class HenselLiftMulti(object):
    """
    A class represents integer polynomials which will be lifted by
    Hensel's lemma.

    If the number of factors is two, you should use HenselLiftPair.
    """
    def __init__(self, target, factors, cofactors, ladder, p, q=None):
        """
        HenselLiftMulti(f, factors, ladder, p, q)

        The parameters satisfy that:
         +  ai's of factors and f are all monic
         +  bi's of cofactors are product of aj's for i < j
         +  f == a1*...*ar (mod q)
         +  ladder is a tuple (sis, tis) and si's of sis and ti's of tis
            satisfy: ai*si + bi*ti == 1 (mod p),
        with positive integers p dividing q.

        If p==q, q can be omit from the argument.
        """
        self.f = target
        self.a = list(factors)
        self.r = len(self.a)
        self.b = list(cofactors)
        self.s = list(ladder[0])
        self.t = list(ladder[1])
        self.p = p
        if q is None:
            self.q = p
        else:
            self.q = q

    @classmethod
    def from_factors(cls, target, factors, p):
        """
        Create and return an instance of HenselLiftMulti.

        HenselLiftMulti.from_factors(f, factors, p)

        The parameters satisfy that:
          f == a1*...*ar (mod p)
        with a prime number p and ai's in factors.
        """
        bi0s = _init_cofactors(factors)
        sis, tis = [], []
        for ai0, bi0 in zip(factors, bi0s):
            s, t = _extgcdp(ai0, bi0, p)[:2]
            sis.append(s)
            tis.append(t)
        return cls(target, factors, bi0s, (sis, tis), p, p)

    def lift_factors(self):
        """
        Find a lifted integer coefficient polynomials such that:
          f = A1*A2*...*Ar (mod q*p),
          Ai = ai (mod q) (i=1,2,...,r),
        then, update q with q*p.

        PRECONDITIONS (automatically satisfied):
          f = a1*a2*...*ar (mod q),
          ai*si + bi*ti = 1 (mod p) (i=1,2,...,r) and
          p | q.
        """
        mini_target = self.f
        for i in range(self.r - 1):
            assert not (self.a[i]*self.s[i] + self.b[i]*self.t[i] - the_one).coefficients_map(lambda c: c % self.p), \
                   "a = %s\ns = %s\nb = %s\nt = %s" %(self.a[i], self.s[i], self.b[i], self.t[i])
            mini_pair = HenselLiftPair(mini_target,
                                       self.a[i], self.b[i],
                                       self.s[i], self.t[i],
                                       self.p, self.q)
            mini_pair.lift_factors()
            self.a[i] = mini_pair.a1
            mini_target = self.b[i] = mini_pair.a2
        self.a[self.r - 1] = mini_target
        self.q = mini_pair.q

    def lift_ladder(self):
        """
        Update ladder si's and ti's with Si's and Ti's:
          ai*Si + bi*Ti == 1 (mod p**2) (i = 1, 2, ..., r),
          Si == si (mod p) (i = 1, 2, ..., r) and
          Ti == ti (mod p) (i = 1, 2, ..., r)
        then, update p with p**2.

        PRECONDITIONS (automatically satisfied):
          ai*si + bi*ti == 1 (mod p)
        """
        for i in range(self.r - 1):
            mini_pair = HenselLiftPair(self.f,
                                       self.a[i], self.b[i],
                                       self.s[i], self.t[i],
                                       self.p, self.q)
            mini_pair.lift_ladder()
            self.s[i], self.t[i] = mini_pair.u1, mini_pair.u2
        self.p *= self.p

    def lift(self):
        """
        The lift.
        """
        self.lift_factors()
        self.lift_ladder()

    def _get_factors(self):
        """
        getter for factors of target mod q
        """
        return list(self.a)

    factors = property(_get_factors)


class HenselLiftSimultaneously(object):
    """
    INVARIANTS:
      ai's, pi's and gi's are monic
      f == g1*g2*...*gr (mod p)
      f == d0 + d1*p + d2*p**2 +...+ dk*p**k
      hi == g(i+1)*...*gr
      1 == gi*si + hi*ti (mod p) (i = 1, 2,..., r)
      deg(si) < deg(hi), deg(ti) < deg(gi)
      p | q
      f == l1*l2*...*lr (mod q/p)
      f == a1*a2*...*ar (mod q)
      ui == ai*yi + bi*zi (mod p) (i = 1, 2,..., r)

    REFERENCE:
    G.E.Collins & M.J.Encarnaci'on.
    Improved Techniques for factoring Univariate Polynomials,
    J.Symb.Comp. 21, 313--327 (1996).
    """

    def __init__(self, target, factors, cofactors, bases, p):
        """
        Prepare attributes
        """
        # invariants
        self.f = target
        self.gis = tuple(factors)
        self.his = tuple(cofactors)
        self.r = len(self.gis)
        self.sis = tuple(bases[0])
        self.tis = tuple(bases[1])
        self.p = p
        self.dis = []

        # q will be a power of p
        self.q = self.p

        # these are containers
        self.ais, self.bis = {}, {} # used later
        self.uis = self.yis = self.zis = None # used later

        # functions
        self.modp = lambda c: c % self.p

    def _init_dis(self):
        """
        Return p-adic expansion of target polynomial.

        This method is private for __init__.
        """
        self.dis.append(self.f.coefficients_map(self.modp))
        rest = self.f - self.dis[-1]
        q = self.p
        while any(abs(x) >= q for x in self.f.itercoefficients()):
            self.dis.append(self.f.coefficients_map(lambda c: (c // q) % self.p))
            rest -= self.dis[-1].scalar_mul(q)
            q *= self.p
        if self.f != rest:
            self.dis.append(rest.scalar_exact_division(q))

    @classmethod
    def from_factors(cls, target, factors, p, ubound=sys.maxint):
        """
        Create and return an instance of HenselLiftSimultaneously,
        whose factors are lifted by HenselLiftMulti upto ubound (if it
        is smaller than sys.maxint, or upto sys.maxint).

        HenselLiftSimultaneously.from_factors(f, factors, p)

        The parameters satisfy that:
          f == a1*...*ar (mod p)
        with a prime number p and ai's in factors.
        """
        lifter = HenselLiftMulti.from_factors(target, factors, p)
        interbound = min(ubound, sys.maxint)
        while lifter.p < interbound:
            lifter.lift_factors()
            lifter.lift_ladder()
        new = cls(target, lifter.a, lifter.b, (lifter.s, lifter.t), lifter.p)
        new.first_lift()
        return new

    def first_lift(self):
        """
        Start lifting.
            f == l1*l2*...*lr (mod p**2)
        Initialize di's, ui's, yi's and zi's.
        Update ai's, bi's.
        Then, update q with p**2.
        """
        assert self.p == self.q # q has not been raised yet
        self._assertEqualModulo(self.f, arith1.product(self.gis), self.q)
        for i in range(self.r - 1):
            self._assertEqualModulo(self.f, arith1.product(self.gis[:i+1])*self.his[i], self.q)
        assert self.gis[-1] == self.his[self.r - 2]

        self._init_dis()
        if len(self.dis) > 1:
            mini_target = self.dis[0] + self.dis[1].scalar_mul(self.p)
            self._assertEqualModulo(self.f, mini_target, self.p**2)
        else:
            mini_target = self.dis[0]
        aj, bj = [], []
        self.uis, self.yis, self.zis = [], {}, {}
        for i in xrange(self.r - 1):
            dividend = mini_target - self.gis[i] * self.his[i]
            self.uis.append(dividend.scalar_exact_division(self.p))
            self.yis[i], self.zis[i] = self._solve_yz(i)
            aj.append(self.gis[i] + self.zis[i].scalar_mul(self.q))
            bj.append(self.his[i] + self.yis[i].scalar_mul(self.q))
            self._assertEqualModulo(mini_target, aj[i]*bj[i], self.q*self.p)
            mini_target = bj[i]
        self._assertEqualModulo(self.gis[-1], mini_target, self.q)
        aj.append(bj[-1])
        self.q *= self.p

        for i in range(self.r - 1):
            self._assertEqualModulo(self.f, arith1.product(aj[:i+1])*bj[i], self.q, "f != l%d * m%d" % (i, i))
            self._assertEqualModulo(self.gis[i], aj[i], self.p)
        self._assertEqualModulo(self.f, arith1.product(aj), self.q)

        self.ais[-1], self.ais[0] = self.gis, tuple(aj)
        self.bis[-1], self.bis[0] = self.his, tuple(bj)

    def general_lift(self):
        """
        Continue lifting.
            f == a1*a2*...*ar (mod p*q)
        Update ai's, bi's, ui's, yi's and zi's.
        Then, update q with p*q.
        """
        j = arith1.vp(self.q, self.p)[0]
        if len(self.dis) > j:
            mini_target = self.dis[j]
        else:
            mini_target = the_zero

        self.ais[-2], self.ais[-1] = self.ais[-1], self.ais[0]
        self.bis[-2], self.bis[-1] = self.bis[-1], self.bis[0]
        aj, bj = [], []
        for i in xrange(self.r - 1):
            yi, zi = self.yis[i], self.zis[i]
            dividend = self.ais[-2][i]*yi + self.bis[-2][i]*zi - self.uis[i]
            v_j = dividend.scalar_exact_division(self.p)
            if j == 2:
                self.uis[i] = mini_target - v_j - yi*zi
            else:
                self.uis[i] = mini_target - v_j - (yi*zi).scalar_mul(self.p**(j - 2))
            self.yis[i], self.zis[i] = self._solve_yz(i)
            aj.append(self.ais[-1][i] + self.zis[i].scalar_mul(self.q))
            bj.append(self.bis[-1][i] + self.yis[i].scalar_mul(self.q))
            self._assertEqualModulo(self.f, arith1.product(aj)*bj[i], self.q*self.p, (i, j))
            mini_target = self.yis[i]
        aj.append(bj[-1])
        self.q *= self.p

        for i in range(self.r - 1):
            self._assertEqualModulo(self.f, arith1.product(aj[:i+1])*bj[i], self.q, "f != a%d * b%d" % (i, i))
            self._assertEqualModulo(self.gis[i], aj[i], self.p)
        self._assertEqualModulo(self.f, arith1.product(aj), self.q)

        self.ais[0] = tuple(aj)
        self.bis[0] = tuple(bj)

    def lift(self):
        """
        The lift
        """
        self.general_lift()

    def _solve_yz(self, i):
        """
        Solve the equation
            gi Y + hi Z = ui (mod p)
        satisfying conditions:
        (1) deg(Y) < deg(hi) and
        (2) deg(Z) < deg(gi),
        and return a tuple(Y, Z).
        The method needs the following conditions:
        (I) deg(ui) <= deg(gi)deg(hi),
        (II)  gi si + hi ti = 1 (mod p).
        """
        umodp = self.uis[i].coefficients_map(self.modp)
        quot, rem = (self.tis[i] * umodp).pseudo_divmod(self.gis[i])
        y = (self.sis[i] * umodp + self.his[i] * quot).coefficients_map(self.modp)
        z = rem.coefficients_map(self.modp)

        assert y.degree() < self.his[i].degree()
        self._assertEqualModulo(self.gis[i] * y + self.his[i] * z,
                                self.uis[i], self.p)

        return y, z

    def _get_factors(self):
        """
        getter for factors of target mod q
        """
        if self.ais:
            return list(self.ais[0])
        else:
            return list(self.gis)

    factors = property(_get_factors)

    def _assertEqualModulo(self, expected, actual, modulus, message=None):
        """
        assert expected == actual (mod modulus)
        """
        if message is None:
            for d, c in actual - expected:
                assert 0 == c % modulus, "degree %d\nactual = %s" % (d, str(actual))
        else:
            for d, c in actual - expected:
                assert 0 == c % modulus, "degree %d\nactual = %s\n%s" % (d, str(actual), message)


def lift_upto(target, factors, p, bound):
    """
    Hensel lift factors mod p of target upto bound and return factors
    mod q and q.

    PRECONDITIONS:
    target == product(factors) mod p
    POSTCONDITIONS:
    there exist k s.t. q == p**k >= bound and
    target == product(result) mod q
    """
    if len(factors) == 2:
        lifter = HenselLiftPair.from_factors(target, factors[0], factors[1], p)
    elif bound < sys.maxint:
        lifter = HenselLiftMulti.from_factors(target, factors, p)
    else:
        lifter = HenselLiftSimultaneously.from_factors(target, factors, p)
    while lifter.q < bound:
        lifter.lift()
    return lifter.factors, lifter.q
