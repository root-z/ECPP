""" Elliptic Curves over finite field.
"""

from __future__ import division
import logging
import random

import nzmath.poly.uniutil as uniutil
import nzmath.poly.multiutil as multiutil
import nzmath.poly.termorder as termorder
import nzmath.arith1 as arith1
import nzmath.bigrandom as bigrandom
import nzmath.gcd as gcd
import nzmath.prime as prime
import nzmath.rational as rational
import nzmath.ring as ring
import nzmath.factor.methods as factor_methods
import nzmath.factor.misc as factor_misc
import nzmath.finitefield as finitefield
import nzmath.compatibility

_log = logging.getLogger('nzmath.elliptic')

# symbol aliases
uniutil.special_ring_table[finitefield.FinitePrimeField] = uniutil.FinitePrimeFieldPolynomial
MultiVarPolynomial = multiutil.MultiVariableSparsePolynomial

# polynomial wrapper
def _UniVarPolynomial(dict,coeffring=None):
    return uniutil.polynomial(dict,coeffring)

# string format wrapper
def _strUniPoly(poly, symbol="x", asc=True):
    """return format string of UniutilPolynomial"""
    return termorder.ascending_order.format(poly, symbol, asc)

def _strMultiPoly(poly, symbol=["x","y"], asc=False):
    """return format string of MultiVarPolynomial for EC"""
    def Weierstrasscmp(left, right):
        sum_left, sum_right = sum(left), sum(right)
        if sum_left != sum_right:
            return -cmp(sum_left, sum_right)
        return -cmp(right, left)

    return termorder.MultivarTermOrder(Weierstrasscmp).format(MultiVarPolynomial(poly,symbol), symbol, asc)


def _PolyMod(f, g):
    """
    return f (mod g)
    """
    return f % g

def _PolyGCD(f, g):
    # other cases
    return f.gcd(g)

def _PolyPow(f, d, g):
    """
    returns (f^d)%g
    """
    return g.mod_pow(f, d)

def _PolyMulRed(multipliees, poly):
    """
    multipliees[*] is (OneSparsePoly,int,long)
    poly is OneSparsePoly
    """
    if poly.degree() < 1:
        return poly.getRing().zero
    product = multipliees.pop()
    for factor in multipliees:
        #if factor.degree() >= poly.degree():
        #factor = PolyMod(factor, poly)
        #if factor == 0:
        #    return 0
        product = product * factor
        if product.degree() >= poly.degree():
            product = _PolyMod(product, poly)
            if not product:
                break
    return product


class ECGeneric:
    """
    Definition of Elliptic curves over generic field.
    this class is fundamental class, normally only called sub class.
    """
    def __init__(self, coefficient, basefield=None):
        """
        Initialize an elliptic curve. If coefficient has 5 elements,
        it represents E:y**2+a1*x*y+a3*y=x**3+a2*x**2+a4*x+a6 or 2
        elements, E:y*2=x*3+a*x+b.
        """

        try:
            character = basefield.getCharacteristic()
            self.basefield = basefield
        except:
            # backward compatibility support
            if isinstance(basefield, rational.RationalField) or (not basefield):
                character = 0
                self.basefield = rational.theRationalField
            elif isinstance(basefield, (int,long)):
                character = basefield
                if character == 1 or character < 0:
                    raise ValueError("basefield characteristic must be 0 or prime.")
                self.basefield = finitefield.FinitePrimeField.getInstance(character)
            else:
                raise ValueError("basefield must be FiniteField.")

        self.ch = character
        self.infpoint = [self.basefield.zero]
        if isinstance(coefficient, list):
            self.coefficient = coefficient
            if self.ch == 0:
                if len(self) == 5:
                    self.a1 = self.coefficient[0]
                    self.a2 = self.coefficient[1]
                    self.a3 = self.coefficient[2]
                    self.a4 = self.coefficient[3]
                    self.a6 = self.coefficient[4]
                    self.b2 = self.a1**2+4*self.a2
                    self.b4 = self.a1*self.a3+2*self.a4
                    self.b6 = self.a3**2+4*self.a6
                    self.b8 = self.a1**2*self.a6+4*self.a2*self.a6-self.a1*self.a3*self.a4+self.a2*self.a3**2-self.a4**2
                    self.c4 = self.b2**2-24*self.b4
                    self.c6 = -self.b2**3+36*self.b2*self.b4-216*self.b6
                    self.disc = -self.b2**2*self.b8-8*self.b4**3-27*self.b6**2+9*self.b2*self.b4*self.b6
                elif len(self) == 2:
                    self.a = self.coefficient[0]
                    self.b = self.coefficient[1]
                    self.a1 = 0
                    self.a2 = 0
                    self.a3 = 0
                    self.a4 = self.coefficient[0]
                    self.a6 = self.coefficient[1]
                    self.b2 = 0
                    self.b4 = 2*self.a
                    self.b6 = 4*self.b
                    self.b8 = -self.a**2
                    self.c4 = -48*self.a
                    self.c6 = -864*self.b
                    self.disc = (self.c4**3-self.c6**2)/1728
                else:
                    raise ValueError("coefficient is less or more, can't defined EC.")
                if self.disc == 0:
                    raise ValueError("this curve is singular.")
                self.j = (self.c4**3)/self.disc
                self.cubic = _UniVarPolynomial({0:self.a6, 1:self.a4,
                                               3:self.basefield.one},
                                              self.basefield)
            else:
                pass # support for subclass
        else:
            raise ValueError("parameters must be (coefficient, basefield)")
    def __len__(self):
        return len(self.coefficient)

    def __repr__(self):
        """
        return represantation form.
        only defined generic representation.
        """
        return "EC(%s, %s)" % (str(self.coefficient), repr(self.basefield))

    def __str__(self):
        dictL = {(0, 2):self.basefield.one}
        dictR = {3:self.basefield.one}
        if self.a1:
            dictL[(1, 1)] = self.a1
        if self.a2:
            dictR[2] = self.a2
        if self.a3:
            dictL[(0, 1)] = self.a3
        if self.a4:
            dictR[1] = self.a4
        if self.a6:
            dictR[0] = self.a6
        return ("%s = %s" %
                (str(_strMultiPoly(dictL)),
                 str(_strUniPoly(_UniVarPolynomial(dictR, self.basefield)))))

    def simple(self):
        """
        Transform
          E:y^2+a1*x*y+a3*y = x^3+a2*x^2+a4*x+a6
        to
          E':Y^2 = X^3+(-27*c4)*X+(-54*c6),
        if ch is not 2 or 3
        """
        if len(self) == 2 or not self.a1 and not self.a2 and not self.a3:
            return self
        elif self.ch == 2 or self.ch == 3:
            return self
        else:
            return EC([-27*self.c4, -54*self.c6], self.basefield)

    def changeCurve(self, V):
        """
        this transforms E to E' by V=[u,r,s,t]
          x->u^2*x'+r,y->u^3*y'+s*u^2*x'+t
        """
        # FIXME: is "(not V[0] != 0)" a correct condition?
        if (not isinstance(V, list)) or (not len(V) == 4) or (not V[0] != 0):
            raise ValueError("([u, r, s, t]) with u must not be 0.")

        # FIXME: this is dependent rational.Rational
        if self.ch == 0:
            return EC([rational.Rational(self.a1+2*V[2], V[0]),
                       rational.Rational(self.a2-V[2]*self.a1+3*V[1]-V[2]**2, V[0]**2),
                       rational.Rational(self.a3+V[1]*self.a1+2*V[3], V[0]**3),
                       rational.Rational(self.a4-V[2]*self.a3+2*V[1]*self.a2-(V[3]+V[1]*V[2])*self.a1+3*V[1]**2-2*V[2]*V[3], V[0]**4),
                       rational.Rational(self.a6+V[1]*self.a4+V[1]**2*self.a2+V[1]**3-V[3]*self.a3-V[3]**2-V[1]*V[3]*self.a1, V[0]**6)])
        else:
            for v in V:
                if not isinstance(v, (int, long)) and not (v in self.basefield):
                    raise ValueError("transform V must be integer sequence.")
            v = self.basefield.createElement(V[0]).inverse()
            return EC([(self.a1+2*V[2])*v,
                       ((self.a2-V[2]*self.a1+3*V[1]-V[2]**2)*v**2),
                       ((self.a3+V[1]*self.a1+2*V[3])*v**3),
                       ((self.a4-V[2]*self.a3+2*V[1]*self.a2-(V[3]+V[1]*V[2])*self.a1+3*V[1]**2-2*V[2]*V[3])*v**4),
                       ((self.a6+V[1]*self.a4+V[1]**2*self.a2+V[1]**3-V[3]*self.a3-V[3]**2-V[1]*V[3]*self.a1)*v**6)], self.basefield)

    def changePoint(self, P, V):
        """
        this transforms P to P' by V=[u,r,s,t] ; x->u^2*x'+r,y->u^3*y'+s*u^2*x'+t
        """
        if (not (isinstance(P, list) and isinstance(V, list))) or \
               not (len(P) == 2 and len(V) == 4 and V[0] != 0):
            raise ValueError("(P,V) must be ([px, py], [u, r, s, t]) with u != 0.")

        if self.ch == 0:
            Q0 = rational.IntegerIfIntOrLong(P[0]-V[1])/rational.IntegerIfIntOrLong(V[0]**2)
            Q1 = rational.IntegerIfIntOrLong(P[1]-V[2]*(P[0]-V[1])-V[3])/rational.IntegerIfIntOrLong(V[0]**3)
        else:
            v = self.basefield.createElement(V[0]).inverse()
            Q0 = ((P[0]-V[1])*v**2)
            Q1 = ((P[1]-V[2]*(P[0]-V[1])-V[3])*v**3)
        Q = [Q0, Q1]
        return Q

    def coordinateY(self, x):
        """
        this returns the y(P)>0,x(P)==x
        """
        if self.ch > 2:
            y1 = self.a1*x + self.a3
            y2 = x**3 + self.a2*x**2 + self.a4*x + self.a6
            if len(self) == 2 or (self.a1 == self.a2 == self.a3 == 0):
                if self.basefield.Legendre(y2) == 1:
                    return self.basefield.sqrt(y2)
                else:
                    return False
            else:
                if y1**2 + 4*y2 >= 0:
                    d = y1**2 + 4*y2
                    return (-y1 - self.basefield.sqrt(d)) // (2*self.basefield.one)
                else:
                    return False
        elif self.ch == 2:
            raise NotImplementedError("This is not implemented.")
        y1 = self.a1*x + self.a3
        y2 = x**3 + self.a2*x**2 + self.a4*x + self.a6
        Y = y1**2 + 4*y2
        if Y >= 0:
            if isinstance(Y, rational.Rational):
                yn = arith1.issquare(Y.numerator)
                yd = arith1.issquare(Y.denominator)
                if yn and yd:
                    return ((-1)*y1 + rational.Rational(yn, yd)) / 2
            else:
                Z = arith1.issquare(Y)
                if Z:
                    return rational.Rational((-1)*y1 + Z, 2)
        return False

    def whetherOn(self, P):
        """
        Determine whether P is on curve or not.
        Return True if P is on, return False otherwise.
        """
        if isinstance(P, list):
            if len(P) == 2:
                if self.ch == 0:
                    if P[1]**2+self.a1*P[0]*P[1]+self.a3*P[1] == P[0]**3+self.a2*P[0]**2+self.a4*P[0]+self.a6:
                        return True
                    else:
                        #_log.debug(str(P[1]**2+self.a1*P[0]*P[1]+self.a3*P[1]))
                        #_log.debug(str(P[0]**3+self.a2*P[0]**2+self.a4*P[0]+self.a6))
                        return False
                else:
                    if P[1]**2+self.a1*P[0]*P[1]+self.a3*P[1] == P[0]**3+self.a2*P[0]**2+self.a4*P[0]+self.a6:
                        return True
                    else:
                        return False
            elif P == [self.basefield.zero]:
                return True
        raise ValueError("point P must be [px, py] or [0].")

    def add(self, P, Q):
        """
        return point addition P+Q
        """
        if not (isinstance(P, list) and isinstance(Q, list)):
            raise ValueError("point P (resp. Q) must be [px, py] (resp. [qx, qy])")
        #if not (self.whetherOn(P) and self.whetherOn(Q)):
        #    raise ValueError("either points must not be point on curve.")

        if (P == self.infpoint) and (Q != self.infpoint):
            return Q
        elif (P != self.infpoint) and (Q == self.infpoint):
            return P
        elif (P == self.infpoint) and (Q == self.infpoint):
            return self.infpoint

        if self.ch == 0:
            # FIXME
            if P[0] == Q[0]:
                if P[1]+Q[1]+self.a1*Q[0]+self.a3 == 0:
                    return self.infpoint
                else:
                    s = (3*P[0]**2+2*self.a2*P[0]+self.a4-self.a1*P[1])/(2*P[1]+self.a1*P[0]+self.a3)
                    t = (-P[0]**3+self.a4*P[0]+2*self.a6-self.a3*P[1])/(2*P[1]+self.a1*P[0]+self.a3)
            else:
                s = (Q[1]-P[1])/(Q[0]-P[0])
                t = (P[1]*Q[0]-Q[1]*P[0])/(Q[0]-P[0])
            x3 = s**2+self.a1*s-self.a2-P[0]-Q[0]
            y3 = -(s+self.a1)*x3-t-self.a3
            R = [x3, y3]
            return R
        else:
            if not (P[0] - Q[0]):
                # FIXME: the condition is P[0] == Q[0] intuitively,
                #        but sometimes there are int vs FiniteFieldElement
                #        comparisons ...
                if not (P[1]+Q[1]+self.a1*Q[0]+self.a3):
                    return self.infpoint
                else:
                    s = (3*P[0]**2+2*self.a2*P[0]+self.a4-self.a1*P[1])/(2*P[1]+self.a1*P[0]+self.a3)
                    t = (-P[0]**3+self.a4*P[0]+2*self.a6-self.a3*P[1])/(2*P[1]+self.a1*P[0]+self.a3)
            else:
                s = (Q[1] - P[1]*self.basefield.one) / (Q[0] - P[0])
                t = (P[1]*Q[0] - Q[1]*P[0]*self.basefield.one)/ (Q[0] - P[0])
            x3 = s**2+self.a1*s-self.a2-P[0]-Q[0]
            y3 = -(s+self.a1)*x3-t-self.a3
            R = [x3, y3]
            return R

    def sub(self, P, Q):
        """
        return point subtraction P-Q
        """
        if not (isinstance(P, list) and isinstance(Q, list)):
            raise ValueError("point P (resp. Q) must be [px, py] (resp. [qx, qy])")
        #if not (self.whetherOn(P) and self.whetherOn(Q)):
        #    raise ValueError("either points must not be point on curve.")

        if (P != self.infpoint) and (Q == self.infpoint):
            return P
        elif (P == self.infpoint) and (Q == self.infpoint):
            return self.infpoint

        x = Q[0]
        y = -Q[1]-self.a1*Q[0]-self.a3
        R = [x, y]
        if (P == self.infpoint) and (Q != self.infpoint):
            return R
        else:
            return self.add(P, R)

    def mul(self, k, P):
        """
        this returns [k]*P
        """
        if k >= 0:
            l = arith1.expand(k, 2)
            Q = self.infpoint
            for j in range(len(l)-1, -1, -1):
                Q = self.add(Q, Q)
                if l[j] == 1:
                    Q = self.add(Q, P)
            return Q
        else:
            l = arith1.expand(-k, 2)
            Q = self.infpoint
            for j in range(len(l)-1, -1, -1):
                Q = self.add(Q, Q)
                if l[j] == 1:
                    Q = self.add(Q, P)
            return self.sub(self.infpoint, Q)

    def divPoly(self, Number=None):
        """ Return division polynomial
        """
        def heart(q):
            """
            this is for Schoof's, internal function
            """
            l = []
            j = 1
            bound = 4 * arith1.floorsqrt(q)
            for p in prime.generator():
                if p != q:
                    l.append(p)
                    j *= p
                    if j > bound:
                        break
            return l

        # divPoly mainblock
        one = self.basefield.one
        Kx = ring.getRing(self.cubic)
        f = {-1: -Kx.one, 0: Kx.zero, 1: Kx.one, 2:Kx.one}

        # initialize f[3], f[4] and e
        if not Number:
            if self.ch <= 3:
                raise ValueError("You must input (Number)")
            H = heart(card(self.basefield))
            loopbound = H[-1] + 2
            E = self.simple()
            e = 4 * E.cubic
            f[3] = _UniVarPolynomial({4:3*one,
                                     2:6*E.a,
                                     1:12*E.b,
                                     0:-E.a**2},
                                    self.basefield)
            f[4] = _UniVarPolynomial({6:2*one,
                                     4:10*E.a,
                                     3:40*E.b,
                                     2:-10*E.a**2,
                                     1:-8*E.a*E.b,
                                     0:-2*(E.a**3 + 8*E.b**2)},
                                    self.basefield)
        else:
            loopbound = Number + 1
            if self.ch == 0:
                E = self.simple()
                e = E.cubic
                f[3] = _UniVarPolynomial({4:3*one,
                                         2:6*E.a,
                                         1:12*E.b,
                                         0:-E.a**2},
                                        self.basefield)
                f[4] = _UniVarPolynomial({6:2*one,
                                         4:10*E.a,
                                         3:40*E.b,
                                         2:-10*E.a**2,
                                         1:-8*E.a*E.b,
                                         0:-2*(E.a**3 + 8*E.b**2)},
                                        self.basefield)
            else:
                e = _UniVarPolynomial({3:4*one,
                                      2:self.b2,
                                      1:2*self.b4,
                                      0:self.b6},
                                     self.basefield)
                f[3] = _UniVarPolynomial({4:3*one,
                                         3:self.b2,
                                         2:3*self.b4,
                                         1:3*self.b6,
                                         0:self.b8},
                                        self.basefield)
                f[4] = _UniVarPolynomial({6:2*one,
                                         5:self.b2,
                                         4:5*self.b4,
                                         3:10*self.b6,
                                         2:10*self.b8,
                                         1:self.b2*self.b8 - self.b4*self.b6,
                                         0:self.b4*self.b8 - self.b6**2},
                                        self.basefield)

        # recursive calculation
        i = 5
        while i < loopbound:
            if i & 1:
                j = (i - 1) >> 1
                if j & 1:
                    f[i] = f[j+2]*f[j]**3 - e**2*f[j-1]*f[j+1]**3
                else:
                    f[i] = e**2*f[j+2]*f[j]**3 - f[j-1]*f[j+1]**3
            else:
                j = i >> 1
                f[i] = (f[j+2]*f[j-1]**2 - f[j-2]*f[j+1]**2) * f[j]
            i += 1

        # final result
        if not Number:
            for i in range(2, loopbound, 2):
                f[i] = 2*f[i]
            return (f, H)
        else:
            return f[Number]

class ECoverQ(ECGeneric):
    """
    Elliptic curves over Q.
    """
    def __init__(self, coefficient):
        field = rational.RationalField()
        coeffs_list = []
        if isinstance(coefficient, list):
            for c in coefficient:
                if isinstance(c, (int, long)):
                    coeff = field.createElement(c)
                elif c in field:
                    coeff = c
                else:
                    raise ValueError("coefficient not in basefield.")
                coeffs_list.append(coeff)

        ECGeneric.__init__(self, coeffs_list, field)

    def __repr__(self):
        if len(self) == 2 or self.a1 == self.a2 == self.a3 == 0:
            return "ECoverQ(["+repr(self.a4)+","+repr(self.a6)+"])"
        else:
            return "ECoverQ(["+repr(self.a1)+","+repr(self.a2)+","+repr(self.a3)+","+repr(self.a4)+","+repr(self.a6)+"])"

    def point(self, limit=1000):
        """
        this returns a random point on eliiptic curve over Q.
        limit set maximal find time.
        """
        i = 9
        while i < limit:
            s = random.randrange(1, i)
            t = random.randrange(1, i)
            x = rational.Rational(s, t);
            y = self.coordinateY(x)
            if y != False:
                return [x, y]
            i = i+10
        raise ValueError("Times exceeded for limit.")


class ECoverGF(ECGeneric):
    """
    Elliptic curves over Galois Field.
    """
    def __init__(self, coefficient, basefield=None):
        """ create ECoverGF object.
        coefficient must be length 5 or 2 list(represent as Weierstrass form),
        any coefficient is in basefield.
        basefield must be FiniteField subclass object
        (i.e. FinitePrimeField or FiniteExtendedField object.)
        """

        # parameter parse
        try:
            character = basefield.getCharacteristic()
            field = basefield
        except AttributeError:
            # backward compatibility
            if isinstance(basefield, (int, long)):
                field = finitefield.FinitePrimeField.getInstance(basefield)
                character = basefield
            else:
                raise ValueError("basefield must be FiniteField object.")

        coeffs_list = []
        if isinstance(coefficient, list):
            for c in coefficient:
                if isinstance(c, (int, long)):
                    coeff = field.createElement(c)
                elif c in field:
                    coeff = c
                else:
                    raise ValueError("coefficient not in basefield.")
                coeffs_list.append(coeff)

        # general initialize
        ECGeneric.__init__(self, coeffs_list, field)

        zero = self.basefield.zero
        one = self.basefield.one

        # format attribute
        if self.ch == 2:
            if len(self) == 5:
                # FIXME
                if coeffs_list[0] & 1 == one and coeffs_list[2] & 1 == coeffs_list[3] & 1 == zero and coeffs_list[4]:
                    self.a1 = one
                    self.a2 = coeffs_list[1]
                    self.a3 = zero
                    self.a4 = zero
                    self.a6 = coeffs_list[4]
                    self.b2 = one
                    self.b4 = zero
                    self.b6 = zero
                    self.b8 = self.a6
                    self.c4 = one
                    self.c6 = one
                    self.disc = self.a6
                    self.j = self.disc.inverse()
                elif coeffs_list[0] & 1 == coeffs_list[1] & 1 == zero and coeffs_list[2]:
                    self.a1 = zero
                    self.a2 = zero
                    self.a3 = coeffs_list[2]
                    self.a4 = coeffs_list[3]
                    self.a6 = coeffs_list[4]
                    self.b2 = zero
                    self.b4 = zero
                    self.b6 = self.a3**2
                    self.b8 = self.a4**2
                    self.c4 = zero
                    self.c6 = zero
                    self.disc = self.a3**4
                    self.j = zero
                else:
                    raise ValueError("coefficient may be not representation of EC.")
            else:
                raise ValueError("coefficient may only use full Weierstrass form for characteristic 2.")
        elif self.ch == 3: # y^2=x^3+a2*x^2+a6 or y^2=x^3+a4*x+a6
            # FIXME
            if len(self) == 5:
                if coeffs_list[0] % 3 == coeffs_list[2] % 3 == coeffs_list[3] % 3 == 0 and coeffs_list[1] and coeffs_list[4]:
                    self.a1 = zero
                    self.a2 = coeffs_list[1]
                    self.a3 = zero
                    self.a4 = zero
                    self.a6 = coeffs_list[4]
                    self.b2 = self.a2
                    self.b4 = zero
                    self.b6 = self.a6
                    self.b8 = self.a2*self.a6
                    self.c4 = self.b2**2
                    self.c6 = 2*self.b2**3
                    self.disc = -self.a2**3*self.a6
                    self.j = (-self.a2**3)*self.a6.inverse()
                elif coeffs_list[0] == coeffs_list[1] == coeffs_list[2] == 0 and coeffs_list[3]:
                    self.a1 = zero
                    self.a2 = zero
                    self.a3 = zero
                    self.a4 = coeffs_list[3]
                    self.a6 = coeffs_list[4]
                    self.b2 = zero
                    self.b4 = 2*self.a4
                    self.b6 = self.a6
                    self.b8 = 2*self.a4**2
                    self.c4 = zero
                    self.c6 = zero
                    self.disc = -self.a4**3
                    self.j = zero
                else:
                    raise ValueError("can't defined EC.")
                if not self.disc:
                    raise ValueError("this curve is singular.")
            else:
                raise ValueError("coefficient is less or more, can't defined EC.")
        else:
            if len(self) == 5:
                self.a1 = coeffs_list[0]
                self.a2 = coeffs_list[1]
                self.a3 = coeffs_list[2]
                self.a4 = coeffs_list[3]
                self.a6 = coeffs_list[4]
                self.b2 = self.a1**2+4*self.a2
                self.b4 = self.a1*self.a3+2*self.a4
                self.b6 = self.a3**2+4*self.a6
                self.b8 = self.a1**2*self.a6+4*self.a2*self.a6-self.a1*self.a3*self.a4+self.a2*self.a3**2-self.a4**2
                self.c4 = self.b2**2-24*self.b4
                self.c6 = -self.b2**3+36*self.b2*self.b4-216*self.b6
                self.disc = -self.b2**2*self.b8-8*self.b4**3-27*self.b6**2+9*self.b2*self.b4*self.b6
                if self.disc:
                    self.j = self.c4**3*self.disc.inverse()
                else:
                    raise ValueError("coefficients creates singular curve.")
            elif len(self) == 2:
                self.a = coeffs_list[0]
                self.b = coeffs_list[1]
                self.a1 = zero
                self.a2 = zero
                self.a3 = zero
                self.a4 = self.a
                self.a6 = self.b
                self.b2 = zero
                self.b4 = 2*self.a
                self.b6 = 4*self.b
                self.b8 = -(self.a**2)
                self.c4 = -48*self.a
                self.c6 = -864*self.b
                self.disc = -self.b2**2*self.b8-8*self.b4**3-27*self.b6**2+9*self.b2*self.b4*self.b6
                if self.disc:
                    self.j = self.c4**3*self.disc.inverse()
                else:
                    raise ValueError("coefficients creates singular curve.")
            else:
                raise ValueError("coefficient is less or more, can't defined EC.")

        self.ord = None
        self.abelian = None
        self.cubic = _UniVarPolynomial({0:self.a6, 1:self.a4, 2:self.a2, 3:one},
                                      self.basefield)

    def point(self):
        """
        Return a random point on eliiptic curve over ch(field)>3
        """
        bfsize = card(self.basefield)
        one = self.basefield.one
        t = self.basefield.zero
        if len(self) == 2 or (self.a1 == self.a2 == self.a3 == self.basefield.zero):
            while self.basefield.Legendre(t) != 1:
                s = self.basefield.createElement(bigrandom.randrange(bfsize))
                t = self.cubic(s)
                if not t:
                    return [s, t]
            t = self.basefield.sqrt(t)
            r = bigrandom.randrange(2)
            if r:
                return [s, -t]
            return [s, t]
        elif self.ch != 2 and self.ch != 3:
            sform = self.simple()
            while sform.basefield.Legendre(t) != 1:
                s = sform.basefield.createElement(bigrandom.randrange(bfsize))
                t = (s**3+sform.a*s+sform.b)
            x = (s-3*self.b2) // (36*one)
            y = (sform.basefield.sqrt(t) // (108*one)-self.a1*x-self.a3)//(2*one)
            return [x, y]
        elif self.ch == 3:
            while sform.basefield.Legendre(t) != 1:
                s = self.basefield.createElement(bigrandom.randrange(bfsize))
                t = (s**3+self.a2*s**2+self.a4*s+self.a6)
            return [s, self.basefield.sqrt(t)]
        else:
            raise NotImplementedError("This is not implemented.")

    def Schoof(self):
        """
        Return t = q + 1 - #E(F_q). q is basefield order >>1.
        """
        if self.ch in (2, 3):
            raise NotImplementedError("characteristic should be >> 1")

        if len(self) != 2:
            return self.simple().Schoof()

        traces = []
        D, prime_moduli = self.divPoly()
        self.division_polynomials = D # temporary attribute

        # the main loop
        for l in prime_moduli:
            _log.debug("l = %d" % l)
            traces.append(self._Schoof_mod_l(l))

        del self.division_polynomials # temporary attribute deleted

        trace = arith1.CRT(traces)
        modulus = arith1.product(prime_moduli)
        if trace > (modulus >> 1):
            trace -= modulus
        assert abs(trace) <= 2*arith1.floorsqrt(card(self.basefield)), str(self)
        return trace

    def _Schoof_mod2(self):
        """
        Return (#E(F_q) mod 2, 2).  char(F_q) > 3 is required.

        For odd characteristic > 3, t = #E(F_q) mod 2 is determined by
        gcd(self.cubic, X^q - X) == 1 <=> t = 1 mod 2.
        """
        if not self.b:
            result = 0
            _log.debug("(%d, 2) #" % result)
        else:
            linearfactors = _UniVarPolynomial({card(self.basefield):self.basefield.one, 1:-self.basefield.one}, self.basefield)
            if _PolyGCD(self.cubic, linearfactors).degree() == 0:
                result = 1
                _log.debug("(%d, 2) ##" % result)
            else:
                result = 0
                _log.debug("(%d, 2) ###" % result)
        return (result, 2)

    def _Schoof_mod_l(self, l):
        """
        Return q + 1 - #E(F_q) mod l.
        """
        if l == 2:
            return self._Schoof_mod2()
        E = self.cubic
        D = self.division_polynomials
        lth_div = self.division_polynomials[l]
        field = self.basefield
        bfsize = card(field)
        x = _UniVarPolynomial({1:field.one}, field)
        k = bfsize % l
        x_frob = _PolyPow(x, bfsize, lth_div) #x_frob=x^q
        x_frobfrob = _PolyPow(x_frob, bfsize, lth_div) #x_frobfrob=x^{q^2}

        # test for x^{q^2} - x
        f, P = self._sub1(k, x_frobfrob - x, lth_div)
        f0, f3 = f[0], f[3]

        if _PolyGCD(lth_div, P).degree() > 0:
            if arith1.legendre(k, l) == -1:
                _log.debug("%s $" % str((0, l)))
                return (0, l)

            # arith1.legendre(k, l) == 1 <=> k is QR
            w = arith1.modsqrt(k, l)
            f, P = self._sub1(w, x_frob - x, lth_div)

            if _PolyGCD(lth_div, P).degree() == 0: # coprime
                _log.debug("%s $$$$" % str((0, l)))
                return (0, l)

            # there exist non trivial common divisors
            g0 = _PolyPow(E, (bfsize - 1) >> 1, lth_div) #y^(q-1)
            P = self._sub2(w, g0, f[3], lth_div)

            if _PolyGCD(lth_div, P).degree() > 0:
                _log.debug("%s $$" % str((2*w % l, l)))
                return (2*w % l, l)
            else:
                _log.debug("%s $$$" % str((-2*w % l, l)))
                return (-2*w % l, l)

        else: # coprime (GCD(P, lth_div).degree() == 0)
            Y = x - x_frobfrob
            g0 = _PolyPow(E, (bfsize - 1) >> 1, lth_div) #y^(q-1)
            g1 = _PolyPow(g0, bfsize + 1, lth_div) #y^(q^2-1)
            f = -self._sub2(k, g1, f3, lth_div)
            h1 = _PolyMulRed([f, f], lth_div)
            if k & 1 == 0:
                g = (_PolyMulRed([Y, E, f3], lth_div) - f0) * 4
                h0 = _PolyMulRed([g, g], lth_div)
                aux1 = _PolyMulRed([f0, h0], lth_div) + h1
                X_d = _PolyMulRed([E, f3, h0], lth_div)
            else:
                g = (_PolyMulRed([Y, f3], lth_div) - _PolyMulRed([E, f0], lth_div)) * 4
                h0 = _PolyMulRed([g, g], lth_div)
                aux1 = _PolyMulRed([E, _PolyMulRed([f0, h0], lth_div) + h1], lth_div)
                X_d = _PolyMulRed([f3, h0], lth_div)
            X_n = _PolyMulRed([X_d, x_frobfrob + x_frob + x], lth_div) - aux1

            # loop of t
            e_q = _PolyPow(self.cubic, bfsize, lth_div)
            for t in range(1, ((l - 1) >> 1) + 1):
                Z_d_x, Z_n_x = self._Z_x(t, D, e_q, bfsize, lth_div)
                # X_n * Z_d_x == X_d * Z_n_x (mod lth_div)?
                if not _PolyMod(X_n * Z_d_x - X_d * Z_n_x, lth_div):
                    break
            else: # loop of t exhausted
                _log.debug("%s @@@" % str((0, l)))
                return (0, l)

            # found: X_n * Z_d_x == X_d * Z_n_x (mod lth_div)
            y0 = _PolyMulRed([-2*x_frobfrob - x, X_d], lth_div) + aux1
            if k & 1 == 0:
                Y_d = _PolyMulRed([E, D[k], g, X_d], lth_div)
            else:
                Y_d = _PolyMulRed([D[k], g, X_d], lth_div)
            Y_n = -_PolyMulRed([g1, Y_d], lth_div) - _PolyMulRed([f, y0], lth_div)
            Z_d_y, Z_n_y = self._Z_y(t, D, g0, bfsize, lth_div)

            # Y_n * Z_d_y == Y_d * Z_n_y (mod lth_div)?
            if _PolyMod(Y_n * Z_d_y - Y_d * Z_n_y, lth_div):
                _log.debug("%s @@" % str((l-t, l)))
                return (l-t, l)
            else:
                _log.debug("%s @" % str((t, l)))
                return (t, l)

    def _sub1(self, k, t, lth_div):
        """
        Compute in F_q[X]/(D_l):
        P = t * D_{k}^2 * f_E + D_{k-1} * D_{k+1}  (k is even)
        P = t * D_{k}^2 + D{k-1} * D{k+1} * f_E    (k is odd)
        and
        return a dict f and P, where f includes
        f[0] = D{k-1} * D{k+1} and f[3] = D_{k}^2
        for later computation.

        The aim is to compute x^q - x mod D_{l} or x^{q^2} - x mod D_{l}.
        """
        D = self.division_polynomials
        f = {}
        f[0] = _PolyMulRed([D[k-1], D[k+1]], lth_div)
        f[3] = _PolyMulRed([D[k], D[k]], lth_div)
        f[4] = _PolyMulRed([t, f[3]], lth_div)
        if k & 1 == 0:
            P = _PolyMulRed([f[4], self.cubic], lth_div) + f[0]
        else:
            P = f[4] + _PolyMulRed([f[0], self.cubic], lth_div)
        return f, P

    def _sub2(self, k, g, t, lth_div):
        """
        Compute in F_q[X]/(D[l]):
        P = 4 * g * t * D_{k} * f_E^2 - f_1 + f_2  (k is even)
        P = 4 * g * t * D_{k} - f_1 + f_2    (k is odd)
        where
        f_1 = D_{k-1}^2 * D_{k+2}
        f_2 = D_{k-2} * D_{k+1}^2
        and return P.

        The aim is to compute y^q - y mod D_{l} or y^{q^2} - y mod D_{l}.
        """
        D = self.division_polynomials
        fk = _PolyMulRed([4 * g, t, D[k]], lth_div)
        f1 = _PolyMulRed([D[k-1], D[k-1], D[k+2]], lth_div)
        f2 = _PolyMulRed([D[k-2], D[k+1], D[k+1]], lth_div)
        if k & 1 == 0:
            cubic = self.cubic
            P = _PolyMulRed([fk, cubic, cubic], lth_div) - f1 + f2
        else:
            P = fk - f1 + f2
        return P

    def _Z_x(self, t, D, e_q, bfsize, lth_div):
        d = _PolyPow(D[t], 2 * bfsize, lth_div)
        n = _PolyPow(_PolyMulRed([D[t-1], D[t+1]], lth_div), bfsize, lth_div)
        if t & 1 == 0:
            d = _PolyMulRed([e_q, d], lth_div)
        else:
            n = _PolyMulRed([e_q, n], lth_div)
        return d, n

    def _Z_y(self, t, D, g0, bfsize, lth_div):
        E = self.cubic
        if t & 1 == 0:
            d = _PolyPow(_PolyMulRed([4 * E, E, D[t], D[t], D[t]], lth_div), bfsize, lth_div)
        else:
            d = _PolyPow(_PolyMulRed([4 * D[t], D[t], D[t]], lth_div), bfsize, lth_div)
        z0 = _PolyPow(_PolyMulRed([D[t-1], D[t-1], D[t+2]], lth_div) - _PolyMulRed([D[t-2], D[t+1], D[t+1]], lth_div), bfsize, lth_div)
        n = _PolyMulRed([g0, z0], lth_div)
        return d, n

    def _step(self, P, W):
        """
        Return three component list [A, B, C] used in Shanks_Mestre,
        where A is a list of x-coordinates of [q+1]P to [q+1+W]P,
        B is a list of x-coordinates of [0]P, [W]P to [W*W]P and
        C is the intersection set of A and B.
        """
        Qa = self.mul(self.ch + 1, P)
        if len(Qa) == 1:
            _log.debug("[%d]P is zero" % (self.ch + 1))
        Qb = R = self.mul(W, P)
        A = [Qa[0]]
        B = [0, Qb[0]] # 0 = [0][0] ((infinity)_x)
        for i in range(1, W):
            Qa = self.add(Qa, P)
            Qb = self.add(Qb, R)
            A.append(Qa[0])
            B.append(Qb[0])
            if len(Qa) == 1:
                _log.debug("[%d]P is zero" % (self.ch + 1 + i))
                break
        return [A, B, set(A).intersection(set(B))]

    def Shanks_Mestre(self):
        """
        Return t = p + 1 - #E(F_p) for prime field F_p.

        Use in the range 229 <= self.ch <=10**5+o(1).

        This method is using
        Algorithm 7.5.3(Shanks-Mestre assessment of curve order)
        Crandall & Pomerance, PRIME NUMBERS
        """
        if type(self.basefield) != finitefield.FinitePrimeField:
            raise NotImplementedError("FIXME: this is only implemented for FinitePrimeField.")

        if self.ch <= 229:
            return self.naive()
        else: #E.ch>229
            if len(self) != 2:
                return self.simple().Shanks_Mestre()
            g = 0
            while arith1.legendre(g, self.ch) != -1:
                g = bigrandom.randrange(2, self.ch)
            W = arith1.floorsqrt(2 * (arith1.floorsqrt(self.ch) + 1)) + 1

            c, d = g**2*self.a, g**3*self.b
            f = self.cubic
            used = set()
            while True:
                x = bigrandom.randrange(self.ch)
                while x in used:
                    x = bigrandom.randrange(self.ch)
                used.add(x)
                y2 = f(x)
                if not y2:
                    continue
                if arith1.legendre(y2.n, self.ch) == 1:
                    E = self
                    cg = 1
                else: #arith1.legendre(f(x),self.ch)==-1
                    E = EC([c, d], self.ch)
                    cg = -1
                    x = g*x % self.ch
                x = self.basefield.createElement(x)
                P = [x, E.coordinateY(x)]
                A, B, S = E._step(P, W)
                _log.debug("#S = %d" % len(S))
                if len(S) == 1:
                    s = S.pop()
                    if B.count(s) <= 2:
                        break
            aa = A.index(s)
            bb = B.index(s)
            t = aa - bb*W
            if E.mul(self.ch + 1 + t, P) == [0]:
                return -cg*t
            else:
                t = aa + bb*W
                return -cg*t

    def naive(self):
        """
        Return t = p + 1 - #E(Fp).

        This method computes t by counting up Legendre symbols and thus
        needs exponential time.

        RESTRICTIONS:
        The field cannot be of characteristic two nor three.
        """
        if self.ch in (2, 3):
            raise TypeError("cannot be defined over characteristic two/three.")

        sform = self.simple()

        k = 0
        f = sform.cubic
        for i in range(card(sform.basefield)):
            x = sform.basefield.createElement(i)
            k += sform.basefield.Legendre(f(x))
        return -k

    def _order_2(self):
        """
        Return #E(F_2).

        E is in long Weierstrass form:
        Y^2 + a1 XY + a3 Y = X^3 + a2 X^2 + a4 X + a6
        """
        assert self.ch == 2
        result = 1 # for infinite place
        if not self.a6: # (0, 0)
            result += 1
        if self.a3 != self.a6: # (0, 1)
            result += 1
        if self.a2 + self.a4 + self.a6: # (1, 0)
            result += 1
        if self.a1 + self.a3 == self.a2 + self.a4 + self.a6: # (1, 1)
            result += 1
        return result

    def _order_3(self):
        """
        Return #E(F_3).

        E is in long Weierstrass form:
        Y^2 + a1 XY + a3 Y = X^3 + a2 X^2 + a4 X + a6
        """
        assert self.ch == 3
        one = self.basefield.one
        result = 1 # for infinite place
        if not self.a6: # (0, 0)
            result += 1
        if one + self.a3 == self.a6: # (0, 1)
            result += 1
        if one - self.a3 == self.a6: # (0, 2)
            result += 1
        if not (one + self.a2 + self.a4 + self.a6): # (1, 0)
            result += 1
        if self.a1 + self.a3 == self.a2 + self.a4 + self.a6: # (1, 1)
            result += 1
        if -(self.a1 + self.a3) == self.a2 + self.a4 + self.a6: # (1, 2)
            result += 1
        if one == self.a2 - self.a4 + self.a6: # (2, 0)
            result += 1
        if -self.a1 + self.a3 == one + self.a2 - self.a4 + self.a6: # (2, 1)
            result += 1
        if self.a1 - self.a3 == one + self.a2 - self.a4 + self.a6: # (2, 2)
            result += 1
        return result

    def _order_to_trace(self, order):
        """
        Return the trace calculated from the given order:
          t = q + 1 - #E(F_q)
        """
        return card(self.basefield) + 1 - order

    def _trace_to_order(self, trace):
        """
        Return the order calculated from the given trace:
          #(E_q) = q + 1 - t
        """
        return card(self.basefield) + 1 - trace

    def trace(self, index=None):
        """
        Return the Frobenius trace t = q + 1 - #E(F_q),
        where q is basefield cardinality.

        If index is an integer greater than 1, then return the trace
        t = q^r + 1 - #E(F_q^r) for a subfield curve defined over F_q.
        """
        bfsize = card(self.basefield)

        if not self.ord:
            if bfsize == self.ch: # prime field
                # special cases
                if bfsize == 2 or bfsize == 3:
                    trace = self._order_to_trace(self.order())
                # trace main block
                elif bfsize < 10**4:
                    trace = self.naive()
                elif bfsize < 10**30:
                    trace = self.Shanks_Mestre()
                else: # self.ch>=10**30
                    trace = self.Schoof()
            else:
                if self.ch in (2, 3):
                    error_message = "no E/F_{%d} trace" % bfsize
                    raise NotImplementedError(error_message)
                else:
                    trace = self.Schoof()

            self.ord = self._trace_to_order(trace) # cached
        else:
            trace = self._order_to_trace(self.ord)

        # final result
        if index is not None:
            # for subfield curve
            basetrace = trace
            trace, oldtrace = basetrace, 2
            for i in range(2, index + 1):
                trace, oldtrace = basetrace*trace - bfsize*oldtrace, trace

        return trace

    def order(self, index=None):
        """
        Return #E(F_q) or #E(F_{q^r}).

        E is defined over F_q.
        If the method is called as E.order(), the result is #E(F_q).
        If the method is called as E.order(r), the result is #E(F_{q^r}).

        RESTRICTIONS:
        F_q cannot be q = 2^k or q = 3^k for k > 1.
        """
        bfsize = card(self.basefield)

        if not self.ord:
            if self.ch in (2, 3):
                if bfsize == self.ch == 2:
                    self.ord = self._order_2()
                elif bfsize == self.ch == 3:
                    self.ord = self._order_3()
                else:
                    error_message = "no E/F_{%d} order" % bfsize
                    raise NotImplementedError(error_message)
            else:
                self.ord = self._trace_to_order(self.trace())

        # final result
        if index:
            # for subfield curve
            basetrace = self._order_to_trace(self.ord)
            trace, oldtrace = basetrace, 2
            for i in range(2, index + 1):
                trace, oldtrace = basetrace*trace - bfsize*oldtrace, trace
            return bfsize ** index + 1 - trace

        return self.ord

    def line(self, P, Q, R):
        """
        this use to compute weil pairing.
        return line_{P,Q}(R).
        """
        if not R or R == self.infpoint:
            raise ValueError("R must not be zero")
        if P == Q == self.infpoint:
            return self.basefield.one
        if P == self.infpoint:
            return R[0] - Q[0]
        if Q == self.infpoint:
            return R[0] - P[0]
        if P[0] != Q[0]:
            return (Q[0] - P[0]) * R[1] - (Q[1] - P[1]) * \
                     R[0]- Q[0] * P[1] + P[0] * Q[1]
        if P[1] != Q[1]:
            return R[0] - P[0]
        return (3 * P[0] ** 2 + 2 * self.a2 * P[0] + self.a4 - self.a1 * P[1] ) * \
                R[0] - (2 * P[1] + self.a1 * P[0] + self.a3 ) * R[1] - \
                (P[0] ** 3) + self.a4 * P[0] + 2 * self.a6 - self.a3 * P[1]

    def pointorder(self, P, ord_factor=None):
        """
        find point order of P and return order.
        """
        # parameter ord_factor is extension for structre.
        if ord_factor is not None:
            N = int(ord_factor)
        else:
            N = self.order()
            ord_factor = factor_misc.FactoredInteger(N)
        o = 1
        for p, e in ord_factor:
            B = self.mul(N//(p**e), P)
            while B != self.infpoint:
                o = o*p
                B = self.mul(p, B)
        return o

    def Miller(self, P, m, Q, R):
        """
        return f_P(D_Q)
        this is used for Tate/Weil pairing

        suppose that P be an m-torsion point .
        support point R must be in neither groups generated by P nor Q.

        return False if parameters lack any conditions above.

        NOTE: algorithm limitations forced that assume R is not P-Q .
        """
        # check order
        if m < 2 or not (m % self.ch):
            raise ValueError("order more than 1 and not be divisible by characteristic")

        O = self.infpoint

        # support point must not be P-Q
        S = self.add(R, Q)
        if S == O:
            return False

        # j = 1
        jP = P
        v = self.basefield.one
        for k in arith1.expand(m, 2)[-2::-1]:
            j2P = self.mul(2, jP)
            denominator = self.line(jP, jP, R) * self.line(j2P, O, S)
            if not denominator:
                return False
            numerator = self.line(jP, jP, S) * self.line(j2P, O, R)
            if not numerator:
                return False
            f = numerator / denominator
            v = v**2 * f
            # j *= 2
            jP = j2P
            if k:
                kjP = self.add(P, jP)
                denominator = self.line(P, jP, R) * self.line(kjP, O, S)
                if not denominator:
                    return False
                numerator = self.line(P, jP, S) * self.line(kjP, O, R)
                if not numerator:
                    return False
                f = numerator / denominator
                v = v * f
                # j += 1
                jP = kjP
        # now j == m
        return v

    def TatePairing(self, m, P, Q):
        """
        computing the Tate-Lichetenbaum pairing with Miller's algorithm.
        parameters satisfies that mul(m,P)==[0].
        """
        O = self.infpoint
        if self.mul(m, P) != O:
            raise ValueError("sorry, not mP=[0].")

        if P == O or Q == O:
            return self.basefield.one

        forbidden = [O, P, self.mul(-1, Q), self.sub(P, Q)]
        R = self.add(P, Q)
        T = False
        while (not T):
            while R in forbidden:
                R = self.point()
            forbidden.append(R)
            T = self.Miller(P, m, Q, R)
        return T

    def TatePairing_Extend(self, m, P, Q):
        """
        computing the Tate-Lichtenbaum pairing extended original Tate Pairing.
        """
        return self.TatePairing(m, P, Q)**((card(self.basefield)-1)//m)

    def WeilPairing(self, m, P, Q):
        """
        computing the Weil pairing with Miller's algorithm.
        we assume point P and Q that be in different m-tortion group .
        """
        O = self.infpoint
        if self.mul(m, P) != O or self.mul(m, Q) != O:
            raise ValueError("sorry, not mP=[0] or mQ=[0].")

        if P == O or Q == O or P == Q:
            return self.basefield.one

        T = U = False
        forbidden = [O, P, Q, self.sub(Q, P)]
        R = self.sub(P,Q) # assume Q not in group generated P
        while (not T) or (not U):
            while R in forbidden:
                R = self.point()
            T = self.Miller(Q, m, P, R)
            U = self.Miller(P, m, Q, self.mul(-1, R))
            forbidden.append(R)
        F = U/T
        return F

    def BSGS(self, P):
        """
        returns order of P such that kP=[0]
        refered to Washington 4.3.4.
        """
        if P == self.infpoint:
            return 1

        bfsize = card(self.basefield)

        Q = self.mul(bfsize + 1, P)
        m = arith1.floorpowerroot(bfsize, 4) + 1
        Plist = [self.infpoint]
        R = P
        j = 1
        while j <= m:
            Plist.append(R)
            R = self.add(R, P)
            j += 1
        R = self.mul(2*m, P)
        k = -m
        Plist_rev = map(self.mul, [-1]*(m+1), Plist) # make reverse point mapping
        j = 0
        while k <= m:
            S = self.add(Q, self.mul(k, R))
            if S in Plist:
                j = Plist.index(S)
                break
            elif S in Plist_rev:
                j = -Plist_rev.index(S)
                break
            k += 1
        M = self.ch + 1 + 2*m*k - j
        Flist = factor_methods.factor(M)
        for p, e in Flist:
            for i in range(e):
                if self.mul(M//p, P) == self.infpoint:
                    M = M//p
        return M

    def DLP_BSGS(self, n, P, Q):
        """
        returns k such that Q=kP
        P,Q is the elements of the same group
        """
        B = []
        m = arith1.floorsqrt(n) + 1
        R = Q
        B.append(Q)
        i = 1
        while i <= m:
            R = self.sub(R, P)
            if R == [0]:
                return i
            B.append(R)
            i += 1
        R = self.mul(m, P)
        P = R
        if R in B:
            return m+B.index(R)
        else:
            i = 2
            while i <= m:
                R = self.add(R, P)
                if R in B:
                    return i*m+B.index(R)
                i = i+1
            if self.mul(n, P) == Q:
                return n
            else:
                return False

    def structure(self):
        """
        returns group structure E(K)=Z_d x Z_m with d|m|#E(K)
        """
        if self.abelian:
            return self.abelian

        # step 1. find order E/F_p.
        simplified = self.simple()
        N = simplified.order()
        if prime.primeq(N):
            return (1, N)

        # step 2. decompose N.
        r = gcd.gcd(card(simplified.basefield) - 1, N)
        _log.debug("r = %d, N = %d" % (r, N))
        r_factor = factor_misc.FactoredInteger(r)
        N1, N2 = 1, N
        for p in r_factor.prime_divisors():
            k, N2 = arith1.vp(N2, p=p)
            N1 *= p**k

        _log.debug("loop")
        while 1:
            P1 = self.infpoint
            while P1 == self.infpoint:
                P1 = simplified.point()
            P2 = self.infpoint
            while P2 == self.infpoint:
                P2 = simplified.point()
            P1, P2 = simplified.mul(N2, P1), simplified.mul(N2, P2)
            s = simplified.pointorder(P1, r_factor)
            t = simplified.pointorder(P2, r_factor)
            m = gcd.lcm(s, t)
            if m > 1:
                e = simplified.WeilPairing(m, P1, P2)
                if e != self.basefield.one:
                    d = e.order()
                else:
                    d = 1
                if m*d == N1:
                    _log.debug("N1 = %d" % N1)
                    _log.debug("P1 = %s (pointorder=%d)" % (P1, s))
                    _log.debug("P2 = %s (pointorder=%d)" % (P2, t))
                    assert (not (N//d) % d), d
                    self.abelian = (d, N//d)
                    return self.abelian

    def issupersingular(self):
        """
        returns supersingularities.
        """
        if self.order() % self.ch == 1:
            return True
        else:
            return False

def EC(coefficient, basefield=None):
    """
    generate new elliptic curve class.
    """
    try:
        character = basefield.getCharacteristic()
        field = basefield
    except:
        # backward compatiblity
        if isinstance(basefield, (int, long)):
            field = finitefield.FinitePrimeField(basefield)
            character = basefield
        elif isinstance(basefield, rational.RationalField) or not basefield:
            character = 0 # necessary to be defined
        else:
            raise ValueError("basefield must be RationalFieid or FiniteField.")

    if isinstance(coefficient, list):
        if not character:
            return ECoverQ(coefficient)
        else:
            return ECoverGF(coefficient, field)
