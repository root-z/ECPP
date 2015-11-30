"""
Utilities for univar.

The module provides higher level interfaces to univar classes and
functions.
"""

from __future__ import division
import logging
import nzmath.arith1 as arith1
import nzmath.bigrandom as bigrandom
import nzmath.gcd as gcd
import nzmath.rational as rational
import nzmath.ring as ring
import nzmath.poly.univar as univar
import nzmath.poly.termorder as termorder
import nzmath.poly.ring as poly_ring
import nzmath.poly.ratfunc as ratfunc
import nzmath.compatibility


_log = logging.getLogger("nzmath.poly.uniutil")


_MIXIN_MSG = "%s is mix-in"


class OrderProvider(object):
    """
    OrderProvider provides order and related operations.
    """
    def __init__(self, order=termorder.ascending_order):
        """
        Do not instantiate OrderProvider.
        This initializer should be called from descendant:
          OrderProvider.__init__(self, order)
        where order is default to termorder.ascending_order.
        """
        if type(self) is OrderProvider:
            raise NotImplementedError(_MIXIN_MSG % self.__class__.__name__)
        self.order = order

    def shift_degree_to(self, degree):
        """
        Return polynomial whose degree is the given degree.

        More precisely, let f(X) = a_0 + ... + a_n * X**n, then
        order.shift_degree_to(f, m) returns:
        - zero polynomial, if f is zero polynomial
        - a_(n-m) + ... + a_n * X**m, if 0 <= m < n
        - a_0 * X**(m-n) + ... + a_n * X**m, if m >= n
        """
        original_degree = self.order.degree(self)
        difference = degree - original_degree
        if difference == 0:
            return self
        elif difference < 0:
            return self.construct_with_default([(d + difference, c) for (d, c) in self if d + difference >= 0])
        else:
            return self.upshift_degree(difference)

    def split_at(self, degree):
        """
        Return tuple of two polynomials, which are splitted at the
        given degree.  The term of the given degree, if exists,
        belongs to the lower degree polynomial.
        """
        lower, upper = [], []
        for d, c in self:
            if self.order.cmp(d, degree) <= 0:
                lower.append((d, c))
            else:
                upper.append((d, c))
        return (self.construct_with_default(lower),
                self.construct_with_default(upper))


class DivisionProvider(object):
    """
    DivisionProvider provides all kind of divisions for univariate
    polynomials.  It is assumed that the coefficient ring of the
    polynomials is a field.

    The class should be used as a mix-in.

    REQUIRE: OrderProvider
    """
    def __init__(self):
        """
        Do not instantiate DivisionProvider.
        This initializer should be called from descendant.
        """
        if type(self) is DivisionProvider:
            raise NotImplementedError(_MIXIN_MSG % self.__class__.__name__)
        self._reduced = {}

    def __divmod__(self, other):
        """
        divmod(self, other)

        The method does, as the built-in divmod, return a tuple of
        (self // other, self % other).
        """
        degree, lc = self.order.leading_term(other)
        quotient, remainder = [], self
        while self.order.degree(remainder) >= degree:
            rdegree, rlc = self.order.leading_term(remainder)
            q = rdegree - degree, rlc / lc
            remainder = remainder - other.term_mul(q)
            quotient.append(q)
        quotient = self.construct_with_default(quotient)
        return quotient, remainder

    def __floordiv__(self, other):
        """
        self // other
        """
        degree, lc = self.order.leading_term(other)
        quotient, remainder = [], self
        while self.order.degree(remainder) >= degree:
            rdegree, rlc = self.order.leading_term(remainder)
            q = rdegree - degree, rlc / lc
            remainder = remainder - other.term_mul(q)
            quotient.append(q)
        return self.construct_with_default(quotient)

    def __mod__(self, other):
        """
        self % other
        """
        degree, lc = self.order.leading_term(other)
        remainder = self
        rdegree, rlc = self.order.leading_term(remainder)
        if rdegree > degree + 5 and degree >= 1:
            return other.mod(remainder)
        ilc = ring.inverse(lc)
        while rdegree >= degree:
            q = rdegree - degree, rlc * ilc
            remainder = remainder - other.term_mul(q)
            rdegree, rlc = self.order.leading_term(remainder)
        return remainder

    def mod(self, dividend):
        """
        Return dividend % self.

        self should have attribute _reduced to cache reduced monomials.
        """
        degree, lc = self.order.leading_term(self)
        div_deg = self.order.degree(dividend)
        if div_deg < degree:
            return dividend
        upperbound = min(degree * 2, div_deg) + 1
        populate = False
        if not self._reduced:
            populate = True
        else:
            i = min(self._reduced.keys())
            while i in self._reduced.keys():
                i += 1
            if i < upperbound:
                populate = True
        if populate:
            self._populate_reduced(degree, lc, upperbound)
        if div_deg > degree * 2:
            dividend_degrees = sorted(dividend.iterbases(), reverse=True)
            dividend_degrees = [deg for deg in dividend_degrees if deg > degree]
            self._populate_reduced_more(dividend_degrees)
        accum = self.construct_with_default(())
        lowers = []
        for d, c in dividend:
            if c:
                if d < degree:
                    lowers.append((d, c))
                else:
                    accum += self._reduced[d].scalar_mul(c)
        return accum + self.construct_with_default(lowers)

    def mod_pow(self, polynom, index):
        """
        Return polynom ** index % self.
        """
        if not self:
            raise ZeroDivisionError("polynomial division or modulo by zero.")
        polynom %= self
        if index == 1 or not polynom:
            return polynom
        elif polynom.degree() == 0:
            return self.getRing().createElement([(0, polynom[0]**index)])
        acoefficient = polynom.itercoefficients().next()
        one = ring.getRing(acoefficient).one
        power_product = self.construct_with_default([(0, one)])
        if index:
            power_of_2 = polynom
            while index:
                if index & 1:
                    power_product = self.mod(power_product * power_of_2)
                power_of_2 = self.mod(power_of_2.square())
                index //= 2
        return power_product

    def _populate_reduced(self, degree, lc, upperbound):
        """
        Populate self._reduced.

        degree, lc is of self, and self._reduced is populated up to
        the given upperbound.
        """
        one = ring.getRing(self.itercoefficients().next()).one
        if not self._reduced:
            minimum = degree
            redux = self.construct_with_default([(degree - 1, one)])
        else:
            minimum = max(self._reduced.keys()) + 1
            redux = self._reduced[minimum - 1]
        moniced = self.scalar_mul(one / lc)
        for i in range(minimum, upperbound):
            redux = redux.term_mul((1, 1))
            if self.order.degree(redux) == degree:
                redux -= moniced.scalar_mul(self.order.leading_coefficient(redux))
            self._reduced[i] = redux

    def _populate_reduced_more(self, degrees):
        """
        Populate self._reduced more for much higher degree dividend.
        This method has to be called after _populate_reduced.

        self._reduced is populated so that it will include all of
        degrees > degree * 2.  The degrees are recommended to be
        sorted in descending order.
        """
        degree = self.order.degree(self)
        if not [deg for deg in degrees if deg not in self._reduced]:
            # no need to populate more
            return
        # self._reduced has every key from degree up to maxreduced.
        # The default (prepared by _populate_reduced) is degree*2.
        maxreduced = degree * 2
        while (maxreduced + 1) in self._reduced:
            maxreduced += 1
        # use binary chain multiplied by a stride
        stride = maxreduced
        if len(degrees) > 1:
            # try to use wider stride
            common_degree = gcd.gcd_of_list(degrees)[0]
            if common_degree > maxreduced:
                # common_degree is better than maxreduced.
                if common_degree not in self._reduced:
                    self._populate_reduced_more([common_degree])
                stride = common_degree
        redux = self._reduced[stride]
        maximum = max(degrees)
        binary = {stride:redux}
        while stride * 2 <= maximum:
            stride += stride
            if stride not in self._reduced:
                redux = self.mod(redux.square())
            else:
                redux = self._reduced[stride]
            binary[stride] = redux
        binarykeys = sorted(binary.keys(), reverse=True)
        for deg in (d for d in degrees if d > maxreduced):
            pickup = []
            rest = deg
            for key in binarykeys:
                if rest < key:
                    continue
                rest -= key
                pickup.append(key)
                if rest < maxreduced or rest in self._reduced:
                    break
            assert pickup
            total = pickup.pop() # start with the smallest degree picked up
            prod = binary[total]
            for picked in reversed(pickup):
                total += picked
                prod = self.mod(prod * binary[picked])
                self._reduced[total] = prod
            if rest in self._reduced:
                final = self.mod(prod * self._reduced[rest])
            else: # rest < degree
                final = self.mod(prod.term_mul((rest, 1)))
            self._reduced[deg] = final

    def __truediv__(self, other):
        """
        self / other

        The result is a rational function.
        """
        return ratfunc.RationalFunction(self, other)

    def scalar_exact_division(self, scale):
        """
        Return quotient by a scalar which can divide each coefficient
        exactly.
        """
        return self.coefficients_map(lambda c: c.exact_division(scale))

    def gcd(self, other):
        """
        Return a greatest common divisor of self and other.
        Returned polynomial is always monic.
        """
        if self.order.degree(self) < self.order.degree(other):
            divident, divisor = other, self
        else:
            divident, divisor = self, other
        while divisor:
            divident, divisor = divisor, divident % divisor
        lc = self.order.leading_coefficient(divident)
        if lc and lc != ring.getRing(lc).one:
            divident = divident.scalar_exact_division(lc)
        return divident

    def extgcd(self, other):
        """
        Return a tuple (u, v, d); they are the greatest common divisor
        d of two polynomials self and other and u, v such that
        d = self * u + other * v.

        see nzmath.gcd.extgcd
        """
        order = termorder.ascending_order
        polynomial_ring = self.getRing()
        zero, one = polynomial_ring.zero, polynomial_ring.one
        a, b, g, u, v, w = one, zero, self, zero, one, other
        while w:
            q = g // w
            a, b, g, u, v, w = u, v, w, a - q*u, b - q*v, g - q*w
        lc = self.order.leading_coefficient(g)
        if lc and lc != 1:
            linv = lc.inverse()
            a, b, g = linv * a, linv * b, linv * g
        return (a, b, g)


class PseudoDivisionProvider(object):
    """
    PseudoDivisionProvider provides pseudo divisions for univariate
    polynomials.  It is assumed that the coefficient ring of the
    polynomials is a domain.

    The class should be used as a mix-in.
    """
    def pseudo_divmod(self, other):
        """
        self.pseudo_divmod(other) -> (Q, R)

        Q, R are polynomials such that
        d**(deg(self) - deg(other) + 1) * self == other * Q + R,
        where d is the leading coefficient of other.
        """
        order = termorder.ascending_order
        if hasattr(self, 'order'):
            assert self.order is order
        degree, lc = order.leading_term(other)
        # step 1
        quotient, remainder = self.construct_with_default(()), self
        rdegree, rlc = order.leading_term(remainder)
        e = rdegree - degree + 1
        if e <= 0:
            return quotient, remainder
        while rdegree >= degree:
            # step 3
            canceller = self.construct_with_default([(rdegree - degree, rlc)])
            quotient = quotient.scalar_mul(lc) + canceller
            remainder = remainder.scalar_mul(lc) - canceller * other
            e -= 1
            rdegree, rlc = order.leading_term(remainder)
        # step 2
        q = lc ** e
        quotient = quotient.scalar_mul(q)
        remainder = remainder.scalar_mul(q)
        return quotient, remainder

    def pseudo_floordiv(self, other):
        """
        self.pseudo_floordiv(other) -> Q

        Q is a polynomial such that
        d**(deg(self) - deg(other) + 1) * self == other * Q + R,
        where d is the leading coefficient of other and R is a
        polynomial.
        """
        order = termorder.ascending_order
        if hasattr(self, 'order'):
            assert self.order is order
        degree, lc = order.leading_term(other)
        # step 1
        quotient, remainder = self.construct_with_default(()), self
        rdegree, rlc = order.leading_term(remainder)
        e = order.degree(remainder) - degree + 1
        if e <= 0:
            return quotient
        while rdegree >= degree:
            # step 3
            canceller = self.construct_with_default([(rdegree - degree, rlc)])
            quotient = quotient.scalar_mul(lc) + canceller
            remainder = remainder.scalar_mul(lc) - canceller * other
            e -= 1
            rdegree, rlc = order.leading_term(remainder)
        # step 2
        return quotient.scalar_mul(lc ** e)

    def pseudo_mod(self, other):
        """
        self.pseudo_mod(other) -> R

        R is a polynomial such that
        d**(deg(self) - deg(other) + 1) * self == other * Q + R,
        where d is the leading coefficient of other and Q a
        polynomial.
        """
        order = termorder.ascending_order
        if hasattr(self, 'order'):
            assert self.order is order
        degree, lc = order.leading_term(other)
        # step 1
        remainder = self
        rdegree, rlc = order.leading_term(remainder)
        e = order.degree(remainder) - degree + 1
        if e <= 0:
            return remainder
        while rdegree >= degree:
            # step 3
            canceller = self.construct_with_default([(rdegree - degree, rlc)])
            remainder = remainder.scalar_mul(lc) - canceller * other
            e -= 1
            rdegree, rlc = order.leading_term(remainder)
        # step 2
        return remainder.scalar_mul(lc ** e)

    def __truediv__(self, other):
        """
        self / other

        Return the result as a rational function.
        """
        order = termorder.ascending_order
        return ratfunc.RationalFunction(self, other)

    def exact_division(self, other):
        """
        Return quotient of exact division.
        """
        order = termorder.ascending_order
        if hasattr(self, 'order'):
            assert self.order is order
        quotient, remainder = self.pseudo_divmod(other)
        if not remainder:
            deg_other, lc = order.leading_term(other)
            deg_self = order.degree(self)
            extra_factor = lc ** (deg_self - deg_other + 1)
            return quotient.scalar_exact_division(extra_factor)
        raise ValueError("division is not exact")

    def scalar_exact_division(self, scale):
        """
        Return quotient by a scalar which can divide each coefficient
        exactly.
        """
        return self.coefficients_map(lambda c: c.exact_division(scale))

    def monic_divmod(self, other):
        """
        self.monic_divmod(other) -> (Q, R)

        Q, R are polynomials such that
        self == other * Q + R.

        The leading coefficient of other MUST be one.
        """
        order = termorder.ascending_order
        if hasattr(self, 'order'):
            assert self.order is order

        degree = order.degree(other)
        quotient, remainder = self.construct_with_default(()), self
        rdegree, rlc = order.leading_term(remainder)
        if rdegree < degree:
            return quotient, remainder
        while rdegree >= degree:
            # step 3
            canceller = self.construct_with_default([(rdegree - degree, rlc)])
            quotient += canceller
            remainder -= canceller * other
            rdegree, rlc = order.leading_term(remainder)

        return quotient, remainder

    def monic_floordiv(self, other):
        """
        self.monic_floordiv(other) -> Q

        Q is a polynomial such that
        self == other * Q + R,
        where R is a polynomial whose degree is smaller than other's.

        The leading coefficient of other MUST be one.
        """
        order = termorder.ascending_order
        if hasattr(self, 'order'):
            assert self.order is order
        degree = order.degree(other)
        # step 1
        quotient, remainder = self.construct_with_default(()), self
        rdegree, rlc = order.leading_term(remainder)
        if rdegree < degree:
            return quotient
        while rdegree >= degree:
            # step 3
            canceller = self.construct_with_default([(rdegree - degree, rlc)])
            quotient += canceller
            remainder -= canceller * other
            rdegree, rlc = order.leading_term(remainder)

        return quotient

    def monic_mod(self, other):
        """
        self.monic_mod(other) -> R

        R is a polynomial such that
        self == other * Q + R,
        where Q is another polynomial.

        The leading coefficient of other MUST be one.
        """
        order = termorder.ascending_order
        if hasattr(self, 'order'):
            assert self.order is order
        degree = order.degree(other)
        remainder = self
        rdegree, rlc = order.leading_term(remainder)
        if rdegree < degree:
            return remainder
        while rdegree >= degree:
            canceller = self.construct_with_default([(rdegree - degree, rlc)])
            remainder -= canceller * other
            rdegree, rlc = order.leading_term(remainder)

        return remainder

    def monic_pow(self, index, mod):
        """
        Return self**index % mod.

        The leading coefficient of mod MUST be one.
        """
        if index < 0:
            raise ValueError("negative index is not allowed.")
        # special indices
        elif index == 0:
            return self.getRing().one
        elif index == 1:
            return self.monic_mod(mod)
        elif index == 2:
            return self.square().monic_mod(mod)
        # special polynomials
        if not self:
            return self
        # general
        power_product = self.getRing().one
        power_of_2 = self.monic_mod(mod)
        while index:
            if index & 1:
                power_product = (power_product * power_of_2).monic_mod(mod)
            index //= 2
            if index:
                power_of_2 = power_of_2.square().monic_mod(mod)
        return power_product


class ContentProvider(object):
    """
    ContentProvider provides content and primitive part.

    The coefficient ring must be a unique factorization domain.
    """
    def content(self):
        """
        Return content of the polynomial.
        """
        coeffring = None
        isquotfield = None
        cont = 0
        num, den = 0, 1
        for c in self.itercoefficients():
            if isquotfield is None:
                coeffring = ring.getRing(c)
                if coeffring.isfield() and isinstance(coeffring, ring.QuotientField):
                    isquotfield = True
                    basedomain = coeffring.basedomain
                else:
                    isquotfield = False
            if isquotfield:
                num = basedomain.gcd(num, c.numerator)
                den = basedomain.lcm(den, c.denominator)
            else:
                if not cont:
                    cont = c
                else:
                    cont = coeffring.gcd(cont, c)
        if isquotfield:
            cont = coeffring.createElement(num, den)
        return cont

    def primitive_part(self):
        """
        Return the primitive part of the polynomial.
        """
        return self.scalar_exact_division(self.content())


class SubresultantGcdProvider(object):
    """
    SubresultantGcdProvider provides gcd method using subresultant
    algorithm.

    REQUIRE: PseudoDivisionProvider, ContentProvider
    """
    def resultant(self, other):
        """
        Return the resultant of self and other.
        """
        order = termorder.ascending_order
        f, g = self, other
        negative = False
        if order.degree(f) < order.degree(g):
            f, g = g, f
            if (order.degree(f) & 1) and (order.degree(g) & 1):
                negative = not negative
        coeff_ring = self.getCoefficientRing()
        a = b = coeff_ring.one
        while order.degree(g) > 0:
            delta = order.degree(f) - order.degree(g)
            if (order.degree(f) & 1) and (order.degree(g) & 1):
                negative = not negative
            h = f.pseudo_mod(g)
            f, g = g, h.scalar_exact_division(a * (b ** delta))
            a = order.leading_coefficient(f)
            b = ((a ** delta) * b).exact_division(b ** delta)
        if not g:
            return coeff_ring.zero
        scalar = g[0]
        degree_f = order.degree(f)
        if negative:
            return (-b * scalar ** degree_f).exact_division(b ** degree_f)
        else:
            return (b * scalar ** degree_f).exact_division(b ** degree_f)

    def subresultant_gcd(self, other):
        """
        Return the greatest common divisor of given polynomials.  They
        must be in the polynomial ring and its coefficient ring must
        be a UFD.

        Reference: Algorithm 3.3.1 of Cohen's book
        """
        order = termorder.ascending_order
        divident = self
        divisor = other
        polynomial_ring = self.getRing()
        coeff_ring = polynomial_ring.getCoefficientRing()
        one = coeff_ring.one
        # step 1
        if order.degree(divisor) > order.degree(divident):
            divident, divisor = divisor, divident
        if not divisor:
            return divident
        content_gcd = coeff_ring.gcd(divident.content(), divisor.content())
        divident = divident.primitive_part()
        divisor = divisor.primitive_part()
        g = h = one

        while 1:
            # step 2
            delta = order.degree(divident) - order.degree(divisor)
            quotient, remainder = divident.pseudo_divmod(divisor)
            if not remainder:
                return divisor.primitive_part().scalar_mul(content_gcd)
            if order.degree(remainder) == 0:
                return polynomial_ring.createElement(content_gcd)

            # step 3
            divident, divisor = divisor, remainder
            if delta == 0 and g != one:
                divisor = divisor.scalar_exact_division(g)
            elif delta == 1 and (g != one or h != one):
                divisor = divisor.scalar_exact_division(g * h)
            elif delta > 1 and (g != one or h != one):
                divisor = divisor.scalar_exact_division(g * h**delta)
            g = divident.leading_coefficient()
            if delta == 1 and h != g:
                h = g
            elif delta > 1 and (g != one or h != one):
                h = g * (h * g) ** (delta - 1)

    def subresultant_extgcd(self, other):
        """
        Return (A, B, P) s.t. A*self+B*other=P,
        where P is the greatest common divisor of given polynomials.
        They must be in the polynomial ring and its coefficient ring must
        be a UFD.

        Reference: Kida's paper p.18
        """
        order = termorder.ascending_order
        coeff_ring = self.getCoefficientRing()
        P_1, P_2 = self, other
        A_1, A_2 = coeff_ring.one, coeff_ring.zero
        if order.degree(P_1) < order.degree(P_2):
            P_1, P_2 = P_2, P_1
            A_1, A_2 = A_2, A_1
        if not P_2:
            return (A_1, A_2, P_1)
        a = b = coeff_ring.one
        while P_2:
            delta = order.degree(P_1) - order.degree(P_2)
            b = a * (b ** delta)
            Q, R = P_1.pseudo_divmod(P_2)
            a = order.leading_coefficient(P_2)
            ad = a ** delta
            P_1, P_2 = P_2, R.scalar_exact_division(b)
            A_1, A_2 = A_2, (ad * a * A_1 - Q * A_2).scalar_exact_division(b)
            b = (ad * b) // (b ** delta)
        B_1 = (P_1 - A_1 * self).exact_division(other)
        return (A_1, B_1, P_1)


class PrimeCharacteristicFunctionsProvider(object):
    """
    PrimeCharacteristicFunctionsProvider provides efficient powering
    and factorization for polynomials whose coefficient ring has prime
    characteristic.

    - A client of this mix-in class should use DivisionProvider also.
    - A client of this mix-in class must have attribute ch, which
      stores the prime characteristic of the coefficient ring.
    """
    def __init__(self, ch):
        """
        Do not instantiate PrimeCharacteristicFunctionsProvider.
        This initializer should be called from descendant.
        """
        if type(self) is PrimeCharacteristicFunctionsProvider:
            raise NotImplementedError(_MIXIN_MSG % self.__class__.__name__)
        self.ch = ch

    def __pow__(self, index, mod=None):
        """
        self ** index

        A client of the mix-in class should write the name in front of
        the main class so that the __pow__ defined here would be
        preceded in method resolution order.  For example:

        class Client (PrimeCharacteristicFunctionsProvider
                      DivisionProvider,
                      Polynomial):
        """
        if not isinstance(index, (int, long)):
            raise TypeError("index must be an integer.")
        if index < 0:
            raise ValueError("index must be a non-negative integer.")
        if mod is not None:
            return mod.mod_pow(self, index)
        if index > 0:
            if not self:
                return self
            q = 1
            while index % self.ch == 0:
                q *= self.ch
                index //= self.ch
            if q > 1:
                powered = self.construct_with_default([(d * q, c ** q) for (d, c) in self])
            else:
                powered = self
        if index == 1:
            return powered
        acoefficient = self.itercoefficients().next()
        one = ring.getRing(acoefficient).one
        power_product = self.construct_with_default([(0, one)])
        if index:
            power_of_2 = powered
            while index:
                if index & 1:
                    power_product *= power_of_2
                power_of_2 = power_of_2.square()
                index //= 2
        return power_product

    def mod_pow(self, polynom, index):
        """
        Return polynom ** index % self.
        """
        if not self:
            raise ZeroDivisionError("polynomial division or modulo by zero.")
        polynom %= self
        if index == 1 or not polynom:
            return polynom
        elif polynom.degree() == 0:
            return self.getRing().createElement([(0, polynom[0]**index)])
        acoefficient = polynom.itercoefficients().next()
        cardfq = card(ring.getRing(acoefficient))
        final_product = self.getRing().one
        qpow = polynom
        # q-th power expansion
        while index:
            index, small = divmod(index, cardfq)
            if small == 1:
                final_product *= qpow
            elif small:
                final_product *= self._small_index_mod_pow(qpow, small)
            # c ** q = c for any Fq element c, thus
            qpow = self.mod(qpow.bases_map(lambda d: d * cardfq))

        return final_product

    def _small_index_mod_pow(self, polynom, index):
        """
        Return polynom ** index % self for small index (0 < index < q).
        """
        # binary powering
        power_product = self.getRing().one
        if index > 0:
            power_of_2 = polynom
            while index:
                if index & 1:
                    power_product = self.mod(power_product * power_of_2)
                    if index == 1:
                        break
                power_of_2 = self.mod(power_of_2.square())
                index //= 2
        else:
            raise ValueError("this private method requires index > 0")
        return power_product

    def squarefree_decomposition(self):
        """
        Return the square free decomposition of the polynomial.  The
        return value is a dict whose keys are integers and values are
        corresponding powered factors.  For example, if
        A = A1 * A2**2,
        the result is {1: A1, 2: A2}.

        gcd method is required.
        """
        result = {}
        if self.order.degree(self) == 1:
            return {1: self}
        f = self
        df = f.differentiate()
        if df:
            b = f.gcd(df)
            a = f.exact_division(b)
            i = 1
            while self.order.degree(a) > 0:
                c = a.gcd(b)
                b = b.exact_division(c)
                if a != c:
                    r = a.exact_division(c)
                    if self.order.degree(r) > 0:
                        result[i] = r
                    a = c
                i += 1
            f = b
        if self.order.degree(f) > 0:
            f = f.pthroot()
            subresult = f.squarefree_decomposition()
            for i, g in subresult.iteritems():
                result[i*self.ch] = g
        return result

    def pthroot(self):
        """
        Return a polynomial obtained by sending X**p to X, where p is
        the characteristic.  If the polynomial does not consist of
        p-th powered terms only, result is nonsense.
        """
        return self.construct_with_default([(d // self.ch, c) for (d, c) in self if c])

    def distinct_degree_factorization(self):
        """
        Return the distinct degree factorization of the polynomial.
        The return value is a dict whose keys are integers and values
        are corresponding product of factors of the degree.  For
        example, if A = A1 * A2, and all irreducible factors of A1
        having degree 1 and all irreducible factors of A2 having
        degree 2, then the result is:
          {1: A1, 2: A2}.

        The given polynomial must be square free, and its coefficient
        ring must be a finite field.
        """
        Fq = ring.getRing(self.itercoefficients().next())
        q = card(Fq)
        f = self
        x = f.construct_with_default([(1, Fq.one)])
        w = x
        i = 1
        result = {}
        while i*2 <= self.order.degree(f):
            w = pow(w, q, f)
            result[i] = f.gcd(w - x)
            if self.order.degree(result[i]) > 0:
                f = f.exact_division(result[i])
                w = w % f
            else:
                del result[i]
            i += 1
        if self.order.degree(f) != 0:
            result[self.order.degree(f)] = f
        return result

    def split_same_degrees(self, degree):
        """
        Return the irreducible factors of the polynomial.  The
        polynomial must be a product of irreducible factors of the
        given degree.
        """
        r = self.order.degree(self) // degree
        Fq = ring.getRing(self.itercoefficients().next())
        q = card(Fq)
        p = Fq.getCharacteristic()
        if degree == 1:
            result = []
            X = self.construct_with_default([(1, Fq.one)])
            f = self
            while not f[0]:
                f = f // X
                result.append(X)
            if self.order.degree(f) >= 1:
                result.append(f)
        else:
            result = [self]
        while len(result) < r:
            # g is a random polynomial
            rand_coeff = {}
            for i in range(2 * degree):
                rand_coeff[i] = Fq.createElement(bigrandom.randrange(q))
            if not rand_coeff[2 * degree - 1]:
                rand_coeff[2 * degree - 1] = Fq.one
            randpoly = self.construct_with_default(rand_coeff)
            if p == 2:
                G = self.construct_with_default(())
                for i in range(degree):
                    G = G + randpoly
                    randpoly = self.mod(randpoly.square())
            else:
                one = self.construct_with_default([(0, Fq.one)])
                G = pow(randpoly, (q**degree - 1) >> 1, self) - one
            subresult = []
            while result:
                h = result.pop()
                if self.order.degree(h) == degree:
                    subresult.append(h)
                    continue
                z = h.gcd(G)
                if 0 < self.order.degree(z) < self.order.degree(h):
                    subresult.append(z)
                    subresult.append(h.exact_division(z))
                else:
                    subresult.append(h)
            result = subresult
        return result

    def factor(self):
        """
        Factor the polynomial.

        The returned value is a list of tuples whose first component
        is a factor and second component is its multiplicity.
        """
        result = []
        lc = self.order.leading_coefficient(self)
        if lc != ring.getRing(lc).one:
            self = self.scalar_exact_division(lc)
            result.append((self.getRing().one*lc, 1))
        squarefreefactors = self.squarefree_decomposition()
        for m, f in squarefreefactors.iteritems():
            distinct_degree_factors = f.distinct_degree_factorization()
            for d, g in distinct_degree_factors.iteritems():
                if d == self.order.degree(g):
                    result.append((g, m))
                else:
                    for irred in g.split_same_degrees(d):
                        result.append((irred, m))
        return result

    def isirreducible(self):
        """
        f.isirreducible() -> bool

        If f is irreducible return True, otherwise False.
        """
        if 0 <= self.order.degree(self) <= 1:
            # degree 1 polynomials and constants are irreducible
            return True
        elif not self[0]:
            # degree >= 2 polynomials are reducible if constant term is zero
            return False
        squareFreeFactors = self.squarefree_decomposition()
        if len(squareFreeFactors) != 1:
            return False
        m, f = squareFreeFactors.popitem()
        if m != 1:
            return False
        distinctDegreeFactors = f.distinct_degree_factorization()
        if len(distinctDegreeFactors) != 1:
            return False
        d, g = distinctDegreeFactors.popitem()
        if d != self.order.degree(g):
            return False
        return True


class KaratsubaProvider(object):
    """
    define Karatsuba method multiplication and squaring.

    REQUIRE: OrderProvider
    """
    def ring_mul_karatsuba(self, other):
        """
        Multiplication of two polynomials in the same ring.

        Computation is carried out by Karatsuba method.
        """
        # zero
        if not self or not other:
            return self.construct_with_default(())
        # monomial
        if len(self) == 1:
            return other.term_mul(self)
        if len(other) == 1:
            return self.term_mul(other)
        # binomial
        if len(self) == 2:
            p, q = [other.term_mul(term) for term in self]
            return p + q
        if len(other) == 2:
            p, q = [self.term_mul(term) for term in other]
            return p + q
        # suppose self is black and other is red.
        black_least_degree = self.order.tail_degree(self)
        black_most_degree = self.order.degree(self)
        red_least_degree = self.order.tail_degree(other)
        red_most_degree = self.order.degree(other)
        least_degree = min(black_least_degree, red_least_degree)
        most_degree = max(black_most_degree, red_most_degree)
        assert least_degree < most_degree
        half_degree = (least_degree + most_degree) >> 1

        if black_least_degree > half_degree:
            return self.downshift_degree(black_least_degree).ring_mul_karatsuba(other).upshift_degree(black_least_degree)
        if red_least_degree > half_degree:
            return self.ring_mul_karatsuba(other.downshift_degree(red_least_degree)).upshift_degree(red_least_degree)

        black = self.downshift_degree(least_degree)
        red = other.downshift_degree(least_degree)
        club, spade = black.split_at(half_degree - least_degree)
        dia, heart = red.split_at(half_degree - least_degree)
        weaker = club.ring_mul_karatsuba(dia)
        stronger = spade.ring_mul_karatsuba(heart)
        karatsuba = (club + spade).ring_mul_karatsuba(dia + heart) - weaker - stronger
        return (weaker.upshift_degree(least_degree * 2) +
                karatsuba.upshift_degree(least_degree + half_degree) +
                stronger.upshift_degree(half_degree * 2))

    def square_karatsuba(self):
        """
        Return the square of self.
        """
        # zero
        if not self:
            return self

        polynomial = self.construct_with_default
        data_length = len(self)
        # monomial
        if data_length == 1:
            d, c = iter(self).next()
            if d:
                return polynomial([(d*2, c**2)])
            else:
                return polynomial([(0, c**2)])
        # binomial
        if data_length == 2:
            (d1, c1), (d2, c2) = [(d, c) for (d, c) in self]
            if "_sorted" in self._init_kwds:
                del self._init_kwds["_sorted"]
            return polynomial({d1*2:c1**2, d1+d2:c1*c2*2, d2*2:c2**2}, _sorted=False)
        # general (Karatsuba)
        least_degree = self.order.tail_degree(self)
        if least_degree:
            chopped = self.downshift_degree(least_degree)
        else:
            chopped = self
        half_degree = self.order.degree(self) >> 1
        fst, snd = chopped.split_at(half_degree)
        fst_squared = fst.square()
        snd_squared = snd.square()
        karatsuba = (fst + snd).square() - fst_squared - snd_squared
        result = (fst_squared +
                  karatsuba.upshift_degree(half_degree) +
                  snd_squared.upshift_degree(half_degree * 2))
        if least_degree:
            return result.upshift_degree(least_degree)
        else:
            return result


class VariableProvider(object):
    """
    VariableProvider provides the variable name and other cariable
    related methods.
    """
    def __init__(self, varname):
        """
        Do not instantiate VariableProvider.
        This initializer should be called from descendant.
        """
        if type(self) is VariableProvider:
            raise NotImplementedError(_MIXIN_MSG % self.__class__.__name__)
        self.variable = varname

    def getVariable(self):
        """
        Get variable name
        """
        return self.variable

    def getVariableList(self):
        """
        Get variable name as list.
        """
        return [self.variable]


class RingElementProvider(ring.CommutativeRingElement):
    """
    Provides interfaces for ring.CommutativeRingElement.
    """
    def __init__(self):
        """
        Do not instantiate RingElementProvider.
        This initializer should be called from descendant:
          RingElementProvider.__init__(self)
        """
        if type(self) is RingElementProvider:
            raise NotImplementedError(_MIXIN_MSG % self.__class__.__name__)
        ring.CommutativeRingElement.__init__(self)
        self._coefficient_ring = None
        self._ring = None

    def getRing(self):
        """
        Return an object of a subclass of Ring, to which the element
        belongs.
        """
        if self._coefficient_ring is None or self._ring is None:
            myring = None
            for c in self.itercoefficients():
                cring = ring.getRing(c)
                if not myring or myring != cring and myring.issubring(cring):
                    myring = cring
                elif not cring.issubring(myring):
                    myring = myring.getCommonSuperring(cring)
            if not myring:
                myring = rational.theIntegerRing
            self.set_coefficient_ring(myring)
        return self._ring

    def set_coefficient_ring(self, coeffring):
        if self._coefficient_ring is None:
            self._coefficient_ring = coeffring
            # variable names are ignored now
            self._ring = poly_ring.PolynomialRing.getInstance(self._coefficient_ring)


class RingPolynomial(OrderProvider,
                     univar.SortedPolynomial,
                     RingElementProvider):
    """
    General polynomial with commutative ring coefficients.
    """
    def __init__(self, coefficients, coeffring=None, _sorted=False, **kwds):
        """
        Initialize the polynomial.

        - coefficients: initializer for polynomial coefficients
        - coeffring: commutative ring
        """
        if coeffring is None:
            raise TypeError("argument `coeffring' is required")
        coefficients = dict(coefficients)
        if coefficients and coefficients.itervalues().next() not in coeffring:
            coefficients = [(d, coeffring.createElement(c)) for (d, c) in coefficients.iteritems()]
            _sorted = False
        kwds["coeffring"] = coeffring
        univar.SortedPolynomial.__init__(self, coefficients, _sorted, **kwds)
        OrderProvider.__init__(self, termorder.ascending_order)
        RingElementProvider.__init__(self)
        self.set_coefficient_ring(coeffring)

    def getRing(self):
        """
        Return an object of a subclass of Ring, to which the element
        belongs.
        """
        # short-cut self._ring is None case
        return self._ring

    def getCoefficientRing(self):
        """
        Return an object of a subclass of Ring, to which the all
        coefficients belong.
        """
        # short-cut self._coefficient_ring is None case
        return self._coefficient_ring

    def __repr__(self):
        return "%s(%s, %s)" % (self.__class__.__name__, repr(self.terms()), repr(self._coefficient_ring))

    def __add__(self, other):
        try:
            return univar.SortedPolynomial.__add__(self, other)
        except AttributeError:
            one = self.getRing().one
            try:
                return univar.SortedPolynomial.__add__(self, other * one)
            except Exception:
                return NotImplemented

    def __radd__(self, other):
        one = self.getRing().one
        try:
            return other * one + self
        except Exception:
            return NotImplemented

    def __sub__(self, other):
        try:
            return univar.SortedPolynomial.__sub__(self, other)
        except AttributeError:
            one = self.getRing().one
            try:
                return univar.SortedPolynomial.__sub__(self, other * one)
            except Exception:
                return NotImplemented

    def __rsub__(self, other):
        one = self.getRing().one
        try:
            return other * one - self
        except Exception:
            return NotImplemented

    def ismonic(self):
        """
        Return True if the polynomial is monic, i.e. its leading
        coefficient is one.  False otherwise.
        """
        return self.leading_coefficient() == self._coefficient_ring.one

    def __getitem__(self, degree):
        """
        Return the coefficient of specified degree.
        If there is no term of degree, return 0.
        """
        result = univar.SortedPolynomial.__getitem__(self, degree)
        if result is 0:
            result = self._coefficient_ring.zero
        return result


class DomainPolynomial(PseudoDivisionProvider,
                       RingPolynomial):
    """
    Polynomial with domain coefficients.
    """
    def __init__(self, coefficients, coeffring=None, _sorted=False, **kwds):
        """
        Initialize the polynomial.

        - coefficients: initializer for polynomial coefficients
        - coeffring: domain
        """
        if coeffring is None:
            raise TypeError("argument `coeffring' is required")
        elif not coeffring.isdomain():
            raise TypeError("coefficient ring has to be a domain")
        RingPolynomial.__init__(self, coefficients, coeffring, _sorted, **kwds)
        PseudoDivisionProvider.__init__(self)

    def discriminant(self):
        """
        Return discriminant of the polynomial.
        """
        df = self.differentiate()
        lc = self.leading_coefficient()
        if self.degree() & 3 in (2, 3):
            sign = -1
        else:
            sign = 1
        if self.getCoefficientRing().getCharacteristic() == 0:
            if lc != 1:
                return sign * self.resultant(df) / lc
            else:
                return sign * self.resultant(df)
        else:
            return sign * self.resultant(df) * lc**(self.degree() - df.degree() - 2)

    def to_field_polynomial(self):
        """
        Return a FieldPolynomial object obtained by embedding the
        polynomial ring over the domain D to over the quatient field
        of D.
        """
        field = self.getCoefficientRing().getQuotientField()
        return polynomial([(d, field.createElement(c)) for d, c in self.iterterms()], field)

    def __pow__(self, index, mod=None):
        """
        self ** index (% mod)

        It overrides the method from SortedPolynomial.  The mod MUST
        be monic, otherwise the method raises ValueError.
        """
        if mod is None:
            return RingPolynomial.__pow__(self, index)
        elif mod.ismonic():
            return self.monic_pow(index, mod)
        else:
            raise ValueError("non-monic modulus")


class UniqueFactorizationDomainPolynomial(SubresultantGcdProvider,
                                          ContentProvider,
                                          DomainPolynomial):
    """
    Polynomial with unique factorization domain coefficients.
    """
    def __init__(self, coefficients, coeffring=None, _sorted=False, **kwds):
        """
        Initialize the polynomial.

        - coefficients: initializer for polynomial coefficients
        - coeffring: unique factorization domain
        """
        if coeffring is None:
            raise TypeError("argument `coeffring' is required")
        elif not coeffring.isufd():
            raise TypeError("coefficient ring has to be a UFD")
        DomainPolynomial.__init__(self, coefficients, coeffring, _sorted, **kwds)
        ContentProvider.__init__(self)
        SubresultantGcdProvider.__init__(self)


class IntegerPolynomial(UniqueFactorizationDomainPolynomial):
    """
    Polynomial with integer coefficients.

    This class is required because special initialization must be done
    for built-in int/long.
    """
    def __init__(self, coefficients, coeffring=None, _sorted=False, **kwds):
        dc = dict(coefficients)
        coefficients = [(d, rational.IntegerIfIntOrLong(c)) for (d, c) in dc.iteritems()]
        UniqueFactorizationDomainPolynomial.__init__(self, coefficients, coeffring, _sorted, **kwds)

    def normalize(self):
        """
        returns the unique normalized polynomial g
        which is associate to self (so g=u*self for some unit in coeffring).
        For IntegerPolynomial, g is positive.
        """
        if self.leading_coefficient() < 0:
            return -self
        return self

class FieldPolynomial(DivisionProvider,
                      ContentProvider,
                      RingPolynomial):
    """
    Polynomial with field coefficients.
    """
    def __init__(self, coefficients, coeffring=None, _sorted=False, **kwds):
        """
        Initialize the polynomial.

        - coefficients: initializer for polynomial coefficients
        - coeffring: field
        """
        if coeffring is None:
            raise TypeError("argument `coeffring' is required")
        elif not coeffring.isfield():
            raise TypeError("coefficient ring has to be a field")
        RingPolynomial.__init__(self, coefficients, coeffring, _sorted, **kwds)
        DivisionProvider.__init__(self)
        ContentProvider.__init__(self)

    def __pow__(self, index, mod=None):
        """
        self ** index (% mod).
        """
        if mod is None:
            return RingPolynomial.__pow__(self, index)

        if index < 0:
            raise ValueError("negative index is not allowed.")
        # special indices
        elif index == 0:
            return self.getRing().one
        elif index == 1:
            return self % mod
        elif index == 2:
            return self.square() % mod
        # special polynomials
        if not self:
            return self
        # general
        power_product = self.getRing().one
        power_of_2 = self % mod
        while index:
            if index & 1:
                power_product = mod.mod(power_product * power_of_2)
            index //= 2
            if index:
                power_of_2 = mod.mod(power_of_2.square())
        return power_product

    def resultant(self, other):
        """
        Return the resultant of self and other.
        """
        order = termorder.ascending_order
        f, g = self, other
        negative = False
        if order.degree(f) < order.degree(g):
            f, g = g, f
            if (order.degree(f) & 1) and (order.degree(g) & 1):
                negative = not negative
        coeff_ring = self.getCoefficientRing()
        a = b = coeff_ring.one
        while order.degree(g) > 0:
            delta = order.degree(f) - order.degree(g)
            if (order.degree(f) & 1) and (order.degree(g) & 1):
                negative = not negative
            h = f % g
            h *= order.leading_coefficient(g)**(order.degree(f) - order.degree(g) + 1)
            f, g = g, h.scalar_exact_division(a * (b ** delta))
            a = order.leading_coefficient(f)
            b = ((a ** delta) * b) // (b ** delta)
        if not g:
            return coeff_ring.zero
        scalar = g[0]
        degree_f = order.degree(f)
        if negative:
            return -(b * scalar ** degree_f) // (b ** degree_f)
        else:
            return (b * scalar ** degree_f) // (b ** degree_f)

    def discriminant(self):
        """
        Return discriminant of the polynomial.
        """
        df = self.differentiate()
        lc = self.leading_coefficient()
        if self.degree() & 3 in (2, 3):
            sign = -1
        else:
            sign = 1
        if self.getCoefficientRing().getCharacteristic() == 0:
            return sign * self.resultant(df) / lc
        else:
            return sign * self.resultant(df) * lc**(self.degree() - df.degree() - 2)


class FiniteFieldPolynomial(PrimeCharacteristicFunctionsProvider,
                            FieldPolynomial):
    """
    Fq polynomial
    """
    def __init__(self, coefficients, coeffring=None, _sorted=False, **kwds):
        """
        Initialize the polynomial.

        - coefficients: initializer for polynomial coefficients
        - coeffring: finite prime field
        """
        if coeffring is None:
            raise TypeError("argument `coeffring' is required")
        FieldPolynomial.__init__(self, coefficients, coeffring, _sorted, **kwds)
        PrimeCharacteristicFunctionsProvider.__init__(self, coeffring.getCharacteristic())


class FinitePrimeFieldPolynomial(FiniteFieldPolynomial):
    """
    Fp polynomial
    """
    def __init__(self, coefficients, coeffring=None, _sorted=False, **kwds):
        """
        Initialize the polynomial.

        - coefficients: initializer for polynomial coefficients
        - coeffring: finite prime field
        """
        FiniteFieldPolynomial.__init__(self, coefficients, coeffring, _sorted, **kwds)


    # IMPLEMENTATION REMARK:
    # The most time consuming part of computation is bunch of
    # object creations.  Thus, here, the creation of finite field
    # elements is avoided during summing up coefficients.

    def ring_mul(self, other):
        """
        Multiplication of two polynomials in the same ring.
        """
        try:
            mul_coeff = {}
            ## stripped to bare integer
            stripped_self = [(ds, cs.n) for (ds, cs) in self]
            stripped_other = [(do, co.n) for (do, co) in other]
            for ds, cs in stripped_self:
                for do, co in stripped_other:
                    term_degree = ds + do
                    if term_degree in mul_coeff:
                        mul_coeff[term_degree] = (mul_coeff[term_degree] + cs*co ) % self.ch
                    else:
                        mul_coeff[term_degree] = cs*co
            ## back to decorated datatype automatically
            return self.construct_with_default(mul_coeff)
        except AttributeError:
            # maybe fail due to lacking attribute .n sometimes
            _log.debug("fall back to good old ring_mul")
            return univar.PolynomialInterface.ring_mul(self, other)

    def square(self):
        """
        Return the square of self.
        """
        # zero
        if not self:
            return self
        # monomial, binomial
        elif len(self.sorted) <= 2:
            return FieldPolynomial.square(self)

        # general
        try:
            result = {}
            # stripped to bare integer
            stripped_self = [(ds, cs.n) for (ds, cs) in self]
            for i, term in enumerate(stripped_self):
                di, ci = term[0], term[1]
                result[di*2] = result.get(di*2, 0) + pow(ci, 2, self.ch)
                for j in range(i + 1, len(stripped_self)):
                    dj, cj = self.sorted[j][0], self.sorted[j][1]
                    result[di + dj] = result.get(di + dj, 0) + 2 * ci * cj
            # back to decorated datatype automatically
            return self.construct_with_default(result)
        except AttributeError:
            _log.debug("fall back to old square")
            return FieldPolynomial.square(self)


def inject_variable(polynom, variable):
    """
    Inject variable into polynom temporarily.  The variable name will
    be lost after any arithmetic operations on polynom, though the
    class name of polynom will remain prefixed with 'Var'.  If one need
    variable name permanently, he/she should define a class inheriting
    VariableProvider.
    """
    baseclasses = polynom.__class__.__bases__
    if VariableProvider not in baseclasses:
        polynom.__class__ = type("Var" + polynom.__class__.__name__,
                                 (polynom.__class__, VariableProvider,),
                                 {})
    polynom.variable = variable


special_ring_table = {rational.IntegerRing: IntegerPolynomial}


def polynomial(coefficients, coeffring):
    """
    Return a polynomial.
    - coefficients has to be a initializer for dict, whose keys are
      degrees and values are coefficients at degrees.
    - coeffring has to be an object inheriting ring.Ring.

    One can override the way to choose a polynomial type from a
    coefficient ring, by setting:
    special_ring_table[coeffring_type] = polynomial_type
    before the function call.
    """
    if type(coeffring) in special_ring_table:
        poly_type = special_ring_table[type(coeffring)]
    elif coeffring.isfield():
        poly_type = FieldPolynomial
    elif coeffring.isufd():
        poly_type = UniqueFactorizationDomainPolynomial
    elif coeffring.isdomain():
        poly_type = DomainPolynomial
    else:
        poly_type = RingPolynomial
    return poly_type(coefficients, coeffring)


def OneVariableDensePolynomial(coefficient, variable, coeffring=None):
    """
    OneVariableDensePolynomial(coefficient, variable [,coeffring])

    - coefficient has to be a sequence of coefficients in ascending order
    of degree.
    - variable has to be a character string.
    - coeffring has to be, if specified, an object inheriting ring.Ring.

    This function is provided for backward compatible way of defining
    univariate polynomial.  The argument variable is ignored.
    """
    _coefficients = dict(enumerate(coefficient))
    if coeffring is None:
        coeffring = init_coefficient_ring(_coefficients)
    return polynomial(_coefficients, coeffring)

def OneVariableSparsePolynomial(coefficient, variable, coeffring=None):
    """
    OneVariableSparsePolynomial(coefficient, variable [,coeffring])

    - coefficient has to be a dictionary of degree-coefficient pairs.
    - variable has to be a character string.
    - coeffring has to be, if specified, an object inheriting ring.Ring.

    This function is provided for backward compatible way of defining
    univariate polynomial.  The argument variable is ignored.
    """
    _coefficients = dict(coefficient)
    if coeffring is None:
        coeffring = init_coefficient_ring(_coefficients)
    return polynomial(_coefficients, coeffring)

def init_coefficient_ring(coefficients):
    """
    Return a ring to which all coefficients belong.  The argument
    coefficients is a dictionary whose values are the coefficients.
    """
    myRing = None
    for c in coefficients.itervalues():
        cring = ring.getRing(c)
        if not myRing or myRing != cring and myRing.issubring(cring):
            myRing = cring
        elif not cring.issubring(myRing):
            myRing = myRing.getCommonSuperring(cring)
    return myRing
