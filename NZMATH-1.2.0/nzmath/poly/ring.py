"""
base classes for polynomial rings and rational function fields.
"""

from __future__ import division
import nzmath.ring as ring
import nzmath.poly.termorder as termorder
import nzmath.poly.univar as univar
import nzmath.poly.multivar as multivar


class PolynomialRing(ring.CommutativeRing):
    """
    The class of uni-/multivariate polynomial ring.
    There's no need to specify the variable names.
    """

    _instances = {}

    def __init__(self, coeffring, number_of_variables=1):
        """
        PolynomialRing(coeffring)
        creates a polynomial ring for univariate polynomials, while
        PolynomialRing(coeffring, n)
        creates a polynomial ring for multivariate polynomials.
        """
        if not isinstance(coeffring, ring.Ring):
            raise TypeError("%s should not be passed as ring" % coeffring.__class__)
        ring.CommutativeRing.__init__(self)
        self._coefficient_ring = coeffring
        if number_of_variables == 1 and self._coefficient_ring.isfield():
            self.properties.setIseuclidean(True)
        if self._coefficient_ring.isufd():
            self.properties.setIsufd(True)
        if self._coefficient_ring.isnoetherian():
            self.properties.setIsnoetherian(True)
        elif self._coefficient_ring.isdomain() in (True, False):
            self.properties.setIsdomain(self._coefficient_ring.isdomain())
        self.number_of_variables = number_of_variables

    def getCoefficientRing(self):
        """
        Return the coefficient ring.
        """
        return self._coefficient_ring

    def getQuotientField(self):
        """
        Return the quotient field of the ring if coefficient ring has
        its quotient field.  Otherwise, an exception will be raised.
        """
        coefficientField = self._coefficient_ring.getQuotientField()
        return RationalFunctionField(coefficientField, self.number_of_variables)

    def __eq__(self, other):
        """
        equality test
        """
        if self is other:
            return True
        if (isinstance(other, PolynomialRing) and
            self._coefficient_ring == other._coefficient_ring and
            self.number_of_variables == other.number_of_variables):
            return True
        return False

    def __ne__(self, other):
        """
        not equal
        """
        return not self.__eq__(other)

    def __hash__(self):
        """
        hash(self)
        """
        return (hash(self._coefficient_ring) ^ (self.number_of_variables * hash(self.__class__.__name__) + 1)) & 0x7fffffff

    def __repr__(self):
        """
        Return 'PolynomialRing(Ring, #vars)'
        """
        if self.number_of_variables == 1:
            return "%s(%s)" % (self.__class__.__name__, repr(self._coefficient_ring))
        return "%s(%s, %d)" % (self.__class__.__name__, repr(self._coefficient_ring), self.number_of_variables)

    def __str__(self):
        """
        Return R[][]
        """
        return str(self._coefficient_ring) + "[]" * self.number_of_variables

    def __contains__(self, element):
        """
        `in' operator is provided for checking the element be in the
        ring.
        """
        if element in self._coefficient_ring:
            return True
        elem_ring = ring.getRing(element)
        if elem_ring is not None and elem_ring.issubring(self):
            return True
        return False

    def issubring(self, other):
        """
        reports whether another ring contains this polynomial ring.
        """
        if self is other:
            return True
        elif isinstance(other, PolynomialRing):
            if (self._coefficient_ring.issubring(other.getCoefficientRing()) and
                self.number_of_variables <= other.number_of_variables):
                return True
        elif other.issubring(self._coefficient_ring):
            return False
        try:
            return other.issuperring(self)
        except RuntimeError:
            # reach max recursion by calling each other
            return False

    def issuperring(self, other):
        """
        reports whether this polynomial ring contains another ring.
        """
        if self is other:
            return True
        elif self._coefficient_ring.issuperring(other):
            return True
        elif isinstance(other, PolynomialRing):
            return (self._coefficient_ring.issuperring(other.getCoefficientRing()) and
                    self.number_of_variables >= other.number_of_variables)
        try:
            return other.issubring(self)
        except RuntimeError:
            # reach max recursion by calling each other
            return False

    def getCommonSuperring(self, other):
        """
        Return common superring of two rings.
        """
        if self.issuperring(other):
            return self
        elif other.issuperring(self):
            return other
        elif (not isinstance(other, PolynomialRing) and
              other.issuperring(self._coefficient_ring)):
            return self.__class__(other, self.number_of_variables)
        try:
            return other.getCommonSuperring(self)
        except RuntimeError:
            # reached recursion limit by calling on each other
            raise TypeError("no common super ring")
        except AttributeError:
            # other doesn't have getCommonSuperring
            raise TypeError("%s is not a ring" % str(other))

    def getCharacteristic(self):
        """
        Return characteristic of the ring.
        """
        return self._coefficient_ring.getCharacteristic()

    def createElement(self, seed):
        """
        Return an element in the ring made from seed.
        """
        if ring.getRing(seed) == self:
            return seed
        elif not seed:
            return self._zero_polynomial()
        elif seed in self._coefficient_ring:
            return self._constant_polynomial(seed)
        else:
            return self._prepared_polynomial(seed)

    def _zero_polynomial(self):
        """
        Return the zero polynomial in the polynomial ring.
        """
        if self.number_of_variables == 1:
            import nzmath.poly.uniutil as uniutil
            return uniutil.polynomial((), self._coefficient_ring)
        else:
            import nzmath.poly.multiutil as multiutil
            return multiutil.polynomial((), coeffring=self._coefficient_ring, number_of_variables=self.number_of_variables)

    def _constant_polynomial(self, seed):
        """
        Return a constant polynomial made from a constant seed.
        seed should not be zero.
        """
        if self.number_of_variables == 1:
            import nzmath.poly.uniutil as uniutil
            return uniutil.polynomial({0: seed}, self._coefficient_ring)
        else:
            import nzmath.poly.multiutil as multiutil
            const = (0,) * self.number_of_variables
            return multiutil.polynomial({const: seed}, self._coefficient_ring)

    def _prepared_polynomial(self, preparation):
        """
        Return a polynomial from given preparation, which is suited
        for the first argument of uni-/multi-variable polynomials.
        """
        if self.number_of_variables == 1:
            import nzmath.poly.uniutil as uniutil
            return uniutil.polynomial(preparation, self._coefficient_ring)
        else:
            import nzmath.poly.multiutil as multiutil
            return multiutil.polynomial(preparation, self._coefficient_ring)

    def _get_one(self):
        "getter for one"
        if self._one is None:
            self._one = self._constant_polynomial(self._coefficient_ring.one)
        return self._one

    one = property(_get_one, None, None, "multiplicative unit")

    def _get_zero(self):
        "getter for zero"
        if self._zero is None:
            self._zero = self._zero_polynomial()
        return self._zero

    zero = property(_get_zero, None, None, "additive unit")

    def gcd(self, a, b):
        """
        Return the greatest common divisor of given polynomials.
        The polynomials must be in the polynomial ring.
        If the coefficient ring is a field, the result is monic.
        """
        if hasattr(a, "gcd"):
            return a.gcd(b)
        elif hasattr(a, "subresultant_gcd"):
            return a.subresultant_gcd(b)
        raise NotImplementedError("no gcd")

    def extgcd(self, a, b):
        """
        Return the tuple (u, v, d): d is the greatest common divisor
        of given polynomials, and they satisfy d = u*a + v*b. The
        polynomials must be in the polynomial ring.  If the
        coefficient ring is a field, the result is monic.
        """
        if hasattr(a, "extgcd"):
            return a.extgcd(b)
        raise NotImplementedError("no extgcd")

    @classmethod
    def getInstance(cls, coeffring, number_of_variables=1):
        """
        Return an instance of the class with specified coefficient ring
        and number of variables.
        """
        if (coeffring, number_of_variables) not in cls._instances:
            cls._instances[coeffring, number_of_variables] = cls(coeffring, number_of_variables)
        return cls._instances[coeffring, number_of_variables]


class PolynomialIdeal(ring.Ideal):
    """
    A class to represent an ideal of univariate polynomial ring.
    """
    def __init__(self, generators, polyring):
        """
        Initialize an ideal in the given polyring (polynomial ring).

        generators: a generator polynomial or a list of polynomials
        """
        if type(polyring) is not PolynomialRing:
            raise TypeError("polyring has to be an instance of PolynomialRing")
        ring.Ideal.__init__(self, generators, polyring)
        self._normalize_generators()

    def __contains__(self, elem):
        """
        Return whether elem is in the ideal or not.
        """
        if not elem.getRing().issubring(self.ring):
            return False
        if self.generators == [self.ring.zero]:
            return elem == self.ring.zero
        return not self.reduce(elem)

    def __nonzero__(self):
        """
        Report whether the ideal is zero ideal or not.  Of course,
        False is for zero ideal.
        """
        return self.generators and self.generators != [self.ring.zero]

    def __repr__(self):
        """
        Return repr string.
        """
        return "%s(%r, %r)" % (self.__class__.__name__, self.generators, self.ring)

    def __str__(self):
        """
        Return str string.
        """
        return "(%s)%s" % (", ".join([termorder.ascending_order.format(g) for g in self.generators]), self.ring)

    def issubset(self, other):
        """
        Report whether another ideal contains this ideal.
        """
        if self is other:
            return True
        for g in self.generators:
            if g not in other:
                return False
        return True

    def issuperset(self, other):
        """
        Report whether this ideal contains another ideal.
        """
        if self is other:
            return True
        for g in other.generators:
            if g not in self:
                return False
        return True

    def reduce(self, element):
        """
        Reduce the given element by the ideal.  The result is an
        element of the class which represents the equivalent class.
        """
        order = termorder.ascending_order
        if isinstance(element, univar.PolynomialInterface):
            reduced = element
        else:
            reduced = self.ring.createElement(element)
        if not self:
            if not reduced:
                reduced = self.ring.zero
        elif len(self.generators) == 1:
            g = self.generators[0]
            if self.ring.iseuclidean():
                if isinstance(g, univar.PolynomialInterface):
                    reduced %= g
                else: # g is element of a field
                    reduced = self.ring.zero
            elif self.ring.getCoefficientRing().iseuclidean():
                # higher degree to lower degree
                # subtract euclid quotient of lc * generator
                if isinstance(g, univar.PolynomialInterface):
                    g_degree, lc = order.leading_term(g)
                    degree = order.degree(reduced)
                    for d in range(degree, g_degree - 1, -1):
                        q = reduced[d] // lc
                        reduced -= g.term_mul((d - g_degree, q))
                else:
                    reduced = reduced.coefficients_map(lambda c: c % g)
        elif self.ring.getCoefficientRing().iseuclidean():
            # assert that the generators are sorted descending order
            # of degrees.
            for g in self.generators:
                reduced = self._euclidean_reduce(reduced, g)
        else:
            raise NotImplementedError("should be implemented")
        return reduced

    def _euclidean_reduce(self, element, g):
        """
        Reduce an element of the ring by a polynomial g.
        The coefficient ring has to be a Euclidean domain.
        """
        order = termorder.ascending_order
        coeffring = self.ring.getCoefficientRing()
        reduced = element
        g_degree, g_lc = order.leading_term(g)
        degree, lc = order.leading_term(reduced)
        while degree >= g_degree:
            if not (lc % g_lc):
                reduced -= g.term_mul((degree - g_degree, lc // g_lc))
            else:
                u, v, common_divisor = coeffring.extgcd(lc, g_lc)
                reduced = reduced.scalar_mul(u) + g.term_mul((degree - g_degree, v))
                break
            degree, lc = order.leading_term(reduced)
        return reduced

    def _normalize_generators(self):
        """
        Normalize generators
        """
        order = termorder.ascending_order
        if len(self.generators) > 1 and self.ring.ispid():
            g = self.generators[0]
            for t in self.generators:
                g = self.ring.gcd(g, t)
            self.generators = [g]
        elif self.ring.getCoefficientRing().iseuclidean():
            coeffring = self.ring.getCoefficientRing()
            degree = order.degree
            lc = order.leading_coefficient
            cand_stack = list(self.generators)
            tentative = []
            while cand_stack:
                next = cand_stack.pop()
                if not tentative:
                    tentative.append(next)
                    continue
                next_degree = degree(next)
                while tentative:
                    last = tentative.pop()
                    last_degree = degree(last)
                    if last_degree > next_degree:
                        cand_stack.append(last)
                        continue
                    next_lc, last_lc = lc(next), lc(last)
                    if last_degree == next_degree:
                        u, v, d = coeffring.extgcd(next_lc, last_lc)
                        # make new polynomial whose lc = d
                        head = next.scalar_mul(u) + last.scalar_mul(v)
                        next -= head.scalar_mul(next_lc // d)
                        last -= head.scalar_mul(last_lc // d)
                        assert degree(next) < next_degree
                        assert degree(last) < last_degree
                        cand_stack.append(head)
                        if next:
                            cand_stack.append(next)
                        if last:
                            cand_stack.append(last)
                        break
                    elif not (next_lc % last_lc):
                        next -= last.scalar_mul(next_lc // last_lc)
                        cand_stack.append(next)
                        break
                    elif last_lc % next_lc:
                        u, v, d = coeffring.extgcd(next_lc, last_lc)
                        next = next.scalar_mul(u) + last.term_mul((next_degree - last_degree, v))
                    tentative.append(last)
                    tentative.append(next)
                    break
                else:
                    tentative.append(next)
            self.generators = tentative
        else:
            degree = order.degree
            # sort in descending order
            self.generators.sort(cmp=lambda a, b: cmp(degree(b), degree(a)))


class RationalFunctionField(ring.QuotientField):
    """
    The class for rational function fields.
    """
    _instances = {}

    def __init__(self, field, number_of_variables):
        """
        RationalFunctionField(field, number_of_variables)

        field: The field of coefficients.
        number_of_variables: The number of variables.
        """
        if isinstance(field, ring.QuotientField):
            basedomain = PolynomialRing(field.basedomain, number_of_variables)
        else:
            basedomain = PolynomialRing(field, number_of_variables)
        ring.QuotientField.__init__(self, basedomain)
        self.coefficient_field = field
        self.number_of_variables = number_of_variables

    def __repr__(self):
        """
        Return 'RationalFunctionField(Field, #vars)'
        """
        return "%s(%s, %d)" % (self.__class__.__name__, repr(self.coefficient_field), self.number_of_variables)

    def __str__(self):
        """
        Return K()()
        """
        return str(self.coefficient_field) + "()" * self.number_of_variables

    def __eq__(self, other):
        """
        equality test
        """
        if self is other:
            return True
        if not isinstance(other, RationalFunctionField):
            return False
        if (self.coefficient_field == other.coefficient_field and
            self.number_of_variables == other.number_of_variables):
            return True
        elif isinstance(self.coefficient_field, RationalFunctionField):
            return self.unnest() == other
        elif isinstance(other.coefficient_field, RationalFunctionField):
            return self == other.unnest()
        return False

    def __contains__(self, element):
        """
        Return True if an element is in the field.
        """
        if ring.getRing(element).issubring(self):
            return True
        return False

    def __hash__(self):
        """
        Return hash corresponding to equality.
        """
        return hash(self.coefficient_field) + self.number_of_variables

    def getQuotientField(self):
        """
        Return the quotient field (the field itself).
        """
        return self

    def getCoefficientRing(self):
        """
        Return the coefficient field.
        This method is provided for common interface with PolynomialRing.
        """
        return self.coefficient_field

    def issubring(self, other):
        """
        Return True if self is a subring of the other.
        """
        if self is other:
            return True
        elif (isinstance(other, RationalFunctionField) and
              self.coefficient_field.issubring(other.coefficient_field) and
              self.number_of_variables <= other.number_of_variables):
            return True
        return False

    def issuperring(self, other):
        """
        Return True if self is a superring of the other.
        """
        if self is other:
            return True
        elif (isinstance(other, RationalFunctionField) and
              self.coefficient_field.issuperring(other.coefficient_field) and
              self.number_of_variables >= other.number_of_variables):
            return True
        elif self.coefficient_field.issuperring(other):
            return True
        else:
            try:
                if self.issuperring(other.getQuotientField()):
                    return True
            except Exception:
                return False
            return False

    def getCharacteristic(self):
        """
        The characteristic of a rational function field is same as of
        its coefficient field.
        """
        return self.coefficient_field.getCharacteristic()

    def getCommonSuperring(self, other):
        """
        Return common superring.
        """
        if self.issuperring(other):
            return self
        elif other.issuperring(self):
            return other
        elif isinstance(other, (PolynomialRing, RationalFunctionField)):
            return RationalFunctionField(self.coefficient_field.getCommonSuperring(other.getCoefficientRing()), max(self.number_of_variables, other.number_of_variables))
        try:
            return other.getCommonSuperring(self)
        except RuntimeError:
            # reached recursion limit by calling on each other
            raise TypeError("no common super ring")
        except AttributeError:
            # other doesn't have getCommonSuperring
            raise TypeError("%s is not a ring" % str(other))

    def unnest(self):
        """
        if self is a nested RationalFunctionField i.e. its
        coefficientField is also a RationalFunctionField, then the
        function returns one level unnested RationalFunctionField.

        For example:
        RationalFunctionField(RationalFunctionField(Q, 1), 1).unnest()
        returns
        RationalFunctionField(Q, 2).
        """
        return RationalFunctionField(self.coefficient_field.coefficient_field, self.coefficient_field.number_of_variables + self.number_of_variables)

    def createElement(self, *seedarg, **seedkwd):
        """
        Return an element of the field made from seed.
        """
        import nzmath.poly.ratfunc as ratfunc
        return ratfunc.RationalFunction(*seedarg, **seedkwd)

    def _get_one(self):
        """
        getter for one
        """
        if self._one is None:
            import nzmath.poly.ratfunc as ratfunc
            poly_one = self.basedomain.one
            self._one = ratfunc.RationalFunction(poly_one, poly_one)
        return self._one

    one = property(_get_one, None, None, "multiplicative unit.")

    def _get_zero(self):
        "getter for zero"
        if self._zero is None:
            import nzmath.poly.ratfunc as ratfunc
            poly_one = self.basedomain.one
            poly_zero = self.basedomain.zero
            self._zero = ratfunc.RationalFunction(poly_zero, poly_one)
        return self._zero

    zero = property(_get_zero, None, None, "additive unit.")

    @classmethod
    def getInstance(cls, coefffield, number_of_variables):
        """
        Return an instance of the class with specified coefficient ring
        and number of variables.
        """
        if (coefffield, number_of_variables) not in cls._instances:
            cls._instances[coefffield, number_of_variables] = cls(coefffield, number_of_variables)
        return cls._instances[coefffield, number_of_variables]
