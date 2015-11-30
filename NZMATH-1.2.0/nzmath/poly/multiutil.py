"""
Multivariate polynomial extenders.
"""

from __future__ import division
import nzmath.rational as rational
import nzmath.ring as ring
import nzmath.poly.termorder as termorder
import nzmath.poly.ring as poly_ring
import nzmath.poly.uniutil as uniutil
import nzmath.poly.multivar as multivar
import nzmath.poly.ratfunc as ratfunc


_MIXIN_MSG = "%s is mix-in"


class OrderProvider (object):
    """
    OrderProvider provides order and related operations.
    """
    def __init__(self, order):
        """
        Do not instantiate OrderProvider.
        This initializer should be called from descendant:
          OrderProvider.__init__(self, order)
        """
        if type(self) is OrderProvider:
            raise NotImplementedError(_MIXIN_MSG % self.__class__.__name__)
        self.order = order


class NestProvider (object):
    """
    Provide nest/unnest pair to convert a multivar polynomial to a
    univar polynomial of polynomial coefficient and opposite
    direction.
    """
    def leading_variable(self):
        """
        Return the position of the leading variable (the leading term
        among all total degree one terms).

        The leading term varies with term orders, so does the result.
        The term order can be specified via the attribute 'order'.
        """
        if hasattr(self, 'order'):
            order = self.order
            pivot, pvar = (1,) + (0,)*(self.number_of_variables - 1), 0
            for var in range(1, self.number_of_variables):
                vindex = (0,)*var + (1,) + (0,)*(self.number_of_variables - var - 1)
                if order.cmp(pivot, vindex) < 0:
                    pivot, pvar = vindex, var
        else:
            pvar = 0
        return pvar

    def nest(self, outer, coeffring):
        """
        Nest the polynomial by extracting outer variable at the given
        position.
        """
        combined = {}
        if self.number_of_variables == 2:
            itercoeff = lambda coeff: [(i[0], c) for (i, c) in coeff]
            poly = uniutil.polynomial
            polyring = poly_ring.PolynomialRing.getInstance(coeffring)
        elif self.number_of_variables >= 3:
            itercoeff = lambda coeff: coeff
            poly = self.__class__
            polyring = poly_ring.PolynomialRing.getInstance(coeffring, self.number_of_variables - 1)
        else:
            raise TypeError("number of variable is not multiple")
        for index, coeff in self.combine_similar_terms(outer):
            combined[index] = poly(itercoeff(coeff), coeffring=coeffring)
        return uniutil.polynomial(combined, coeffring=polyring)

    def unnest(self, q, outer, coeffring):
        """
        Unnest the nested polynomial q by inserting outer variable at
        the given position.
        """
        q_coeff = {}
        for d, cpoly in q:
            for inner_d, inner_c in cpoly:
                if isinstance(inner_d, (int, long)):
                    inner_d = [inner_d]
                else:
                    inner_d = list(inner_d)
                inner_d.insert(outer, d)
                q_coeff[tuple(inner_d)] = inner_c
        return self.__class__(q_coeff, coeffring=coeffring)


class RingElementProvider (ring.CommutativeRingElement):
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

    def getCoefficientRing(self):
        """
        Return the coefficient ring.
        """
        return self._coefficient_ring

    def set_coefficient_ring(self, coeffring):
        if self._coefficient_ring is None:
            self._coefficient_ring = coeffring
            self._ring = poly_ring.PolynomialRing.getInstance(self._coefficient_ring, self.number_of_variables)


class PseudoDivisionProvider (object):
    """
    PseudoDivisionProvider provides pseudo divisions for multivariate
    polynomials.  It is assumed that the coefficient ring of the
    polynomials is a domain.

    The class should be used with NestProvider, RingElementProvider.
    """
    def pseudo_divmod(self, other):
        """
        self.pseudo_divmod(other) -> (Q, R)

        Q, R are polynomials such that
        d**(deg(self) - deg(other) + 1) * self == other * Q + R,
        w.r.t. a fixed variable, where d is the leading coefficient of
        other.

        The leading coefficient varies with term orders, so does the
        result.  The term order can be specified via the attribute
        'order'.
        """
        var = self.leading_variable()
        coeffring = self.getCoefficientRing()

        s = self.nest(var, coeffring)
        o = other.nest(var, coeffring)
        q, r = s.pseudo_divmod(o)
        qpoly = self.unnest(q, var, coeffring)
        rpoly = self.unnest(r, var, coeffring)
        return qpoly, rpoly

    def pseudo_floordiv(self, other):
        """
        self.pseudo_floordiv(other) -> Q

        Q is a polynomial such that
        d**(deg(self) - deg(other) + 1) * self == other * Q + R,
        where d is the leading coefficient of other and R is a
        polynomial.

        The leading coefficient varies with term orders, so does the
        result.  The term order can be specified via the attribute
        'order'.
        """
        var = self.leading_variable()
        coeffring = self.getCoefficientRing()

        s = self.nest(var, coeffring)
        o = other.nest(var, coeffring)
        q = s.pseudo_floordiv(o)
        return self.unnest(q, var, coeffring)

    def pseudo_mod(self, other):
        """
        self.pseudo_mod(other) -> R

        R is a polynomial such that
        d**(deg(self) - deg(other) + 1) * self == other * Q + R,
        where d is the leading coefficient of other and Q a
        polynomial.

        The leading coefficient varies with term orders, so does the
        result.  The term order can be specified via the attribute
        'order'.
        """
        var = self.leading_variable()
        coeffring = self.getCoefficientRing()

        s = self.nest(var, coeffring)
        o = other.nest(var, coeffring)
        r = s.pseudo_mod(o)
        return self.unnest(r, var, coeffring)

    def __truediv__(self, other):
        """
        self / other

        The result is a rational function.
        """
        return ratfunc.RationalFunction(self, other)

    def exact_division(self, other):
        """
        Return quotient of exact division.
        """
        coeffring = self.getCoefficientRing()
        if other in coeffring:
            new_coeffs = []
            keep_ring = True
            for i, c in self:
                ratio = c / other
                if keep_ring and ratio not in coeffring:
                    keep_ring = False
                    new_coeffring = ratio.getRing()
                new_coeffs.append((i, ratio))
            if keep_ring:
                return self.__class__(new_coeffs, coeffring=coeffring)
            else:
                return self.__class__(new_coeffs, coeffring=new_coeffring)

        var = self.leading_variable()
        coeffring = self.getCoefficientRing()
        s = self.nest(var, coeffring)
        o = other.nest(var, coeffring)
        q = s.exact_division(o)
        return self.unnest(q, var, coeffring)


class GcdProvider (object):
    """
    Provides greatest common divisor for multivariate polynomial.

    The class should be used with NestProvider, RingElementProvider.
    """
    def gcd(self, other):
        """
        Return gcd.
        The nested polynomials' gcd is used.
        """
        var = self.leading_variable()
        coeffring = self.getCoefficientRing()
        s = self.nest(var, coeffring)
        o = other.nest(var, coeffring)
        if hasattr(s, "gcd"):
            g = s.gcd(o)
        elif hasattr(s, "subresultant_gcd"):
            g = s.subresultant_gcd(o)
        else:
            raise TypeError("no gcd method available")
        return self.unnest(g, var, coeffring)


class RingPolynomial (OrderProvider,
                      NestProvider,
                      multivar.BasicPolynomial,
                      RingElementProvider):
    """
    General polynomial with commutative ring coefficients.
    """
    def __init__(self, coefficients, **kwds):
        """
        Initialize the polynomial.

        Required argument:
        - coefficients: initializer for polynomial coefficients

        Keyword arguments should include:
        - coeffring: domain
        - number_of_variables: the number of variables
        """
        if "number_of_variables" not in kwds:
            coefficients = dict(coefficients)
            for i in coefficients.iterkeys():
                kwds["number_of_variables"] = len(i)
                break
        multivar.BasicPolynomial.__init__(self, coefficients, **kwds)
        NestProvider.__init__(self)
        PseudoDivisionProvider.__init__(self)
        RingElementProvider.__init__(self)
        coeffring = None
        if "coeffring" in kwds:
            coeffring = kwds["coeffring"]
        else:
            coeffring = uniutil.init_coefficient_ring(self._coefficients)
        if coeffring is not None:
            self.set_coefficient_ring(coeffring)
        else:
            raise TypeError("argument `coeffring' is required")
        if "order" in kwds:
            order = kwds["order"]
        else:
            order = termorder.lexicographic_order
        OrderProvider.__init__(self, order)

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

    def __repr__(self): # debug
        return "%s(%s)" % (self.__class__.__name__, self._coefficients)

    def __add__(self, other):
        try:
            return multivar.BasicPolynomial.__add__(self, other)
        except (AttributeError, TypeError):
            one = self.getRing().one
            try:
                return multivar.BasicPolynomial.__add__(self, other * one)
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
            return multivar.BasicPolynomial.__sub__(self, other)
        except (AttributeError, TypeError):
            one = self.getRing().one
            try:
                return multivar.BasicPolynomial.__sub__(self, other * one)
            except Exception:
                return NotImplemented

    def __rsub__(self, other):
        one = self.getRing().one
        try:
            return other * one - self
        except Exception:
            return NotImplemented


class DomainPolynomial (PseudoDivisionProvider,
                        RingPolynomial):
    """
    Polynomial with domain coefficients.
    """
    def __init__(self, coefficients, **kwds):
        """
        Initialize the polynomial.

        - coefficients: initializer for polynomial coefficients
        - coeffring: domain
        """
        RingPolynomial.__init__(self, coefficients, **kwds)
        if not self._coefficient_ring.isdomain():
            raise TypeError("coefficient ring has to be a domain")
        PseudoDivisionProvider.__init__(self)


class UniqueFactorizationDomainPolynomial (GcdProvider,
                                           DomainPolynomial):
    """
    Polynomial with unique factorization domain coefficients.
    """
    def __init__(self, coefficients, **kwds):
        """
        Initialize the polynomial.

        - coefficients: initializer for polynomial coefficients
        - coeffring: unique factorization domain
        """
        DomainPolynomial.__init__(self, coefficients, **kwds)
        if not self._coefficient_ring.isufd():
            raise TypeError("coefficient ring has to be a UFD")
        GcdProvider.__init__(self)

    def resultant(self, other, var):
        """
        Return resultant of two polynomials of the same ring, with
        respect to the variable specified by its position var.
        """
        cring = self._coefficient_ring
        return self.nest(var, cring).resultant(other.nest(var, cring))


class PolynomialRingAnonymousVariables (ring.CommutativeRing):
    """
    The class of multivariate polynomial ring.
    There's no need to specify the variable names.
    """

    _instances = {}

    def __init__(self, coeffring, number_of_variables):
        if not isinstance(coeffring, ring.Ring):
            raise TypeError("%s should not be passed as ring" % coeffring.__class__)
        ring.CommutativeRing.__init__(self)
        self._coefficient_ring = coeffring
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
        getQuotientField returns the quotient field of the ring
        if coefficient ring has its quotient field.  Otherwise,
        an exception will be raised.
        """
        coefficientField = self._coefficient_ring.getQuotientField()
        variables = ["x%d" % i for i in range(self.number_of_variables)]
        return ratfunc.RationalFunctionField(coefficientField, variables)

    def __eq__(self, other):
        if self is other:
            return True
        if (isinstance(other, PolynomialRingAnonymousVariables) and
            self._coefficient_ring == other._coefficient_ring and
            self.number_of_variables == other.number_of_variables):
            return True
        return False

    def __repr__(self):
        """
        Return 'PolynomialRingAnonymousVariables(Ring, #vars)'
        """
        return "%s(%s, %d)" % (self.__class__.__name__, repr(self._coefficient_ring), self.number_of_variables)

    def __str__(self):
        """
        Return R[][]
        """
        return str(self._coefficient_ring) + "[]" * self.number_of_variables

    def __hash__(self):
        """
        hash(self)
        """
        return (hash(self._coefficient_ring) ^ (self.number_of_variables * hash(self.__class__.__name__) + 1)) & 0x7fffffff

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
        if isinstance(other, poly_ring.PolynomialRing):
            if (self._coefficient_ring.issubring(other.getCoefficientRing()) and
                self.number_of_variables <= other.number_of_variables):
                return True
        elif isinstance(other, poly_ring.RationalFunctionField):
            if (len(other.vars) >= self.number_of_variables and
                other.coefficientField.issuperring(self._coefficient_ring)):
                return True
        try:
            return other.issuperring(self)
        except RuntimeError:
            # reach max recursion by calling each other
            return False

    def issuperring(self, other):
        """
        reports whether this polynomial ring contains another ring.
        """
        if self._coefficient_ring.issuperring(other):
            return True
        if isinstance(other, poly_ring.PolynomialRing):
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
        elif (not isinstance(other, PolynomialRingAnonymousVariables) and
              other.issuperring(self._coefficient_ring)):
            return self.__class__(other, self.number_of_variables)
        try:
            if hasattr(other, "getCommonSuperring"):
                return other.getCommonSuperring(self)
        except RuntimeError:
            # reached recursion limit by calling on each other
            pass
        raise TypeError("no common super ring")

    def createElement(self, seed):
        """
        Return an element of the polynomial ring made from seed
        overriding ring.createElement.
        """
        if not seed:
            return polynomial((), coeffring=self._coefficient_ring, number_of_variables=self.number_of_variables)
        elif seed in self._coefficient_ring:
            return polynomial([((0,)*self.number_of_variables, seed)], coeffring=self._coefficient_ring)
        # implementation should be replaced later
        raise NotImplementedError("unclear which type of polynomial be chosen")

    def _getOne(self):
        "getter for one"
        if self._one is None:
            self._one = self.createElement(self._coefficient_ring.one)
        return self._one

    one = property(_getOne, None, None, "multiplicative unit")

    def _getZero(self):
        "getter for zero"
        if self._zero is None:
            self._zero = self.createElement(self._coefficient_ring.zero)
        return self._zero

    zero = property(_getZero, None, None, "additive unit")

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
    def getInstance(cls, coeffring, number_of_variables):
        """
        Return an instance of the class with specified coefficient ring
        and number of variables.
        """
        if (coeffring, number_of_variables) not in cls._instances:
            cls._instances[coeffring, number_of_variables] = cls(coeffring, number_of_variables)
        return cls._instances[coeffring, number_of_variables]


class PolynomialIdeal (ring.Ideal):
    """
    Multivariate polynomial ideal.
    """
    def __init__(self, generators, aring):
        """
        Initialize a polynomial ideal.
        """
        ring.Ideal.__init__(self, generators, aring)

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
        return "(%s)%s" % (", ".join([str(g) for g in self.generators]), self.ring)


# factories

special_ring_table = {}

def polynomial(coefficients, coeffring, number_of_variables=None):
    """
    Return a polynomial.
    - coefficients has to be a initializer for dict, whose keys are
      variable indices and values are coefficients at the indices.
    - coeffring has to be an object inheriting ring.Ring.
    - number_of_variables has to be the number of variables.

    One can override the way to choose a polynomial type from a
    coefficient ring, by setting:
    special_ring_table[coeffring_type] = polynomial_type
    before the function call.
    """
    if type(coeffring) in special_ring_table:
        poly_type = special_ring_table[type(coeffring)]
    elif coeffring.isufd():
        poly_type = UniqueFactorizationDomainPolynomial
    elif coeffring.isdomain():
        poly_type = DomainPolynomial
    else:
        poly_type = multivar.BasicPolynomial
    if number_of_variables is None:
        coefficients = dict(coefficients)
        for k in coefficients:
            number_of_variables = len(k)
            break
    return poly_type(coefficients, coeffring=coeffring, number_of_variables=number_of_variables)

def MultiVariableSparsePolynomial(coefficient, variable, coeffring=None):
    """
    MultiVariableSparsePolynomial(coefficient, variable [,coeffring])

    - coefficient has to be a dictionary of form {(i1,...,ik): c}
    - variable has to be a list of character strings.
    - coeffring has to be, if specified, an object inheriting ring.Ring.
    
    This function is provided for backward compatible way of defining
    multivariate polynomial.  The variable names are ignored, but
    their number is used.
    """
    if not isinstance(variable, list) or not isinstance(coefficient, dict):
        raise ValueError("You must input MultiVariableSparsePolynomial(dict, list) but (%s, %s)." % (coefficient.__class__, variable.__class__))
    if coeffring is None:
        coeffring = uniutil.init_coefficient_ring(coefficient)
    return polynomial(coefficient, coeffring=coeffring, number_of_variables=len(variable))

def prepare_indeterminates(names, ctx, coeffring=None):
    """
    From space separated names of indeterminates, prepare variables
    representing the indeterminates.  The result will be stored in ctx
    dictionary.

    The variables should be prepared at once, otherwise wrong aliases
    of variables may confuse you in later calculation.

    If an optional coeffring is not given, indeterminates will be
    initialized as integer coefficient polynomials.

    Example:
    >>> prepare_indeterminates("X Y Z", globals())
    >>> Y
    UniqueFactorizationDomainPolynomial({(0, 1, 0): 1})

    """
    split_names = names.split()
    number_of_variables = len(split_names)
    if coeffring is None:
        coeffring = uniutil.init_coefficient_ring({1:1})
    for i, name in enumerate(split_names):
        e_i = tuple([0] * i + [1] + [0] * (number_of_variables - i - 1))
        ctx[name] = polynomial({e_i: 1}, coeffring, number_of_variables)

