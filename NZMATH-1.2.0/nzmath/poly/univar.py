"""
Univariate Polynomials

The polynomials are immutable data types, under the public API.  If
one tries to manipulate its underlying data attributes, immutability
will be able to be broken.
"""

from __future__ import division
import warnings
import nzmath.ring as _ring
import nzmath.poly.formalsum as formalsum


class PolynomialInterface(formalsum.FormalSumContainerInterface):
    """
    Base class for all univariate polynomials.
    """
    def __init__(self, coefficients, **kwds):
        """
        For any descendants of PolynomialInterface, `coefficients'
        have to be given.  Other arguments, if necessary, should be
        passed as keyword arguments.  Note that in this __init__
        method, no actual initialization happens.
        """
        formalsum.FormalSumContainerInterface.__init__(self)
        self.number_of_variables = 1

    def __eq__(self, other):
        """
        self == other
        """
        if not isinstance(other, PolynomialInterface):
            warnings.warn("comparison falls back to that of formal sum.")
            return formalsum.FormalSumContainerInterface.__eq__(self, other)
        sterms = [t for t in self.iterterms() if t[1]]
        sterms.sort()
        oterms = [t for t in other.iterterms() if t[1]]
        oterms.sort()
        return sterms == oterms

    def __hash__(self):
       val = sum([hash(t) for t in self.iterterms() if t[1]])
       return val

    def __pow__(self, index):
        """
        self ** index
        """
        raise NotImplementedError("should be overridden")

    def ring_mul(self, other):
        """
        Multiplication of two polynomials in the same ring.
        """
        mul_coeff = {}
        for ds, cs in self:
            for do, co in other:
                term_degree = ds + do
                if term_degree in mul_coeff:
                    mul_coeff[term_degree] += cs*co
                else:
                    mul_coeff[term_degree] = cs*co
        return self.construct_with_default(sorted([(d, c) for (d, c) in mul_coeff.iteritems() if c]))

    def scalar_mul(self, scale):
        """
        Return the result of scalar multiplication.
        """
        return self.construct_with_default([(d, c * scale) for (d, c) in self if c])

    def term_mul(self, term):
        """
        Return the result of multiplication with the given term.
        The term can be given as a tuple (degree, coeff) or as a
        Polynomial instance.
        """
        if isinstance(term, PolynomialInterface):
            degree, coeff = iter(term).next()
        else:
            degree, coeff = term
        return self.construct_with_default([(d + degree, c * coeff) for (d, c) in self])

    def square(self):
        """
        Return square of this polynomial.
        """
        raise NotImplementedError("should be overridden")

    def differentiate(self):
        """
        Return the formal differentiation of self.
        """
        result = {}
        for d, c in self:
            if d > 0:
                dc = d * c
                if dc:
                    result[d - 1] = dc
        return self.construct_with_default(result)

    def upshift_degree(self, slide):
        """
        Return the polynomial obtained by shifting upward all terms
        with degrees of 'slide'.

        f.upshift_degree(slide) is equivalent to
        f.term_mul((slide, 1)).
        """
        return self.bases_map(lambda x: x + slide)

    def downshift_degree(self, slide):
        """
        Return the polynomial obtained by shifting downward all terms
        with degrees of 'slide'.

        Be careful that if the least degree term has the degree less
        than 'slide' then the result is not mathematically a
        polynomial.  Even in such a case, the method does not raise an
        exception.

        f.downshift_degree(slide) is equivalent to
        f.upshift_degree(-slide).
        """
        return self.bases_map(lambda x: x - slide)

    def terms_map(self, func):
        """
        Create a new formal sum container by applying func to each
        term.  func must be a function taking 2 arguments.
        """
        terms = []
        for t in self:
            b, c = func(*t)
            if c:
                terms.append((b, c))
        return self.construct_with_default(terms)

    def construct_with_default(self, terms):
        """
        Create a new univariate polynomial of the same class with
        self, with given only the terms and use copy of self's data if
        necessary.
        """
        return self.__class__(terms, **self._init_kwds)


class BasicPolynomial(PolynomialInterface):
    """
    Basic polynomial data type ignoring a variable name and the ring.
    """
    def __init__(self, coefficients, **kwds):
        """
        BasicPolynomial(coefficients)

        coefficients can be any dict initial values.
        """
        PolynomialInterface.__init__(self, coefficients, **kwds)
        self._coefficients = dict(coefficients)
        self._init_kwds = kwds

    def __add__(self, other):
        """
        self + other
        """
        sum_coeff = dict(iter(self))
        for term, coeff in other:
            if term in sum_coeff:
                sum_coeff[term] += coeff
            else:
                sum_coeff[term] = coeff
        return self.construct_with_default([(d, c) for (d, c) in sum_coeff.iteritems() if c])

    def __sub__(self, other):
        """
        self - other
        """
        dif_coeff = dict(iter(self))
        for term, coeff in other:
            if term in dif_coeff:
                dif_coeff[term] -= coeff
            else:
                dif_coeff[term] = -coeff
        return self.construct_with_default([(d, c) for (d, c) in dif_coeff.iteritems() if c])

    def __mul__(self, other):
        """
        self * other

        If type of other is Polynomial, do multiplication in ring.
        Otherwise, do scalar multiplication.
        """
        if isinstance(other, PolynomialInterface):
            return self.ring_mul(other)
        else:
            try:
                return self.scalar_mul(other)
            except TypeError:
                return NotImplemented

    def __rmul__(self, other):
        """
        other * self

        If type of other does not support multiplication with self
        from left, this method is called.  In the context, it is only
        posiible that other be a scalar.
        """
        try:
            return self.scalar_mul(other)
        except TypeError:
            return NotImplemented

    def __neg__(self):
        """
        -self
        """
        return self.construct_with_default([(d, -c) for (d, c) in self])

    def __pos__(self):
        """
        +self
        """
        return self.construct_with_default(self._coefficients)

    def square(self):
        """
        Return the square of self.
        """
        # zero
        if not self:
            return self
        data_length = len(self._coefficients)
        # monomial
        if data_length == 1:
            return self.construct_with_default([(d*2, c**2) for (d, c) in self])
        # binomial
        if data_length == 2:
            (d1, c1), (d2, c2) = self.terms()
            return self.construct_with_default({d1*2:c1**2, d1+d2:c1*c2*2, d2*2:c2**2})
        # general (inefficient)
        items = self._coefficients.items()
        fst, snd = {}, {}
        if data_length & 1:
            b, c = items.pop()
            fst[b] = c
        while items:
            b, c = items.pop()
            fst[b] = c
            b, c = items.pop()
            snd[b] = c
        fst = self.__class__(fst, **self._init_kwds)
        snd = self.__class__(snd, **self._init_kwds)
        mid = fst.ring_mul(snd.scalar_mul(2))
        return fst.square() + mid + snd.square()

    def __pow__(self, index):
        """
        self ** index
        """
        # special indices
        if index < 0:
            raise ValueError("negative index is not allowed.")
        elif index == 0:
            for c in self.itercoefficients():
                if c:
                    one = _ring.getRing(c).one
                    break
            else:
                one = 1
            return self.construct_with_default({0: one})
        elif index == 1:
            return self
        elif index == 2:
            return self.square()
        # special polynomials
        if not self:
            return self
        elif len(self._coefficients) == 1:
            return self.construct_with_default([(d*index, c**index) for (d, c) in self])
        # general
        power_product = self.construct_with_default({0: 1})
        power_of_2 = self
        while index:
            if index & 1:
                power_product *= power_of_2
            index //= 2
            if index:
                power_of_2 = power_of_2.square()
        return power_product

    def __call__(self, val):
        """
        substitution
        """
        items = self._coefficients.items()
        if not items:
            return 0*val
        d, c = items.pop()
        result = c * val**d
        for d, c in items:
            result += c * val**d
        return result

    def iterterms(self):
        """
        iterator for (degree, coefficient) pairs.
        The iterator is equivalent to
          zip(self.iterbases(), self.itercoefficients())
        """
        return self._coefficients.iteritems()

    def itercoefficients(self):
        """
        iterator for coefficients.
        """
        return self._coefficients.itervalues()

    def iterbases(self):
        """
        iterator for degrees.
        """
        return self._coefficients.iterkeys()

    def __getitem__(self, degree):
        """
        Return the coefficient of specified degree.
        If there is no term of degree, return 0.
        """
        return self._coefficients.get(degree, 0)

    def __contains__(self, degree):
        """
        Return True if there is a term of specified degree.
        False otherwise.
        """
        return degree in self._coefficients

    def __len__(self):
        """
        Return the number of data entries.
        """
        return len(self._coefficients)

    def __eq__(self, other):
        """
        self == other
        """
        if self is other:
            return True
        return PolynomialInterface.__eq__(self, other)

    def __hash__(self):
        return PolynomialInterface.__hash__(self)

    def __repr__(self): # for debug
        return "%s(%s)" % (self.__class__.__name__, repr(self._coefficients))


class SortedPolynomial (PolynomialInterface):
    """
    SortedPolynomial stores terms of a polynomial in sorted manner.
    All methods and operations keep it sorted.
    """
    def __init__(self, coefficients, _sorted=False, **kwds):
        """
        SortedPolynomial(coefficients)

        'coefficients' can be any dict initial values.
        Optionally '_sorted' can be True if the coefficients is
        already sorted.
        """
        PolynomialInterface.__init__(self, coefficients, **kwds)
        self.sorted = []
        if not _sorted or isinstance(coefficients, dict):
            for t in dict(coefficients).iteritems():
                self._insort(t)
        else:
            self.sorted = list(coefficients)
        self._init_kwds = kwds
        self._init_kwds['_sorted'] = True

    def _insort(self, term):
        """
        Insert the term into self.seorted list, and keep it sorted.

        - This method is destructive.
        - This method is not in a part of the API.
        """
        self.sorted.insert(self._bisect(term[0]), term)

    def _bisect(self, degree):
        """
        Return the index where to insert term of degree.

        - This method is not in a part of the API.
        - The code is just adapting bisect.bisect_right to the
          context.
        """
        lo, hi = 0, len(self.sorted)
        while lo < hi:
            mid = (lo + hi) >> 1
            if degree < self.sorted[mid][0]:
                hi = mid
            else:
                lo = mid + 1
        return lo

    def __pos__(self):
        """
        +self
        """
        return self.construct_with_default(self.sorted)

    def __neg__(self):
        """
        -self
        """
        return self.construct_with_default([(d, -c) for (d, c) in self])

    def __add__(self, other):
        """
        self + other
        """
        if self.sorted:
            iter_self = iter(self.sorted)
            self_term = iter_self.next()
        else:
            self_term = None
        if other.sorted:
            iter_other = iter(other.sorted)
            other_term = iter_other.next()
        else:
            other_term = None
        sorted = []
        while self_term and other_term:
            compared = cmp(self_term[0], other_term[0])
            if compared < 0:
                sorted.append(self_term)
                try:
                    self_term = iter_self.next()
                except StopIteration:
                    self_term = None
                    break
            elif compared == 0:
                c = self_term[1] + other_term[1]
                if c:
                    sorted.append((self_term[0], c))
                try:
                    self_term = iter_self.next()
                except StopIteration:
                    self_term = None
                try:
                    other_term = iter_other.next()
                except StopIteration:
                    other_term = None
                if self_term is None or other_term is None:
                    break
            else:
                sorted.append(other_term)
                try:
                    other_term = iter_other.next()
                except StopIteration:
                    other_term = None
                    break
        if other_term is None and self_term:
            sorted.append(self_term)
            for term in iter_self:
                sorted.append(term)
        elif self_term is None and other_term:
            sorted.append(other_term)
            for term in iter_other:
                sorted.append(term)
        return self.construct_with_default(sorted)

    def __sub__(self, other):
        """
        self - other

        It is something like merge sort.
        """
        if self.sorted:
            iter_self = iter(self.sorted)
            self_term = iter_self.next()
        else:
            self_term = None
        if other.sorted:
            iter_other = iter(other.sorted)
            other_term = iter_other.next()
        else:
            other_term = None
        sorted = []
        while self_term and other_term:
            compared = cmp(self_term[0], other_term[0])
            if compared < 0:
                sorted.append(self_term)
                try:
                    self_term = iter_self.next()
                except StopIteration:
                    self_term = None
                    break
            elif compared == 0:
                c = self_term[1] - other_term[1]
                if c:
                    sorted.append((self_term[0], c))
                try:
                    self_term = iter_self.next()
                except StopIteration:
                    self_term = None
                try:
                    other_term = iter_other.next()
                except StopIteration:
                    other_term = None
                if self_term is None or other_term is None:
                    break
            else:
                sorted.append((other_term[0], -other_term[1]))
                try:
                    other_term = iter_other.next()
                except StopIteration:
                    other_term = None
                    break
        if other_term is None and self_term:
            sorted.append(self_term)
            for term in iter_self:
                sorted.append(term)
        elif self_term is None and other_term:
            sorted.append((other_term[0], -other_term[1]))
            for term in iter_other:
                sorted.append((term[0], -term[1]))
        return self.construct_with_default(sorted)

    def __mul__(self, other):
        """
        self * other

        If type of other is SortedPolynomial, do multiplication
        in ring.  Otherwise, do scalar multiplication.
        """
        if isinstance(other, PolynomialInterface):
            return self.ring_mul(other)
        else:
            try:
                return self.scalar_mul(other)
            except TypeError:
                return NotImplemented

    def __rmul__(self, other):
        """
        other * self

        If type of other does not support multiplication with self
        from left, this method is called.  In the context, it is only
        posiible that other be a scalar.
        """
        try:
            return self.scalar_mul(other)
        except TypeError:
            return NotImplemented

    def ring_mul_karatsuba(self, other):
        """
        Multiplication of two polynomials in the same ring.

        Computation is carried out by Karatsuba method.
        """
        polynomial = self.construct_with_default
        # zero
        if not self or not other:
            return polynomial(())
        # monomial
        if len(self.sorted) == 1:
            return other.term_mul(self)
        if len(other.sorted) == 1:
            return self.term_mul(other)
        # binomial
        if len(self.sorted) == 2:
            p, q = [other.term_mul(term) for term in self]
            return p + q
        if len(other.sorted) == 2:
            p, q = [self.term_mul(term) for term in other]
            return p + q
        # suppose self is black and other is red.
        black_left_degree, black_right_degree = self.sorted[0][0], self.sorted[-1][0]
        red_left_degree, red_right_degree = other.sorted[0][0], other.sorted[-1][0]
        left_degree = min(black_left_degree, red_left_degree)
        right_degree = max(black_right_degree, red_right_degree)
        # we assert here that order is of ascending. (is it correct?)
        assert left_degree < right_degree
        half_degree = (left_degree + right_degree) >> 1
        black_half_index = self._bisect(half_degree)
        red_half_index = other._bisect(half_degree)
        if not black_half_index:
            return (self.downshift_degree(black_left_degree).ring_mul_karatsuba(other)).upshift_degree(black_left_degree)
        if not red_half_index:
            return (self.ring_mul_karatsuba(other.downshift_degree(red_left_degree))).upshift_degree(red_left_degree)
        club = polynomial([(d - left_degree, c) for (d, c) in self.sorted[:black_half_index]])
        spade = polynomial([(d - half_degree, c) for (d, c) in self.sorted[black_half_index:]])
        dia = polynomial([(d - left_degree, c) for (d, c) in other.sorted[:red_half_index]])
        heart = polynomial([(d - half_degree, c) for (d, c) in other.sorted[red_half_index:]])
        weaker = club.ring_mul_karatsuba(dia)
        stronger = spade.ring_mul_karatsuba(heart)
        karatsuba = (club + spade).ring_mul_karatsuba(dia + heart) - weaker - stronger
        if left_degree:
            return (weaker.upshift_degree(left_degree * 2) +
                    karatsuba.upshift_degree(left_degree + half_degree) +
                    stronger.upshift_degree(half_degree * 2))
        else:
            return (weaker +
                    karatsuba.upshift_degree(half_degree) +
                    stronger.upshift_degree(half_degree * 2))

    def scalar_mul(self, scale):
        """
        Return the result of scalar multiplication.
        """
        keep_ring = True
        if "coeffring" in self._init_kwds:
            new_coeff = []
            coeffring = self._init_kwds["coeffring"]
            for d, c in self:
                if c:
                    scaled = c * scale
                    if keep_ring and scaled not in coeffring:
                        coeffring = coeffring.getCommonSuperring(_ring.getRing(scaled))
                    new_coeff.append((d, scaled))
            self._init_kwds["coeffring"] = coeffring
        else:
            new_coeff = [(d, c * scale) for (d, c) in self if c]
        return self.construct_with_default(new_coeff)

    def square(self):
        """
        Return the square of self.
        """
        # zero
        if not self:
            return self

        polynomial = self.construct_with_default
        data_length = len(self.sorted)
        # monomial
        if data_length == 1:
            d, c = self.sorted[0]
            if d:
                return polynomial([(d*2, c**2)])
            else:
                return polynomial([(0, c**2)])
        # binomial
        if data_length == 2:
            (d1, c1), (d2, c2) = [(d, c) for (d, c) in self]
            return polynomial({d1*2:c1**2, d1+d2:c1*c2*2, d2*2:c2**2})
        # general (Karatsuba)
        right_degree = self.sorted[-1][0]
        left_degree = self.sorted[0][0]
        half_degree = (right_degree + left_degree) >> 1
        half_index = self._bisect(half_degree)
        fst = polynomial([(d - left_degree, c) for (d, c) in self.sorted[:half_index]])
        snd = polynomial([(d - half_degree, c) for (d, c) in self.sorted[half_index:]])
        fst_squared = fst.square()
        snd_squared = snd.square()
        karatsuba = (fst + snd).square() - fst_squared - snd_squared
        if left_degree:
            return (fst_squared.upshift_degree(left_degree * 2) +
                    karatsuba.upshift_degree(left_degree + half_degree) +
                    snd_squared.upshift_degree(half_degree * 2))
        else:
            return (fst_squared +
                    karatsuba.upshift_degree(half_degree) +
                    snd_squared.upshift_degree(half_degree * 2))

    def __pow__(self, index):
        """
        self ** index

        No ternary powering is defined here, because there is no
        modulus operator defined.
        """
        # special indices
        if index < 0:
            raise ValueError("negative index is not allowed.")
        elif index == 0:
            for c in self.itercoefficients():
                if c:
                    one = _ring.getRing(c).one
                    break
            else:
                one = 1
            return self.construct_with_default([(0, one)])
        elif index == 1:
            return self
        elif index == 2:
            return self.square()
        # special polynomials
        if not self:
            return self
        elif len(self.sorted) == 1:
            return self.construct_with_default([(d*index, c**index) for (d, c) in self])
        # general
        power_product = self.construct_with_default([(0, 1)])
        power_of_2 = self
        while index:
            if index & 1:
                power_product *= power_of_2
            index //= 2
            if index:
                power_of_2 = power_of_2.square()
        return power_product

    def degree(self):
        """
        Return the degree of self
        """
        assert self.sorted == sorted(self.sorted)
        for d, c in self.sorted[::-1]:
            if c:
                return d
        return -1

    def leading_coefficient(self):
        """
        Return the leading coefficient
        """
        assert self.sorted == sorted(self.sorted)
        for d, c in self.sorted[::-1]:
            if c:
                return c
        return 0

    def leading_term(self):
        """
        Return the leading term as a tuple (degree, coefficient).
        """
        assert self.sorted == sorted(self.sorted)
        for d, c in self.sorted[::-1]:
            if c:
                return (d, c)
        return (-1, 0)

    def iterterms(self):
        """
        iterator for (base, coefficient) pairs.
        The iterator is equivalent to
          zip(self.iterbases(), self.itercoefficients())
        """
        return iter(self.sorted)

    def iterbases(self):
        """
        iterator for bases.
        """
        for term in self.sorted:
            yield term[0]

    def itercoefficients(self):
        """
        iterator for coefficients.
        """
        for term in self.sorted:
            yield term[1]

    def __getitem__(self, degree):
        """
        Return the coefficient of specified degree.
        If there is no term of degree, return 0.
        """
        if self.sorted:
            rindex = self._bisect(degree)
            if self.sorted[rindex - 1][0] == degree:
                return self.sorted[rindex - 1][1]
        return 0

    def __contains__(self, degree):
        """
        Return True if there is a term of specified degree.
        False otherwise.
        """
        rindex = self._bisect(degree)
        if self.sorted[rindex - 1][0] == degree:
            return True
        return False

    def __len__(self):
        """
        Return the number of data entries.
        """
        return len(self.sorted)

    def __eq__(self, other):
        """
        self == other
        """
        if self is other:
            return True
        if not isinstance(other, SortedPolynomial):
            return PolynomialInterface.__eq__(self, other)
        if [t for t in self if t[1]] == [t for t in other if t[1]]:
            return True
        return False

    def __hash__(self):
        val = sum([hash(t) for t in self if t[1]])
        return val

    def __call__(self, val):
        """
        substitution
        """
        items = list(self.sorted)
        if not items:
            return 0 * val
        val_pow = {1:val}
        d, c = items.pop()
        result = c
        while len(items) > 0:
            (e, f) = items.pop()
            result = result * val_pow.setdefault(d - e, val ** (d - e)) + f
            d = e
        if d:
            return result * val_pow.get(d, val ** d)
        else:
            return result

    def __repr__(self): # debug
        return repr(self.sorted)
