"""
finite fields.
"""

from __future__ import division
import logging
import nzmath.gcd as gcd
import nzmath.bigrandom as bigrandom
import nzmath.bigrange as bigrange
import nzmath.arith1 as arith1
import nzmath.prime as prime
import nzmath.ring as ring
import nzmath.rational as rational
import nzmath.factor.misc as factor_misc
import nzmath.intresidue as intresidue
import nzmath.poly.uniutil as uniutil
import nzmath.matrix as matrix
import nzmath.vector as vector
import nzmath.compatibility

_log = logging.getLogger('nzmath.finitefield')


class FiniteFieldElement(ring.FieldElement):
    """
    The base class for all finite field element.
    """
    def __init__(self):
        # This class is abstract and can not be instantiated.
        if type(self) is FiniteFieldElement:
            raise NotImplementedError("the class FiniteFieldElement is abstract")
        ring.FieldElement.__init__(self)


class FiniteField(ring.Field):
    """
    The base class for all finite fields.
    """
    def __init__(self, characteristic):
        # This class is abstract and can not be instantiated.
        if type(self) is FiniteField:
            raise NotImplementedError("the class FiniteField is abstract")
        ring.Field.__init__(self)
        self.char = characteristic
        self._orderfactor = None  # postpone the initialization

    def card(self):
        "Cardinality of the field"
        raise NotImplementedError("card have to be overridden")

    def getCharacteristic(self):
        """
        Return the characteristic of the field.
        """
        return self.char

    def order(self, elem):
        """
        Find and return the order of the element in the multiplicative
        group of the field.
        """
        if not elem:
            raise ValueError("zero is not in the multiplicative group.")
        if self._orderfactor is None:
            self._orderfactor = factor_misc.FactoredInteger(card(self) - 1)
        o = 1
        for p, e in self._orderfactor:
            b = elem ** (int(self._orderfactor) // (p**e))
            while b != self.one:
                o = o * p
                b = b ** p
        return o

    def random_element(self, *args):
        """
        Return a randomly chosen element og the field.
        """
        return self.createElement(bigrandom.randrange(*args))

    def primitive_element(self):
        """
        Return a primitive element of the field, i.e., a generator of
        the multiplicative group.
        """
        raise NotImplementedError("primitive_element should be overridden")

    def Legendre(self, element):
        """ Return generalize Legendre Symbol for FiniteField.
        """
        if not element:
            return 0

        if element == self.one or self.char == 2:
            return 1 # element of cyclic group order 2**n-1 also always 1

        # generic method:successive squaring
        # generalized Legendre symbol definition:
        #    (self/_ring) := self ** ((card(_ring)-1)/2)
        x = element ** ((card(self)-1) >> 1)
        if x == self.one:
            return 1
        elif x == -self.one:
            return -1
        raise ValueError("element must be not in field")

    def TonelliShanks(self, element):
        """ Return square root of element if exist.
        assume that characteristic have to be more than three.
        """
        if  self.char == 2:
            return self.sqrt(element) # should be error

        if self.Legendre(element) == -1:
            raise ValueError("There is no solution")

        # symbol and code reference from Cohen, CCANT 1.5.1
        (e, q) = arith1.vp(card(self)-1, 2)

        a = element
        n = self.createElement(self.char+1)
        while self.Legendre(n) != -1:
            n = self.random_element(2, card(self)) # field maybe large
        y = z = n ** q
        r = e
        x = a ** ((q-1) >> 1)
        b = a * (x ** 2)
        x = a * x
        while True:
            if b == self.one:
                return x
            m = 1
            while m < r:
                if b ** (2 ** m) == self.one:
                    break
                m = m+1
            if m == r:
                break
            t = y ** (2 ** (r-m-1))
            y = t ** 2
            r = m
            x = x * t
            b = b * y
        raise ValueError("There is no solution")

    def sqrt(self, element):
        """ Return square root if exist.
        """
        if not element or element == self.one:
            return element # trivial case

        # element of characteristic 2 always exist square root
        if  self.char == 2:
            return element ** ((card(self)) >> 1)

        # otherwise,
        return self.TonelliShanks(element)


class FinitePrimeFieldElement(intresidue.IntegerResidueClass, FiniteFieldElement):
    """
    The class for finite prime field element.
    """
    def __init__(self, representative, modulus, modulus_is_prime=True):
        if not modulus_is_prime and not prime.primeq(abs(modulus)):
            raise ValueError("modulus must be a prime.")

        FiniteFieldElement.__init__(self)
        intresidue.IntegerResidueClass.__init__(self, representative, modulus)

        # ring
        self.ring = None

    def __repr__(self):
        return "FinitePrimeFieldElement(%d, %d)" % (self.n, self.m)

    def __str__(self):
        return "%d in F_%d" % (self.n, self.m)

    def getRing(self):
        """
        Return the finite prime field to which the element belongs.
        """
        if self.ring is None:
            self.ring = FinitePrimeField.getInstance(self.m)
        return self.ring

    def order(self):
        """
        Find and return the order of the element in the multiplicative
        group of F_p.
        """
        if self.n == 0:
            raise ValueError("zero is not in the group.")
        if not hasattr(self, "orderfactor"):
            self.orderfactor = factor_misc.FactoredInteger(self.m - 1)
        o = 1
        for p, e in self.orderfactor:
            b = self ** (int(self.orderfactor) // (p**e))
            while b.n != 1:
                o = o*p
                b = b**p
        return o


class FinitePrimeField(FiniteField):
    """
    FinitePrimeField is also known as F_p or GF(p).
    """

    # class variable
    _instances = {}

    def __init__(self, characteristic):
        FiniteField.__init__(self, characteristic)
        self.registerModuleAction(rational.theIntegerRing, self._int_times)
        # mathematically Q_p = Q \ {r/s; gcd(r, s) == 1, gcd(s, p) > 1}
        # is more appropriate.
        self.registerModuleAction(rational.theRationalField, self._rat_times)

    @staticmethod
    def _int_times(integer, fpelem):
        """
        Return k * FinitePrimeFieldElement(n, p)
        """
        return FinitePrimeFieldElement(integer * fpelem.n, fpelem.m)

    @staticmethod
    def _rat_times(rat, fpelem):
        """
        Return Rational(a, b) * FinitePrimeFieldElement(n, p)
        """
        return FinitePrimeFieldElement(rat * fpelem.n, fpelem.m)

    def __eq__(self, other):
        if self is other:
            return True
        if isinstance(other, FinitePrimeField):
            return self.char == other.char
        return False

    def __hash__(self):
        return self.char

    def __ne__(self, other):
        return not (self == other)

    def __str__(self):
        return "F_%d" % self.char

    def __repr__(self):
        return "%s(%d)" % (self.__class__.__name__, self.char)

    def __hash__(self):
        return self.char & 0xFFFFFFFF

    def issubring(self, other):
        """
        Report whether another ring contains the field as a subring.
        """
        if self == other:
            return True
        if isinstance(other, FiniteField) and other.getCharacteristic() == self.char:
            return True
        try:
            return other.issuperring(self)
        except:
            return False

    def issuperring(self, other):
        """
        Report whether the field is a superring of another ring.
        Since the field is a prime field, it can be a superring of
        itself only.
        """
        return self == other

    def __contains__(self, elem):
        if isinstance(elem, FinitePrimeFieldElement) and elem.getModulus() == self.char:
            return True
        return False

    def createElement(self, seed):
        """
        Create an element of the field.

        'seed' should be an integer.
        """
        return FinitePrimeFieldElement(seed, self.char, modulus_is_prime=True)

    def primitive_element(self):
        """
        Return a primitive element of the field, i.e., a generator of
        the multiplicative group.
        """
        if self.char == 2:
            return self.one
        fullorder = card(self) - 1
        if self._orderfactor is None:
            self._orderfactor = factor_misc.FactoredInteger(fullorder)
        for i in bigrange.range(2, self.char):
            g = self.createElement(i)
            for p in self._orderfactor.prime_divisors():
                if g ** (fullorder // p) == self.one:
                    break
            else:
                return g

    def Legendre(self, element):
        """ Return generalize Legendre Symbol for FinitePrimeField.
        """
        if not element:
            return 0
        if element.n == 1 or element.m == 2:
            return 1 # trivial

        return arith1.legendre(element.n, element.m)

    def SquareRoot(self, element):
        """ Return square root if exist.
        """
        if not element or element.n == 1:
            return element # trivial case
        if element.m == 2:
            return element.getRing().one
        return arith1.modsqrt(element.n, element.m)

    def card(self):
        "Cardinality of the field"
        return self.char

    # properties
    def _getOne(self):
        "getter for one"
        if self._one is None:
            self._one = FinitePrimeFieldElement(1, self.char)
        return self._one

    one = property(_getOne, None, None, "multiplicative unit.")

    def _getZero(self):
        "getter for zero"
        if self._zero is None:
            self._zero = FinitePrimeFieldElement(0, self.char)
        return self._zero

    zero = property(_getZero, None, None, "additive unit.")

    @classmethod
    def getInstance(cls, characteristic):
        """
        Return an instance of the class with specified characteristic.
        """
        if characteristic not in cls._instances:
            cls._instances[characteristic] = cls(characteristic)
        return cls._instances[characteristic]


FinitePrimeFieldPolynomial = uniutil.FinitePrimeFieldPolynomial
uniutil.special_ring_table[FinitePrimeField] = FinitePrimeFieldPolynomial


class ExtendedFieldElement(FiniteFieldElement):
    """
    ExtendedFieldElement is a class for an element of F_q.
    """
    def __init__(self, representative, field):
        """
        FiniteExtendedFieldElement(representative, field) creates
        an element of the finite extended field F_q.

        The argument representative has to be a polynomial with
        base field coefficients, i.e., if F_q is over F_{p^k}
        the representative has to F_{p^k} polynomial.

        Another argument field mut be an instance of
        ExtendedField.
        """
        if isinstance(field, ExtendedField):
            self.field = field
        else:
            raise TypeError("wrong type argument for field.")
        if (isinstance(representative, uniutil.FiniteFieldPolynomial) and
            representative.getCoefficientRing() == self.field.basefield):
            self.rep = self.field.modulus.mod(representative)
        else:
            _log.error(representative.__class__.__name__)
            raise TypeError("wrong type argument for representative.")

    def getRing(self):
        """
        Return the field to which the element belongs.
        """
        return self.field

    def _op(self, other, op):
        """
        Do `self (op) other'.
        op must be a name of the special method for binary operation.
        """
        if isinstance(other, ExtendedFieldElement):
            if other.field is self.field:
                result = self.field.modulus.mod(getattr(self.rep, op)(other.rep))
                return self.__class__(result, self.field)
            else:
                fq1, fq2 = self.field, other.field
                emb1, emb2 = double_embeddings(fq1, fq2)
                return getattr(emb1(self), op)(emb2(other))
        if self.field.hasaction(ring.getRing(other)):
            # cases for action ring elements
            embedded = self.field.getaction(ring.getRing(other))(other, self.field.one)
            result = self.field.modulus.mod(getattr(self.rep, op)(embedded.rep))
            return self.__class__(result, self.field)
        else:
            return NotImplemented

    def __add__(self, other):
        """
        self + other

        other can be an element of either F_q, F_p or Z.
        """
        if other is 0 or other is self.field.zero:
            return self
        return self._op(other, "__add__")

    __radd__ = __add__

    def __sub__(self, other):
        """
        self - other

        other can be an element of either F_q, F_p or Z.
        """
        if other is 0 or other is self.field.zero:
            return self
        return self._op(other, "__sub__")

    def __rsub__(self, other):
        """
        other - self

        other can be an element of either F_q, F_p or Z.
        """
        return self._op(other, "__rsub__")

    def __mul__(self, other):
        """
        self * other

        other can be an element of either F_q, F_p or Z.
        """
        if other is 1 or other is self.field.one:
            return self
        return self._op(other, "__mul__")

    __rmul__ = __mul__

    def __truediv__(self, other):
        return self * other.inverse()

    __div__ = __truediv__

    def inverse(self):
        """
        Return the inverse of the element.
        """
        if not self:
            raise ZeroDivisionError("There is no inverse of zero.")
        return self.__class__(self.rep.extgcd(self.field.modulus)[0], self.field)

    def __pow__(self, index):
        """
        self ** index

        pow() with three arguments is not supported.
        """
        if not self:
            return self # trivial
        if index < 0:
            index %= self.field.getCharacteristic()
        power = pow(self.rep, index, self.field.modulus)
        return self.__class__(power, self.field)

    def __neg__(self):
        return self.field.zero - self

    def __pos__(self):
        return self

    def __eq__(self, other):
        try:
            if self.field == other.field:
                if self.rep == other.rep:
                    return True
        except AttributeError:
            pass
        return False

    def __hash__(self):
        return hash(self.rep)

    def __ne__(self, other):
        return not (self == other)

    def __nonzero__(self):
        return bool(self.rep)

    def __repr__(self):
        return "%s(%s, %s)" % (self.__class__.__name__, repr(self.rep), repr(self.field))

    def trace(self):
        """
        Return the absolute trace.
        """
        p = self.field.char
        q = p
        alpha = self
        tr = alpha.rep[0]
        while q < card(self.field):
            q *= p
            alpha **= p
            tr += alpha.rep[0]
        if tr not in FinitePrimeField.getInstance(p):
            tr = FinitePrimeField.getInstance(p).createElement(tr)
        return tr

    def norm(self):
        """
        Return the absolute norm.
        """
        p = self.field.char
        nrm = (self ** ((card(self.field) - 1) // (p - 1))).rep[0]
        if nrm not in FinitePrimeField.getInstance(p):
            nrm = FinitePrimeField.getInstance(p).createElement(nrm)
        return nrm


class ExtendedField(FiniteField):
    """
    ExtendedField is a class for finite field, whose cardinality
    q = p**n with a prime p and n>1. It is usually called F_q or GF(q).
    """
    def __init__(self, basefield, modulus):
        """
        ExtendedField(basefield, modulus)

        Create a field extension basefield[X]/(modulus(X)).

        The modulus has to be an irreducible polynomial with
        coefficients in the basefield.
        """
        FiniteField.__init__(self, basefield.char)
        self.basefield = basefield
        self.modulus = modulus
        if isinstance(self.basefield, FinitePrimeField):
            self.degree = self.modulus.degree()
        else:
            self.degree = self.basefield.degree * self.modulus.degree()
            bf = self.basefield
            while hasattr(bf, "basefield"):
                # all basefields can be used as a scalar field.
                self.registerModuleAction(bf.basefield, self._scalar_mul)
                bf = bf.basefield
        self.registerModuleAction(self.basefield, self._scalar_mul)
        # integer is always a good scalar
        self.registerModuleAction(rational.theIntegerRing, self._scalar_mul)

    @classmethod
    def prime_extension(cls, characteristic, n_or_modulus):
        """
        ExtendedField.prime_extension(p, n_or_modulus) creates a
        finite field extended over prime field.

        characteristic must be prime. n_or_modulus can be:
          1) an integer greater than 1, or
          2) a polynomial in a polynomial ring of F_p with degree
             greater than 1.
        """
        if isinstance(n_or_modulus, (int, long)):
            n = n_or_modulus
            if n <= 1:
                raise ValueError("degree of extension must be > 1.")
            # choose a method among three variants:
            #modulus = cls._random_irriducible(characteristic, n)
            #modulus = cls._small_irriducible(characteristic, n)
            modulus = cls._primitive_polynomial(characteristic, n)
        elif isinstance(n_or_modulus, FinitePrimeFieldPolynomial):
            modulus = n_or_modulus
            if not isinstance(modulus.getCoefficientRing(), FinitePrimeField):
                raise TypeError("modulus must be F_p polynomial.")
            if modulus.degree() <= 1 or not modulus.isirreducible():
                raise ValueError("modulus must be of degree greater than 1.")
        else:
            raise TypeError("degree or modulus must be supplied.")
        return cls(FinitePrimeField.getInstance(characteristic), modulus)

    @staticmethod
    def _random_irriducible(char, degree):
        """
        Return randomly chosen irreducible polynomial of self.degree.
        """
        cardinality = char ** degree
        basefield = FinitePrimeField.getInstance(char)
        seed = bigrandom.randrange(1, char) + cardinality
        cand = uniutil.FiniteFieldPolynomial(enumerate(arith1.expand(seed, cardinality)), coeffring=basefield)
        while cand.degree() < degree or not cand.isirreducible():
            seed = bigrandom.randrange(1, cardinality) + cardinality
            cand = uniutil.FiniteFieldPolynomial(enumerate(arith1.expand(seed, cardinality)), coeffring=basefield)
        _log.debug(cand.order.format(cand))
        return cand

    @staticmethod
    def _small_irriducible(char, degree):
        """
        Return an irreducible polynomial of self.degree with a small
        number of non-zero coefficients.
        """
        cardinality = char ** degree
        basefield = FinitePrimeField.getInstance(char)
        top = uniutil.polynomial({degree: 1}, coeffring=basefield)
        for seed in range(degree - 1):
            for const in range(1, char):
                coeffs = [const] + arith1.expand(seed, 2)
                cand = uniutil.polynomial(enumerate(coeffs), coeffring=basefield) + top
                if cand.isirreducible():
                    _log.debug(cand.order.format(cand))
                    return cand
        for subdeg in range(degree):
            subseedbound = char ** subdeg
            for subseed in range(subseedbound + 1, char * subseedbound):
                if not subseed % char:
                    continue
                seed = subseed + cardinality
                cand = uniutil.polynomial(enumerate(arith1.expand(seed, cardinality)), coeffring=basefield)
                if cand.isirreducible():
                    return cand

    @staticmethod
    def _primitive_polynomial(char, degree):
        """
        Return a primitive polynomial of self.degree.

        REF: Lidl & Niederreiter, Introduction to finite fields and
             their applications.
        """
        cardinality = char ** degree
        basefield = FinitePrimeField.getInstance(char)
        const = basefield.primitive_element()
        if degree & 1:
            const = -const
        cand = uniutil.polynomial({0:const, degree:basefield.one}, basefield)
        maxorder = factor_misc.FactoredInteger((cardinality - 1) // (char - 1))
        var = uniutil.polynomial({1:basefield.one}, basefield)
        while not (cand.isirreducible() and
                   all(pow(var, int(maxorder) // p, cand).degree() > 0 for p in maxorder.prime_divisors())):
            # randomly modify the polynomial
            deg = bigrandom.randrange(1, degree)
            coeff = basefield.random_element(1, char)
            cand += uniutil.polynomial({deg:coeff}, basefield)
        _log.debug(cand.order.format(cand))
        return cand

    @staticmethod
    def _scalar_mul(integer, fqelem):
        """
        Return integer * (Fq element).
        """
        return fqelem.__class__(fqelem.rep * integer, fqelem.field)

    def card(self):
        """
        Return the cardinality of the field
        """
        return self.char ** self.degree

    def createElement(self, seed):
        """
        Create an element of the field.
        """
        if isinstance(seed, (int, long)):
            expansion = arith1.expand(seed, card(self.basefield))
            return ExtendedFieldElement(
                uniutil.FiniteFieldPolynomial(enumerate(expansion), self.basefield),
                self)
        elif seed in self.basefield:
            return ExtendedFieldElement(
                uniutil.FiniteFieldPolynomial([(0, seed)], self.basefield),
                self)
        elif seed in self:
            # seed is in self, return only embedding
            return self.zero + seed
        elif (isinstance(seed, uniutil.FiniteFieldPolynomial) and
              seed.getCoefficientRing() is self.basefield):
            return ExtendedFieldElement(seed, self)
        else:
            try:
                # lastly check sequence
                return ExtendedFieldElement(
                    uniutil.FiniteFieldPolynomial(enumerate(seed), self.basefield),
                    self)
            except TypeError:
                raise TypeError("seed %s is not an appropriate object." % str(seed))

    def __repr__(self):
        return "%s(%d, %d)" % (self.__class__.__name__, self.char, self.degree)

    def __str__(self):
        return "F_%d @(%s)" % (card(self), str(self.modulus))

    def __hash__(self):
        return (self.char ** self.degree) & 0xFFFFFFFF

    def issuperring(self, other):
        """
        Report whether the field is a superring of another ring.
        """
        if self is other:
            return True
        elif isinstance(other, ExtendedField):
            if self.char == other.char and not (self.degree % other.degree):
                return True
            return False
        elif isinstance(other, FinitePrimeField):
            if self.char == other.getCharacteristic():
                return True
            return False
        try:
            return other.issubring(self)
        except:
            return False

    def issubring(self, other):
        """
        Report whether the field is a subring of another ring.
        """
        if self is other:
            return True
        elif isinstance(other, FinitePrimeField):
            return False
        elif isinstance(other, ExtendedField):
            if self.char == other.char and not (other.degree % self.degree):
                return True
            return False
        try:
            return other.issuperring(self)
        except:
            return False

    def __contains__(self, elem):
        """
        Report whether elem is in field.
        """
        if (isinstance(elem, ExtendedFieldElement) and
            elem.getRing().modulus == self.modulus):
            return True
        elif (isinstance(elem, FinitePrimeFieldElement) and
              elem.getRing().getCharacteristic() == self.getCharacteristic()):
            return True
        return False

    def __eq__(self, other):
        """
        Equality test.
        """
        if isinstance(other, ExtendedField):
            return self.char == other.char and self.degree == other.degree
        return False

    def primitive_element(self):
        """
        Return a primitive element of the field, i.e., a generator of
        the multiplicative group.
        """
        fullorder = card(self) - 1
        if self._orderfactor is None:
            self._orderfactor = factor_misc.FactoredInteger(fullorder)
        for i in bigrange.range(self.char, card(self)):
            g = self.createElement(i)
            for p in self._orderfactor.prime_divisors():
                if g ** (fullorder // p) == self.one:
                    break
            else:
                return g

    # properties
    def _getOne(self):
        "getter for one"
        if self._one is None:
            self._one = ExtendedFieldElement(
                uniutil.FiniteFieldPolynomial([(0, 1)], self.basefield),
                self)
        return self._one

    one = property(_getOne, None, None, "multiplicative unit.")

    def _getZero(self):
        "getter for zero"
        if self._zero is None:
            self._zero = ExtendedFieldElement(
                uniutil.FiniteFieldPolynomial([], self.basefield),
                self)
        return self._zero

    zero = property(_getZero, None, None, "additive unit.")


def fqiso(f_q, gfq):
    """
    Return isomorphism function of extended finite fields from f_q to gfq.
    """
    if f_q is gfq:
        return lambda x: x
    if card(f_q) != card(gfq):
        raise TypeError("both fields must have the same cardinality.")

    # find a root of f_q's defining polynomial in gfq.
    p = f_q.getCharacteristic()
    q = card(f_q)
    for i in bigrange.range(p, q):
        root = gfq.createElement(i)
        if not f_q.modulus(root):
            break

    # finally, define a function
    def f_q_to_gfq_iso(f_q_elem):
        """
        Return the image of the isomorphism of the given element.
        """
        if not f_q_elem:
            return gfq.zero
        if f_q_elem.rep.degree() == 0:
            # F_p elements
            return gfq.createElement(f_q_elem.rep)
        return f_q_elem.rep(root)

    return f_q_to_gfq_iso


def embedding(f_q1, f_q2):
    """
    Return embedding homomorphism function from f_q1 to f_q2,
    where q1 = p ** k1, q2 = p ** k2 and k1 divides k2.
    """
    if card(f_q1) == card(f_q2):
        return fqiso(f_q1, f_q2)
    # search multiplicative generators of both fields and relate them.
    # 0. initialize basic variables
    q1, q2 = card(f_q1), card(f_q2)

    # 1. find a multiplicative generator of f_q2
    f_q2_gen = f_q2.primitive_element()
    f_q2_subgen = f_q2_gen ** ((q2 - 1) // (q1 - 1))

    # 2. find a root of defining polynomial of f_q1 in f_q2
    image_of_x_1 = _findroot(f_q1, f_q2, f_q2_subgen)

    # 3. finally, define a function
    def f_q1_to_f_q2_homo(f_q1_elem):
        """
        Return the image of the isomorphism of the given element.
        """
        if not f_q1_elem:
            return f_q2.zero
        if f_q1_elem.rep.degree() == 0:
            # F_p elements
            return f_q2.createElement(f_q1_elem.rep)
        return f_q1_elem.rep(image_of_x_1)

    return f_q1_to_f_q2_homo

def double_embeddings(f_q1, f_q2):
    """
    Return embedding homomorphism functions from f_q1 and f_q2
    to the composite field.
    """
    identity = lambda x: x
    if f_q1 is f_q2:
        return (identity, identity)
    p = f_q2.getCharacteristic()
    k1, k2 = f_q1.degree, f_q2.degree
    if not k2 % k1:
        return (embedding(f_q1, f_q2), identity)
    if not k1 % k2:
        return (identity, embedding(f_q2, f_q1))
    composite = FiniteExtendedField(p, gcd.lcm(k1, k2))
    return (embedding(f_q1, composite), embedding(f_q2, composite))


def _findroot(f_q1, f_q2, f_q2_subgen):
    """
    Find root of the defining polynomial of f_q1 in f_q2
    """
    if card(f_q1) > 50: # 50 is small, maybe
        _log.debug("by affine multiple (%d)" % card(f_q1))
        return affine_multiple_method(f_q1.modulus, f_q2)
    root = f_q2_subgen
    for i in range(1, card(f_q1)):
        if not f_q1.modulus(root):
            image_of_x_1 = root
            break
        root *= f_q2_subgen
    return image_of_x_1

def affine_multiple_method(lhs, field):
    """
    Find and return a root of the equation lhs = 0 by brute force
    search in the given field.  If there is no root in the field,
    ValueError is raised.

    The first argument lhs is a univariate polynomial with
    coefficients in a finite field.  The second argument field is
    an extension field of the field of coefficients of lhs.

    Affine multiple A(X) is $\sum_{i=0}^{n} a_i X^{q^i} - a$ for some
    a_i's and a in the coefficient field of lhs, which is a multiple
    of the lhs.
    """
    polynomial_ring = lhs.getRing()
    coeff_field = lhs.getCoefficientRing()
    q = card(coeff_field)
    n = lhs.degree()

    # residues = [x, x^q, x^{q^2}, ..., x^{q^{n-1}}]
    residues = [lhs.mod(polynomial_ring.one.term_mul((1, 1)))] # x
    for i in range(1, n):
        residues.append(pow(residues[-1], q, lhs)) # x^{q^i}

    # find a linear relation among residues and a constant
    coeff_matrix = matrix.createMatrix(n, n, [coeff_field.zero] * (n**2), coeff_field)
    for j, residue in enumerate(residues):
        for i in range(residue.degree() + 1):
            coeff_matrix[i + 1, j + 1] = residue[i]
    constant_components = [coeff_field.one] + [coeff_field.zero] * (n - 1)
    constant_vector = vector.Vector(constant_components)
    try:
        relation_vector, kernel = coeff_matrix.solve(constant_vector)
        for j in range(n, 0, -1):
            if relation_vector[j]:
                constant = relation_vector[j].inverse()
                relation = [constant * c for c in relation_vector]
                break
    except matrix.NoInverseImage:
        kernel_matrix = coeff_matrix.kernel()
        relation_vector = kernel_matrix[1]
        assert type(relation_vector) is vector.Vector
        for j in range(n, 0, -1):
            if relation_vector[j]:
                normalizer = relation_vector[j].inverse()
                relation = [normalizer * c for c in relation_vector]
                constant = coeff_field.zero
                break

    # define L(X) = A(X) + constant
    coeffs = {}
    for i, relation_i in enumerate(relation):
        coeffs[q**i] = relation_i
    linearized = uniutil.polynomial(coeffs, coeff_field)

    # Fq basis [1, X, ..., X^{s-1}]
    qbasis = [1]
    root = field.createElement(field.char)
    s = arith1.log(card(field), q)
    qbasis += [root**i for i in range(1, s)]
    # represent L as Matrix
    lmat = matrix.createMatrix(s, s, field.basefield)
    for j, base in enumerate(qbasis):
        imagei = linearized(base)
        if imagei.getRing() == field.basefield:
            lmat[1, j + 1] = imagei
        else:
            for i, coeff in imagei.rep.iterterms():
                if coeff:
                    lmat[i + 1, j + 1] = coeff
    # solve L(X) = the constant
    constant_components = [constant] + [coeff_field.zero] * (s - 1)
    constant_vector = vector.Vector(constant_components)
    solution, kernel = lmat.solve(constant_vector)
    assert lmat * solution == constant_vector
    solutions = [solution]
    for v in kernel:
        for i in range(card(field.basefield)):
            solutions.append(solution + i * v)

    # roots of A(X) contains the solutions of lhs = 0
    for t in bigrange.multirange([(card(field.basefield),)] * len(kernel)):
        affine_root_vector = solution
        for i, ti in enumerate(t):
            affine_root_vector += ti * kernel[i]
        affine_root = field.zero
        for i, ai in enumerate(affine_root_vector):
            affine_root += ai * qbasis[i]
        if not lhs(affine_root):
            return affine_root

    raise ValueError("no root found")


def FiniteExtendedField(characteristic, n_or_modulus):
    """
    Return ExtendedField F_{p^n} or F_p[]/(modulus).

    This is a convenience wrapper for backward compatibility.
    """
    return ExtendedField.prime_extension(characteristic, n_or_modulus)

FiniteExtendedFieldElement = ExtendedFieldElement
