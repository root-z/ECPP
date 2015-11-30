"""
base classes for rings.
"""

from __future__ import division


class Ring (object):
    """
    Ring is an abstract class which expresses that
    the derived classes are (in mathematical meaning) rings.

    Definition of ring is as follows:
      Ring is a structure with addition and multiplication.  It is an
      abelian group with addition, and a monoid with multiplication.
      The multiplication obeys the distributive law.
    """

    def __init__(self):
        """
        Initialize _one and _zero for later use for properties 'one'
        and 'zero'.
        """
        # This class is abstract and cannot be instantiated.
        if type(self) is Ring:
            raise NotImplementedError("class Ring is abstract")
        self._one = None
        self._zero = None

    def createElement(self, seed):
        """
        createElement returns an element of the ring with seed.
        """
        raise NotImplementedError("derived class should override")

    def getCharacteristic(self):
        """
        Return the characteristic of the ring.

        The Characteristic of a ring is the smallest positive integer
        n s.t. n * a = 0 for any element a of the ring, or 0 if there
        is no such natural number.
        """
        raise NotImplementedError("derived class should override")

    def issubring(self, other):
        """
        Report whether another ring contains the ring as a subring.
        """
        raise NotImplementedError("derived class should override")

    def issuperring(self, other):
        """
        Report whether the ring is a superring of another ring.
        """
        raise NotImplementedError("derived class should override")

    def getCommonSuperring(self, other):
        """
        Return common super ring of self and another ring.
        """
        if self.issubring(other):
            return other
        elif self.issuperring(other):
            return self
        try:
            return other.getCommonSuperring(self)
        except RuntimeError:
            # reached recursion limit by calling on each other
            raise NotImplementedError("no common super ring")

    def __eq__(self, other):
        """
        Equality test.
        """
        raise NotImplementedError("derived class should override")

    def __hash__(self):
        raise NotImplementedError("derived class should override")

    def __ne__(self, other):
        """
        Inequality test.
        """
        return not self.__eq__(other)


class CommutativeRing (Ring):
    """
    CommutativeRing is an abstract subclass of Ring
    whose multiplication is commutative.
    """

    def __init__(self):
        """
        Initialize 'properties' attribute by an object of
        CommutativeRingProperties.
        """
        # This class is abstract and cannot be instantiated.
        if type(self) is CommutativeRing:
            raise NotImplementedError("class CommutativeRing is abstract")
        Ring.__init__(self)
        self.properties = CommutativeRingProperties()
        self._actions = {}
        self._actions_order = []

    def getQuotientField(self):
        """
        getQuotientField returns the quotient field of the ring
        if available, otherwise raises exception.
        """
        raise NotImplementedError

    def isdomain(self):
        """
        isdomain returns True if the ring is actually a domain,
        False if not, or None if uncertain.
        """
        return self.properties.isdomain()

    def isnoetherian(self):
        """
        isnoetherian returns True if the ring is actually a Noetherian
        domain, False if not, or None if uncertain.
        """
        return self.properties.isnoetherian()

    def isufd(self):
        """
        isufd returns True if the ring is actually a unique
        factorization domain, False if not, or None if uncertain.
        """
        return self.properties.isufd()

    def ispid(self):
        """
        ispid returns True if the ring is actually a principal
        ideal domain, False if not, or None if uncertain.
        """
        return self.properties.ispid()

    def iseuclidean(self):
        """
        iseuclidean returns True if the ring is actually a Euclidean
        domain, False if not, or None if uncertain.
        """
        return self.properties.iseuclidean()

    def isfield(self):
        """
        isfield returns True if the ring is actually a field,
        False if not, or None if uncertain.
        """
        return self.properties.isfield()

    def registerModuleAction(self, action_ring, action):
        """
        Register a ring 'action_ring', which act on the ring through
        'action' so the ring be an 'action_ring' module.
        """
        self._actions[action_ring] = action
        for i in range(len(self._actions_order) - 1, -1, -1):
            if self._actions_order[i].issubring(action_ring):
                self._actions_order.insert(i + 1, action_ring)
                break
        else:
            self._actions_order.insert(0, action_ring)

    def hasaction(self, action_ring):
        """
        Return True if 'action_ring' is registered to provide action.
        """
        if action_ring in self._actions:
            return True
        for action_superring in self._actions_order:
            if action_ring.issubring(action_superring):
                return True
        return False

    def getaction(self, action_ring):
        """
        Return the registered action for 'action_ring'.
        """
        for action_superring in self._actions_order:
            if action_ring.issubring(action_superring):
                return self._actions[action_superring]
        raise KeyError("no action is defined")


class Field (CommutativeRing):
    """
    Field is an abstract class which expresses that
    the derived classes are (in mathematical meaning) fields.
    """

    def __init__(self):
        """
        Set field flag True of 'properties' attribute.
        """
        # This class is abstract and cannot be instantiated.
        if type(self) is Field:
            raise NotImplementedError
        CommutativeRing.__init__(self)
        self.properties.setIsfield(True)

    def createElement(self, *args):
        """
        createElement returns an element of the field.
        """
        raise NotImplementedError

    def isfield(self):
        """
        Field overrides isfield of CommutativeRing.
        """
        return True

    def gcd(self, a, b):
        """
        A field is trivially a ufd and shuold provide gcd.
        """
        if not a and not b:
            return self.zero
        return self.one

    def getQuotientField(self):
        """
        getQuotientField returns the quotient field of the field.
        It is, of course, itself.
        """
        return self


class QuotientField (Field):
    """
    QuotientField is a class of quotient field.
    """

    def __init__(self, domain):
        """
        Create quotient field from given domain.
        Initialize 'basedomain' attribute by the given 'domain'.
        """
        # This class is abstract and cannot be instantiated.
        if type(self) is QuotientField:
            raise NotImplementedError("QuotientField is an abstract class")
        Field.__init__(self)
        self.basedomain = domain

        def baseaction(baseelement, quotient):
            return quotient.__class__(baseelement * quotient.numerator, quotient.denominator)
        
        self.registerModuleAction(self.basedomain, baseaction)


class RingElement (object):
    """
    RingElement is an abstract class for elements of rings.
    """

    def __init__(self, *args, **kwd):
        """
        This class is abstract and cannot be instantiated.
        """
        if type(self) is RingElement:
            raise NotImplementedError("RingElement is an abstract class.")

    def getRing(self):
        """
        getRing returns an object of a subclass of Ring,
        to which the element belongs.
        """
        raise NotImplementedError

    def __eq__(self, other):
        """
        Equality test.
        """
        raise NotImplementedError

    def __hash__(self):
        raise NotImplementedError

    def __ne__(self, other):
        """
        Inequality test.
        """
        return not (self == other)


class CommutativeRingElement (RingElement):
    """
    CommutativeRingElement is an abstract class for elements of
    commutative rings.
    """

    def __init__(self):
        """
        This class is abstract and cannot be instantiated.
        """
        if type(self) is CommutativeRingElement:
            raise NotImplementedError("CommutativeRingElement is an abstract class.")
        RingElement.__init__(self)

    def mul_module_action(self, other):
        """
        Return the result of a module action.
        other must be in one of the action rings of self's ring.
        """
        try:
            self_ring = self.getRing()
            other_ring = getRing(other)
            if self_ring.hasaction(other_ring):
                return self_ring.getaction(other_ring)(other, self)
        except RuntimeError, e:
            #print "mul_module_action", e
            raise
        raise TypeError("no module action with %s" % str(other_ring))

    def exact_division(self, other):
        """
        In UFD, if other divides self, return the quotient as a UFD
        element.  The main difference with / is that / may return the
        quotient as an element of quotient field.

        Simple cases:
        - in a euclidean domain, if remainder of euclidean division
          is zero, the division // is exact.
        - in a field, there's no difference with /.

        If other doesn't divide self, raise ValueError.
        """
        if self.getRing().iseuclidean():
            quotient, remainder = divmod(self, other)
            if not remainder:
                return quotient
            raise ValueError("division is not exact")
        elif self.getRing().isfield():
            return self / other
        raise NotImplementedError("exact division is not defined")


class FieldElement (CommutativeRingElement):
    """
    FieldElement is an abstract class for elements of fields.
    """
    def __init__(self):
        """
        This class is abstract and cannot be instantiated.
        """
        if type(self) is FieldElement:
            raise NotImplementedError("FieldElement is an abstract class.")
        CommutativeRingElement.__init__(self)

    def exact_division(self, other):
        """
        in a field, all divisions are exact (without remainder).
        """
        return self / other

class QuotientFieldElement (FieldElement):
    """
    QuotientFieldElement class is an abstract class to be used as a
    super class of concrete quotient field element classes.
    """

    def __init__(self, numerator, denominator):
        FieldElement.__init__(self)
        self.numerator = numerator
        if denominator == 0:
            raise ZeroDivisionError
        self.denominator = denominator

    def __add__(self, other):
        if hasattr(other, "numerator") and hasattr(other, "denominator"):
            numerator = self.numerator*other.denominator + self.denominator*other.numerator
            denominator = self.denominator*other.denominator
            return self.__class__(numerator, denominator)
        try:
            return self + self.getRing().one.mul_module_action(other)
        except TypeError:
            return NotImplemented

    __radd__ = __add__

    def __sub__(self, other):
        if hasattr(other, "numerator") and hasattr(other, "denominator"):
            numerator = self.numerator*other.denominator - self.denominator*other.numerator
            denominator = self.denominator*other.denominator
            return self.__class__(numerator, denominator)
        try:
            return self - self.getRing().one.mul_module_action(other)
        except TypeError:
            return NotImplemented

    def __rsub__(self, other):
        if hasattr(other, "numerator") and hasattr(other, "denominator"):
            numerator = self.denominator*other.numerator - self.numerator*other.denominator
            denominator = self.denominator*other.denominator
            return self.__class__(numerator, denominator)
        try:
            return self.getRing().one.mul_module_action(other) - self
        except TypeError:
            return NotImplemented

    def __neg__(self):
        return self.__class__(-self.numerator, self.denominator)

    def __mul__(self, other):
        if hasattr(other, "numerator") and hasattr(other, "denominator"):
            numerator = self.numerator * other.numerator
            denominator = self.denominator * other.denominator
            return self.__class__(numerator, denominator)
        try:
            return self.mul_module_action(other)
        except TypeError:
            return NotImplemented

    __rmul__ = __mul__

    def __pow__(self, index):
        return self.__class__(self.numerator ** index, self.denominator ** index)

    def __truediv__(self, other):
        if hasattr(other, "numerator") and hasattr(other, "denominator"):
            numerator = self.numerator * other.denominator
            denominator = self.denominator * other.numerator
            return self.__class__(numerator, denominator)
        try:
            return self * self.getRing().one.mul_module_action(other).inverse()
        except TypeError:
            return NotImplemented

    def __rtruediv__(self, other):
        if hasattr(other, "numerator") and hasattr(other, "denominator"):
            numerator =  other.numerator * self.denominator
            denominator = other.denominator * self.numerator
            return self.__class__(numerator, denominator)
        try:
            return self.getRing().one.mul_module_action(other) * self.inverse()
        except TypeError:
            return NotImplemented

    __div__ = __truediv__

    def inverse(self):
        return self.__class__(self.denominator, self.numerator)

    def __eq__(self, other):
        if hasattr(other, "numerator") and hasattr(other, "denominator"):
            return self.numerator*other.denominator == self.denominator*other.numerator
        try:
            return self == self.getRing().one.mul_module_action(other)
        except TypeError:
            return NotImplemented

    def __hash__(self):
        return hash(self.numerator * self.denominator)

class Ideal (object):
    """
    Ideal class is an abstract class to represent the finitely
    generated ideals.  Because the finitely-generatedness is not a
    restriction for Noetherian rings and in the most cases only
    Noetherian rings are used, it is general enough.
    """

    def __init__(self, generators, aring):
        """
        Ideal(generators, ring) creates an ideal of the ring genarated
        by the generators.  generators must be an element of the ring
        or a list of elements of the ring.
        """
        if type(self) is Ideal:
            raise NotImplementedError("class Ideal is abstract")
        self.ring = aring
        if not generators:
            self.generators = [self.ring.zero]
        elif generators in self.ring:
            self.generators = [generators]
        else:
            self.generators = list(generators)

    def __add__(self, other):
        """
        I + J <=> I.__add__(J)

        where I+J = {i+j | i in I and j in J}
        """
        assert self.ring is other.ring
        if self == other:
            return self
        return self.__class__(self.generators + other.generators, self.ring)

    def __mul__(self, other):
        """
        I * J <=> I.__mul__(J)

        where I*J = {sum of i*j | i in I and j in J}
        """
        assert self.ring is other.ring
        generators = []
        for i in self.generators:
            for j in other.generators:
                generators.append(i * j)
        return self.__class__(generators, self.ring)

    def __eq__(self, other):
        """
        I == J <=> I.__eq__(J)
        """
        assert type(self) is type(other)
        if self is other:
            return True
        if self.ring is not other.ring:
            return False
        return self.issubset(other) and self.issuperset(other)

    def __hash__(self):
        raise NotImplementedError

    def __ne__(self, other):
        """
        I != J <=> I.__ne__(J)
        """
        return not self.__eq__(other)

    def __contains__(self, element):
        """
        e in I  <=>  I.__contains__(e)

        for e in the ring, to which the ideal I belongs.
        """
        raise NotImplementedError("should be overridden")

    def issubset(self, other):
        """
        Report whether another ideal contains this ideal.
        """
        if self is other:
            return True
        raise NotImplementedError("should be overridden")

    def issuperset(self, other):
        """
        Report whether this ideal contains another ideal.
        """
        if self is other:
            return True
        raise NotImplementedError("should be overridden")

    def reduce(self, element):
        """
        Reduce an element with the ideal to simpler representative.
        """
        raise NotImplementedError("should be overridden")


class ResidueClassRing (CommutativeRing):
    """
    A residue class ring R/I,
    where R is a commutative ring and I is its ideal.
    """

    def __init__(self, ring, ideal):
        """
        ResidueClassRing(ring, ideal) creates a resudue class ring.
        The ring should be an instance of CommutativeRing, and ideal
        must be an instance of Ideal, whose ring attribute points the
        same ring with the given ring.
        """
        CommutativeRing.__init__(self)
        self.ring = ring
        self.ideal = ideal
        if self.ring.isnoetherian():
            self.properties.setIsnoetherian(True)

    def __contains__(self, element):
        if isinstance(element, ResidueClass) and element.ideal == self.ideal:
            return True
        return False

    def __eq__(self, other):
        if isinstance(other, ResidueClassRing):
            return self.ideal == other.ideal
        return False

    def __hash__(self):
        return hash(self.ideal)

    # properties
    def _getOne(self):
        "getter for one"
        if self._one is None:
            seed = self.ring.one
            if seed is not None:
                self._one = ResidueClass(seed, self.ideal)
        return self._one

    one = property(_getOne, None, None, "multiplicative unit.")

    def _getZero(self):
        "getter for zero"
        if self._zero is None:
            seed = self.ring.zero
            if seed is not None:
                self._zero = ResidueClass(seed, self.ideal)
        return self._zero

    zero = property(_getZero, None, None, "additive unit.")


class ResidueClass (CommutativeRingElement):
    """
    Element of residue class ring x+I, where I is the modulus ideal
    and x is a representative element.
    """

    def __init__(self, x, ideal):
        CommutativeRingElement.__init__(self)
        self.x = x
        self.ideal = ideal

    def __pos__(self):
        return self.__class__(self.ideal.reduce(self.x ), self.ideal)

    def __add__(self, other):
        assert self.ideal == other.ideal
        return self.__class__(self.ideal.reduce(self.x + other.x), self.ideal)

    def __sub__(self, other):
        assert self.ideal == other.ideal
        return self.__class__(self.ideal.reduce(self.x - other.x), self.ideal)

    def __mul__(self, other):
        assert self.ideal == other.ideal
        return self.__class__(self.ideal.reduce(self.x * other.x), self.ideal)

    def __eq__(self, other):
        if isinstance(other, ResidueClass):
            if self.ideal == other.ideal:
                return (self.x - other.x) in self.ideal
        return False

    def __hash__(self):
        raise NotImplementedError 

    def getRing(self):
        """
        Return a ResidueClassRing object.
        This overrides the method inherited from CommutativeRingElement.
        """
        return ResidueClassRing(self.ideal.ring, self.ideal)


class CommutativeRingProperties (object):
    """
    boolean properties of ring.

    Each property can have one of three values; True, False, or None.
    Of cource True is true and False is false, and None means that the
    property is not set neither directly nor indirectly.
    """

    def __init__(self):
        self._isfield = None
        self._iseuclidean = None
        self._ispid = None
        self._isufd = None
        self._isnoetherian = None
        self._isdomain = None

    def isfield(self):
        """
        Return True/False according to the field flag value being set,
        otherwise return None.
        """
        return self._isfield

    def setIsfield(self, value):
        """
        Set True/False value to the field flag.
        Propergation:
          True -> euclidean
        """
        self._isfield = bool(value)
        if self._isfield:
            self.setIseuclidean(True)

    def iseuclidean(self):
        """
        Return True/False according to the euclidean flag value being
        set, otherwise return None.
        """
        return self._iseuclidean

    def setIseuclidean(self, value):
        """
        Set True/False value to the euclidean flag.
        Propergation:
          True  -> pid
          False -> field
        """
        self._iseuclidean = bool(value)
        if self._iseuclidean:
            self.setIspid(True)
        else:
            self.setIsfield(False)

    def ispid(self):
        """
        Return True/False according to the pid flag value being set,
        otherwise return None.
        """
        return self._ispid

    def setIspid(self, value):
        """
        Set True/False value to the pid flag.
        Propergation:
          True  -> ufd, noetherian
          False -> euclidean
        """
        self._ispid = bool(value)
        if self._ispid:
            self.setIsufd(True)
            self.setIsnoetherian(True)
        else:
            self.setIseuclidean(False)

    def isufd(self):
        """
        Return True/False according to the ufd flag value being set,
        otherwise return None.
        """
        return self._isufd

    def setIsufd(self, value):
        """
        Set True/False value to the ufd flag.
        Propergation:
          True  -> domain
          False -> pid
        """
        self._isufd = bool(value)
        if self._isufd:
            self.setIsdomain(True)
        else:
            self.setIspid(False)

    def isnoetherian(self):
        """
        Return True/False according to the noetherian flag value being
        set, otherwise return None.
        """
        return self._isnoetherian

    def setIsnoetherian(self, value):
        """
        Set True/False value to the noetherian flag.
        Propergation:
          True  -> domain
          False -> pid
        """
        self._isnoetherian = bool(value)
        if self._isnoetherian:
            self.setIsdomain(True)
        else:
            self.setIspid(False)

    def isdomain(self):
        """
        Return True/False according to the domain flag value being
        set, otherwise return None.
        """
        return self._isdomain

    def setIsdomain(self, value):
        """
        Set True/False value to the domain flag.
        Propergation:
          False  -> ufd, noetherian
        """
        self._isdomain = bool(value)
        if not self._isdomain:
            self.setIsufd(False)
            self.setIsnoetherian(False)


def getRingInstance(obj):
    """
    Return a RingElement instance which eqauls 'obj'.

    Mainly for python built-in objects such as int or float.
    """
    if isinstance(obj, RingElement):
        return obj
    elif isinstance(obj, (int, long)):
        import nzmath.rational as rational
        return rational.Integer(obj)
    elif isinstance(obj, float):
        import nzmath.real as real
        return real.Real(obj)
    elif isinstance(obj, complex):
        import nzmath.imaginary as imaginary
        return imaginary.Complex(obj)
    return None

def getRing(obj):
    """
    Return a ring to which 'obj' belongs.

    Mainly for python built-in objects such as int or float.
    """
    try:
        # if obj has its getRing method, use it.
        return obj.getRing()
    except AttributeError:
        if isinstance(obj, (int, long)):
            import nzmath.rational as rational
            return rational.theIntegerRing
        if isinstance(obj, float):
            import nzmath.real as real
            return real.theRealField
        if isinstance(obj, complex):
            import nzmath.imaginary as imaginary
            return imaginary.theComplexField
    return None

def inverse(obj):
    """
    Return the inverse of 'obj'.  The inverse can be in the quotient
    field, if the 'obj' is an element of non-field domain.

    Mainly for python built-in objects such as int or float.
    """
    if hasattr(obj, "inverse"):
        return obj.inverse()
    # special cases
    if isinstance(obj, (int, long)):
        import nzmath.rational as rational
        return rational.Rational(1, obj)
    elif isinstance(obj, (float, complex)):
        return 1 / obj
    # final trial
    try:
        # if true division does not fail, it should return the inverse
        return getRing(obj).one / obj
    except TypeError:
        raise NotImplementedError("no inversion method found")

def exact_division(self, other):
    """
    Return the division of 'self' by 'obj' if the division is exact.

    Mainly for python built-in objects such as int or float.
    """
    return getRingInstance(self).exact_division(other)
