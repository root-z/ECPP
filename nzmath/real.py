"""
real -- real numbers and the real number field.

The module real provides arbitrary precision real numbers and their
utilities.  The functions provided are corresponding to the math
standard module.
"""

from __future__ import division
import itertools
import warnings

import nzmath.arith1 as arith1
import nzmath.rational as rational
import nzmath.ring as ring
from nzmath.plugins import MATHMODULE as math, FLOATTYPE as Float, \
     CHECK_REAL_OR_COMPLEX as check_real_or_complex


class Real(ring.FieldElement):
    """
    Real is a class of real.
    This class is only for consistency for other Ring object.
    """

    convertable = (Float, int, long, rational.Rational)

    def __init__(self, value):
        """
        value will be wrapped in Float.
        """
        ring.FieldElement.__init__(self)
        if isinstance(value, rational.Rational):
            self.data = value.toFloat()
        else:
            self.data = Float(value)

    def __add__(self, other):
        if isinstance(other, Real):
            result = self.data + other.data
        elif isinstance(other, self.convertable):
            result = self.data + other
        else:
            return NotImplemented
        return self.__class__(result)

    def __radd__(self, other):
        if isinstance(other, self.convertable):
            result = other + self.data
        else:
            return NotImplemented
        return self.__class__(result)

    def __sub__(self, other):
        if isinstance(other, Real):
            result = self.data - other.data
        elif isinstance(other, self.convertable):
            result = self.data - other
        else:
            return NotImplemented
        return self.__class__(result)

    def __rsub__(self, other):
        if isinstance(other, self.convertable):
            result = other - self.data
        else:
            return NotImplemented
        return self.__class__(result)

    def __mul__(self, other):
        if isinstance(other, Real):
            result = self.data * other.data
        elif isinstance(other, self.convertable):
            result = self.data * other
        else:
            return NotImplemented
        return self.__class__(result)

    def __rmul__(self, other):
        if isinstance(other, self.convertable):
            result = other * self.data
        else:
            return NotImplemented
        return self.__class__(result)

    def __truediv__(self, other):
        if isinstance(other, Real):
            result = self.data / other.data
        elif isinstance(other, self.convertable):
            result = self.data / other
        else:
            return NotImplemented
        return self.__class__(result)

    __div__ = __truediv__

    def __rtruediv__(self, other):
        if isinstance(other, self.convertable):
            result = other / self.data
        else:
            return NotImplemented
        return self.__class__(result)

    __rdiv__ = __rtruediv__

    def __pow__(self, other):
        if isinstance(other, Real):
            result = math.pow(self.data, other.data)
        elif isinstance(other, self.convertable):
            result = math.pow(self.data, other)
        return result

    def __eq__(self, other):
        if isinstance(other, Real):
            return self.data == other.data
        elif isinstance(other, self.convertable):
            return self.data == other
        else:
            return NotImplemented

    def __hash__(self):
        return hash(self.data)

    def getRing(self):
        """
        Return the real field instance.
        """
        return theRealField


class RealField(ring.Field):
    """
    RealField is a class of the field of real numbers.
    The class has the single instance 'theRealField'.
    """

    def __init__(self):
        ring.Field.__init__(self)
        self._one = Real(1)
        self._zero = Real(0)

    def __str__(self):
        return "R"

    def __repr__(self):
        return "%s()" % (self.__class__.__name__, )

    def __contains__(self, element):
        if isinstance(element, (int, long, float, Float, Real, rational.Rational)):
            return True
        else:
            try:
                if check_real_or_complex(element):
                    return True
            except TypeError:
                if hasattr(element, 'conjugate'):
                    return element == element.conjugate()
                pass
        if hasattr(element, 'getRing') and element.getRing().issubring(self):
            return True
        return False

    def __eq__(self, other):
        return isinstance(other, self.__class__)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return 2

    # property one
    def _getOne(self):
        "getter for one"
        return self._one

    one = property(_getOne, None, None, "multiplicative unit.")

    # property zero
    def _getZero(self):
        "getter for zero"
        return self._zero

    zero = property(_getZero, None, None, "additive unit.")

    def issubring(self, aRing):
        if isinstance(aRing, RealField):
            return True
        elif self.issuperring(aRing):
            return False
        return aRing.issuperring(self)

    def issuperring(self, aRing):
        if isinstance(aRing, RealField):
            return True
        elif rational.theRationalField.issuperring(aRing):
            return True
        return aRing.issubring(self)

    def createElement(self, seed):
        return Float(seed)

    def getCharacteristic(self):
        """
        The characteristic of the real field is zero.
        """
        return 0


def log1piter(xx):
    " iterator for log(1+x)."
    d = 1
    positive = True
    t = rational.Rational(xx)
    yield t
    while True:
        d += 1
        positive = not positive
        t *= xx
        if positive:
            yield (t / d)
        else:
            yield (-t / d)

def floor(x):
    """
    floor(x) returns the integer; if x is an integer then x itself,
    otherwise the biggest integer less than x.
    """
    rx = rational.Rational(x)
    if rx.denominator == 1:
        return rational.Integer(rx.numerator)
    return rx.numerator // rx.denominator

def ceil(x):
    """
    ceil(x) returns the integer; if x is an integer then x itself,
    otherwise the smallest integer greater than x.
    """
    rx = rational.Rational(x)
    if rx.denominator == 1:
        return rational.Integer(rx.numerator)
    return rx.numerator // rx.denominator + 1

def tranc(x):
    """
    tranc(x) returns the integer; if x is an integer then x itself,
    otherwise the nearest integer to x.  If x has the fraction part
    1/2, then bigger one will be chosen.
    """
    rx = rational.Rational(x)
    if rx.denominator == 1:
        return rational.Integer(rx.numerator)
    return floor(x + rational.Rational(1, 2))

def fabs(x):
    """
    returns absolute value of x.
    """
    return abs(rational.Rational(x))

def fmod(x, y):
    """
    returns x - n * y, where n is the quotient of x / y, rounded
    towards zero to an integer.
    """
    fquot = rational.Rational(x) / y
    if fquot < 0:
        n = -floor(-fquot)
    else:
        n = floor(fquot)
    return x - n * y

def frexp(x):
    """
    Return a tuple (m, e) where x = m * 2 ** e, 1/2 <= abs(m) < 1 and
    e is an integer.
    This function is provided as the counter-part of math.frexp, but it
    might not be useful.
    """
    if x == 0:
        return (rational.Rational(0), 0)
    m = rational.Rational(x)
    e = 0
    if x > 0:
        while m >= 1:
            m /= 2
            e += 1
        while m < rational.Rational(1, 2):
            m *= 2
            e -= 1
    else:
        while m <= -1:
            m /= 2
            e += 1
        while m > rational.Rational(-1, 2):
            m *= 2
            e -= 1
    return (m, e)

def ldexp(x, i):
    """
    returns x * 2 ** i.
    """
    return x * 2 ** i

def EulerTransform(iterator):
    """
    Return an iterator which yields terms of Euler transform of the
    given iterator.
    """
    stock = []
    b = rational.Rational(1, 2)
    l = -1
    for term in iterator:
        stock.append(term)
        for i in xrange(l, -1, -1):
            stock[i] += stock[i+1]
        yield b * stock[0]
        b /= 2
        l += 1

# constants
theRealField = RealField()
