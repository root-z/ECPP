from __future__ import division
# standard modules
import itertools
# NZMATH modules
import nzmath.real as real
import nzmath.rational as rational
import nzmath.ring as ring

from nzmath.plugins import MATHMODULE as math, CMATHMODULE as cmath, \
     FLOATTYPE as Float, CHECK_REAL_OR_COMPLEX as check_real_or_complex


class Complex (ring.FieldElement):
    """
    Complex is a class for complex numbers.  Each instance has a coupled
    numbers; real and imaginary part of the number.
    """
    def __init__(self, re, im=None):
        if im:
            self.real = re
            self.imag = im
        elif isinstance(re, complex) or isinstance(re, Complex):
            self.real = re.real
            self.imag = re.imag
        else:
            self.real = re
            self.imag = 0

    def __add__(self, other):
        try:
            re = self.real + other.real
            im = self.imag + other.imag
        except AttributeError:
            if other in real.theRealField:
                re = self.real + other
                im = +self.imag
            else:
                return NotImplemented
        return self.__class__(re, im)

    __radd__ = __add__

    def __sub__(self, other):
        try:
            re = self.real - other.real
            im = self.imag - other.imag
        except AttributeError:
            if other in real.theRealField:
                re = self.real - other
                im = +self.imag
            else:
                return NotImplemented
        return self.__class__(re, im)

    def __rsub__(self, other):
        try:
            re = other.real - self.real
            im = other.imag - self.imag
        except AttributeError:
            if other in real.theRealField:
                re = other - self.real
                im = -self.imag
            else:
                return NotImplemented
        return self.__class__(re, im)

    def __mul__(self, other):
        try:
            re = self.real * other.real - self.imag * other.imag
            im = self.real * other.imag + self.imag * other.real
        except AttributeError:
            if other in real.theRealField:
                re = self.real * other
                im = self.imag * other
            else:
                return NotImplemented
        return self.__class__(re, im)

    __rmul__ = __mul__

    def __div__(self, other):
        try:
            denominator = other.real ** 2 + other.imag ** 2
            re = (self.real * other.real + self.imag * other.imag) / denominator
            im = (self.imag * other.real - self.real * other.imag) / denominator
        except AttributeError:
            if other in real.theRealField:
                re = self.real / other
                im = self.imag / other
            else:
                return NotImplemented
        return self.__class__(re, im)

    __truediv__ = __div__

    def __rdiv__(self, other):
        denominator = self.real ** 2 + self.imag ** 2
        try:
            re = (self.real * other.real + self.imag * other.imag) / denominator
            im = (self.real * other.imag - self.imag * other.real) / denominator
        except AttributeError:
            if other in real.theRealField:
                re = other * self.real / denominator
                im = -self.imag * other / denominator
            else:
                return NotImplemented
        return self.__class__(re, im)

    __rtruediv__ = __rdiv__

    def __pow__(self, other):
        if rational.isIntegerObject(other):
            if other == 0:
                return rational.Integer(1)
            elif other == 1:
                return +self
            elif other < 0:
                return (self**(-other)).inverse()
            elif other == 2:
                return self.__class__(self.real ** 2 - self.imag ** 2, 2 * self.real * self.imag)
            else:
                return rational.Integer(other).actMultiplicative(self)
        return exp(other * log(self))

    def __eq__(self, other):
        try:
            return self.real == other.real and self.imag == other.imag
        except AttributeError:
            if other in real.theRealField:
                return self.imag == 0 and self.real == other
            else:
                return NotImplemented

    def __hash__(self):
        return hash(self.real**2 + self.imag**2)

    def __ne__(self, other):
        try:
            return self.real != other.real or self.imag != other.imag
        except AttributeError:
            if other in real.theRealField:
                return self.imag != 0 or self.real != other
            else:
                return NotImplemented

    def __abs__(self):
        if self.imag == 0:
            return abs(self.real)
        if self.real == 0:
            return abs(self.imag)
        return math.hypot(self.real, self.imag)

    def __pos__(self):
        if self.imag == 0:
            return +self.real
        return self.__class__(+self.real, +self.imag)

    def __neg__(self):
        return self.__class__(-self.real, -self.imag)

    def __nonzero__(self):
        return bool(self.real or self.imag)

    def __repr__(self):
        return "Complex(" + repr(self.real) + ", " + repr(self.imag) + ")"

    def __str__(self):
        return str(self.real) + " + " + str(self.imag) + "j"

    def inverse(self):
        denominator = self.real ** 2 + self.imag ** 2
        re = self.real / denominator
        im = -self.imag / denominator
        return self.__class__(re, im)

    def conjugate(self):
        return self.__class__(self.real, -self.imag)

    def copy(self):
        return self.__class__(self.real, self.imag)

    ## comparisons are prohibited
    def __lt__(self, other):
        raise TypeError("cannot compare complex numbers using <, <=, >, >=")

    def __le__(self, other):
        raise TypeError("cannot compare complex numbers using <, <=, >, >=")

    def __gt__(self, other):
        raise TypeError("cannot compare complex numbers using <, <=, >, >=")

    def __ge__(self, other):
        raise TypeError("cannot compare complex numbers using <, <=, >, >=")

    def arg(self):
        x = self.real
        y = self.imag
        return math.atan2(y,x)

    def __complex__(self):
        return complex(float(self.real), float(self.imag))

    def getRing(self):
        """
        Return the complex field instance.
        """
        return theComplexField


class ComplexField (ring.Field):
    """
    ComplexField is a class of the field of real numbers.
    The class has the single instance 'theComplexField'.
    """

    def __init__(self):
        ring.Field.__init__(self)
        self._one = Complex(1)
        self._zero = Complex(0)

    def __str__(self):
        return "C"

    def __repr__(self):
        return "%s()" % self.__class__.__name__

    def __contains__(self, element):
        reduced = +element
        if reduced in real.theRealField:
            return True
        if isinstance(reduced, complex) or isinstance(reduced, Complex):
            return True
        return False  ## How to know an object be complex ?

    def __eq__(self, other):
        return isinstance(other, ComplexField)

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        return 3

    def createElement(self, seed):
        return Complex(seed)

    def _getOne(self):
        return self._one

    one = property(_getOne, None, None, "multiplicative unit.")

    def _getZero(self):
        return self._zero

    zero = property(_getZero, None, None, "additive unit.")

    def issubring(self, aRing):
        if isinstance(aRing, ComplexField):
            return True
        elif self.issuperring(aRing):
            return False
        return aRing.issuperring(self)

    def issuperring(self, aRing):
        if isinstance(aRing, ComplexField):
            return True
        elif real.theRealField.issuperring(aRing):
            return True
        return aRing.issubring(self)

    def getCharacteristic(self):
        """
        The characteristic of the real field is zero.
        """
        return 0


theComplexField = ComplexField()

j = Complex(0,1)
