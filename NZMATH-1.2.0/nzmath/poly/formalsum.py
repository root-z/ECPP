"""
formalsum --- formal sum, or linear combination.

The formalsum is mathematically a finite sum of terms.
The term consists of two parts: coefficient and base.
The coefficients must be in a ring.  The base is arbitrary.

Two formalsums can be added in the following way.
If there are terms with common base, they are fused into a
new term with the same base and coefficients added.

A coefficient can be looked up from the base.
If the specified base does not appear in the formalsum,
it is null.

abstract data type FS:
  FS: list(tuple(B, C)) -> FS

  +: FS x FS -> FS
     FS -> FS
  -: FS x FS -> FS
     FS -> FS
  *: C x FS -> FS
     FS x C -> FS

  []: FS x B -> C
  bases: FS -> list(B)
  coefficeints: FS -> list(C)
  terms: FS -> list(tuple(B, C))
  ==: FS x FS -> bool
  !=: FS x FS -> bool
  len: FS -> int
  repr: FS -> str
"""

import nzmath.compatibility


class FormalSumContainerInterface (object):
    """
    Interface of formal sum container.
    Do not instantiate.
    """
    def __iter__(self):
        """
        Return the iterator.
        It is an alias of iterterms.
        """
        return self.iterterms()

    def __getitem__(self, base):
        """
        Return the coefficient of specified base.
        If there is no term of the base, return 0.
        """
        raise NotImplementedError("method '__getitem__' must be overridden")

    def __contains__(self, base):
        """
        base in self

        membership test.
        """
        return base in self.bases()

    def __len__(self):
        """
        Return the number of data entries.
        """
        raise NotImplementedError("method '__len__' must be overridden")

    def __eq__(self, other):
        """
        self == other

        This implementaion is not optimal for more structured
        descendants.
        """
        if self is other:
            return True
        if not isinstance(other, FormalSumContainerInterface):
            return False
        self_bases = set(self.iterbases())
        other_bases = set(other.iterbases())
        if self_bases != other_bases:
            return False
        for base in self_bases:
            if self[base] != other[base]:
                return False
        return True

    def __hash__(self):
        val = sum([hash(self[base]) for base in set(self.iterbases())]) 
        return val

    def __ne__(self, other):
        """
        self != other
        """
        return not self.__eq__(other)

    def __nonzero__(self):
        """
        Return True, if self has some nonzero coefficients.
        False, otherwise.
        """
        for c in self.itercoefficients():
            if c:
                return True
        return False

    def __hash__(self):
        """
        hash(self)
        """
        try:
            return sum([hash(b)*hash(c) for (b, c) in self]) & 0x7fff
        except TypeError:
            raise TypeError("unhashable")

    def __add__(self, other):
        """
        self + other
        """
        sum_coeff = dict(self.iterterms())
        for base, coeff in other:
            if base in sum_coeff:
                sum_coeff[base] += coeff
            else:
                sum_coeff[base] = coeff
        return self.construct_with_default([(b, c) for (b, c) in sum_coeff.iteritems() if c])

    def __sub__(self, other):
        """
        self - other
        """
        sub_coeff = dict(self.iterterms())
        for base, coeff in other:
            if base in sub_coeff:
                sub_coeff[base] -= coeff
            else:
                sub_coeff[base] = -coeff
        return self.construct_with_default([(b, c) for (b, c) in sub_coeff.iteritems() if c])

    def __neg__(self):
        """
        -self
        """
        raise NotImplementedError("method '__neg__' should be overridden")

    def __pos__(self):
        """
        +self
        """
        raise NotImplementedError("method '__pos__' should be overridden")

    def __mul__(self, other):
        """
        self * other

        Only scalar multiplication is defined.
        """
        raise NotImplementedError("method '__mul__' should be overridden")

    def __rmul__(self, other):
        """
        other * self

        This method is invoked only when type of other does not
        support multiplication with FormalSumContainerInterface
        """
        raise NotImplementedError("method '__rmul__' should be overridden")

    def iterterms(self):
        """
        iterator for (degree, coefficient) pairs.
        """
        raise NotImplementedError("method 'iterterms' must be overridden")

    def itercoefficients(self):
        """
        iterator for coefficients.
        """
        raise NotImplementedError("method 'itercoefficients' must be overridden")

    def iterbases(self):
        """
        iterator for degrees.
        """
        raise NotImplementedError("method 'iterbases' must be overridden")

    def terms(self):
        """
        Return a list of terms as (base, coefficient) pairs.
        """
        return list(self.iterterms())

    def coefficients(self):
        """
        Return a list of all coefficients.
        """
        return list(self.itercoefficients())

    def bases(self):
        """
        Return a list of all bases.
        """
        return list(self.iterbases())

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

    def coefficients_map(self, func):
        """
        Create a new formal sum container by applying func to each
        coefficient.
        """
        return self.terms_map(lambda x, y: (x, func(y)))

    def bases_map(self, func):
        """
        Create a new formal sum container by applying func to each
        base.
        """
        return self.terms_map(lambda x, y: (func(x), y))

    def construct_with_default(self, maindata):
        """
        Create a new formal sum container of the same class with self,
        with given only the maindata and use copy of self's data if
        necessary.
        """
        # Do not extend the signature of this method.
        raise NotImplementedError("method 'construct_with_default' must be overridden")


class DictFormalSum (FormalSumContainerInterface):
    """
    formalsum implementation based on dict.
    """
    def __init__(self, args, defaultvalue=None):
        """
        DictFormalSum(args)

        args can be any dict initial values.
        It makes a mapping from bases to coefficients.

        The optional argument defaultvalue is the default value for
        __getitem__, i.e., if there is no term with the specified
        base, a look up attempt returns the defaultvalue.
        """
        self._data = dict(args)
        self._defaultvalue = defaultvalue

    def __mul__(self, other):
        """
        self * other

        Only scalar multiplication is defined.
        """
        return self.scalar_mul(other)

    def scalar_mul(self, scale):
        """
        scalar multiplication.

        self * scale
        """
        return self.coefficients_map(lambda c: c * scale)

    def __rmul__(self, other):
        """
        other * self

        This method is invoked only when type of other does not
        support multiplication with FormalSumContainerInterface.
        """
        return self.rscalar_mul(other)

    def rscalar_mul(self, scale):
        """
        r-scalar multiplication (r- means as of r-methods of python
        special methods, where self is the right operand.)

        scale * self
        """
        return self.coefficients_map(lambda c: scale * c)

    def __neg__(self):
        """
        -self
        """
        return self.construct_with_default([(b, -c) for (b, c) in self])

    def __pos__(self):
        """
        +self
        """
        return self.construct_with_default(self._data)

    def __eq__(self, other):
        """
        self == other
        """
        if self is other:
            return True
        if not isinstance(other, DictFormalSum):
            return FormalSumContainerInterface.__eq__(self, other)
        if self._data == other._data:
            return True
        return False

    def __hash__(self):
        val = sum([hash(ele) for ele in self._data]) 
        return val
   
    def __getitem__(self, base):
        """
        self[base]

        no KeyError will be raised. Insteads, it returns defaultvalue
        specified on time of initialization.
        """
        return self._data.get(base, self._defaultvalue)

    def iterterms(self):
        """
        iterator for (base, coefficient) pairs.
        """
        return self._data.iteritems()

    def itercoefficients(self):
        """
        iterator for coefficients.
        """
        return self._data.itervalues()

    def iterbases(self):
        """
        iterator for bases.
        """
        return self._data.iterkeys()

    def __len__(self):
        """
        len(self)

        Return the number of terms.
        """
        return len(self._data)

    def __repr__(self): # for debug
        return "%s(%s)" % (self.__class__.__name__, repr(self._data))

    def construct_with_default(self, maindata):
        """
        Create a new formal sum container of the same class with self,
        with given only the maindata and use copy of self's data if
        necessary.
        """
        return self.__class__(maindata, defaultvalue=self._defaultvalue)


class ListFormalSum (FormalSumContainerInterface):
    """
    formalsum implementation based on List.
    """
    def __init__(self, args, defaultvalue=None):
        """
        ListFormalSum(args)

        args can be sequence of tuples of (base, coefficient) pair.
        It makes a mapping from bases to coefficients.
        """
        self._data = list(args)
        self._defaultvalue = defaultvalue

    def __mul__(self, other):
        """
        self * other

        Only scalar multiplication is defined.
        """
        return self.scalar_mul(other)

    def scalar_mul(self, scale):
        """
        scalar multiplication.

        sum * scale
        """
        return self.coefficients_map(lambda c: c * scale)

    def __rmul__(self, other):
        """
        other * self

        This method is invoked only when type of other does not
        support multiplication with FormalSum.
        """
        return self.rscalar_mul(other)

    def rscalar_mul(self, scale):
        """
        scalar multiplication.

        scale * sum
        """
        return self.coefficients_map(lambda c: scale * c)

    def __neg__(self):
        """
        -self
        """
        return self.construct_with_default([(b, -c) for (b, c) in self])

    def __pos__(self):
        """
        +self
        """
        return self.construct_with_default(self._data)

    def __eq__(self, other):
        """
        self == other
        """
        if self is other:
            return True
        if not isinstance(other, ListFormalSum):
            return FormalSumContainerInterface.__eq__(self, other)
        if set(self._data) == set(other._data):
            return True
        return False

    def __hash__(self):
        val = sum([hash(ele) for ele in set(self._data)]) 
        return val

    def __getitem__(self, base):
        """
        self[base]

        no KeyError will be raised. Insteads, it returns defaultvalue
        specified on time of initialization.
        """
        for b, c in self:
            if b == base:
                return c
        return self._defaultvalue

    def iterterms(self):
        """
        iterator for (base, coefficient) pairs.
        """
        return iter(self._data)

    def itercoefficients(self):
        """
        iterator for coefficients.
        """
        for b, c in self._data:
            yield c

    def iterbases(self):
        """
        iterator for bases.
        """
        for b, c in self._data:
            yield b

    def __len__(self):
        """
        len(self)

        Return the number of terms.
        """
        return len(self._data)

    def __repr__(self): # for debug
        return "%s(%s)" % (self.__class__.__name__, repr(self._data))

    def construct_with_default(self, maindata):
        """
        Create a new formal sum container of the same class with self,
        with given only the maindata and use copy of self's data if
        necessary.
        """
        return self.__class__(maindata, defaultvalue=self._defaultvalue)
