"""
misc functions using factorization.
"""

import nzmath.gcd as gcd
import nzmath.arith1 as arith1
import nzmath.prime as prime
import nzmath.factor.methods as methods


def primePowerTest(n):
    """
    This program using Algo. 1.7.5 in Cohen's book judges whether
    n is of the form p**k with prime p or not.
    If it is True, then (p,k) will be returned,
    otherwise (n,0).
    """
    if n & 1:
        q = n
        while True:
            if not prime.primeq(q):
                a = 2
                while prime.spsp(n, a):
                    a += 1
                d = gcd.gcd(pow(a,q,q) - a, q)
                if d == 1 or d == q:
                    return (n, 0)
                q = d
            else:
                p = q
                break
    else:
        p = 2

    k, q = arith1.vp(n, p)
    if q == 1:
        return (p, k)
    else:
        return (n, 0)


class FactoredInteger(object):
    """
    Integers with factorization information.
    """
    def __init__(self, integer, factors=None):
        """
        FactoredInteger(integer [, factors])

        If factors is given, it is a dict of type {prime:exponent}
        and the product of prime**exponent is equal to the integer.
        Otherwise, factorization is carried out in initialization.
        """
        self.integer = long(integer)
        if factors is None:
            self.factors = dict(methods.factor(self.integer))
        else:
            self.factors = dict(factors)

    @classmethod
    def from_partial_factorization(cls, integer, partial):
        """
        Construct a new FactoredInteger object from partial
        factorization information given as dict of type
        {prime:exponent}.
        """
        partial_factor = 1
        for p, e in partial.iteritems():
            partial_factor *= p**e
        assert not integer % partial_factor, "wrong factorization"
        return cls(integer // partial_factor) * cls(partial_factor, partial)

    def __iter__(self):
        """
        Default iterator
        """
        return self.factors.iteritems()

    def __mul__(self, other):
        if isinstance(other, FactoredInteger):
            integer = self.integer * other.integer
            new_factors = self.factors.copy()
            for p in other.factors:
                new_factors[p] = new_factors.get(p, 0) + other.factors[p]
            return self.__class__(integer, new_factors)
        else:
            return self * FactoredInteger(other)

    __rmul__ = __mul__

    def __pow__(self, other):
        new_integer = self.integer**other
        new_factors = {}
        for p in self.factors:
            new_factors[p] = self.factors[p] * other
        return self.__class__(new_integer, new_factors)

    def __pos__(self):
        return self.copy()

    def __str__(self):
        return str(self.integer)

    def __eq__(self, other):
        return self.integer == int(other)

    def __hash__(self):
        return hash(self.integer)

    def __ne__(self, other):
        return self.integer != int(other)

    def __long__(self):
        return self.integer

    __int__ = __long__

    def copy(self):
        return self.__class__(self.integer, self.factors.copy())

    def __mod__(self, other):
        # maybe you want self.is_divisible_by(other)
        if int(other) in self.factors:
            return 0
        return self.integer % int(other)

    def is_divisible_by(self, other):
        """
        Return True if other divides self.
        """
        if int(other) in self.factors:
            # other is prime and divides
            return True
        return not self.integer % int(other)

    def exact_division(self, other):
        """
        Divide by a factor.
        """
        divisor = int(other)
        quotient = self.copy()
        if divisor in quotient.factors:
            if quotient.factors[divisor] == 1:
                del quotient.factors[divisor]
            else:
                quotient.factors[divisor] -= 1
        elif not isinstance(other, FactoredInteger):
            dividing = divisor
            for p, e in self.factors.iteritems():
                while not dividing % p:
                    dividing //= p
                    if quotient.factors[p] == 1:
                        del quotient.factors[p]
                        assert dividing % p, dividing
                    else:
                        quotient.factors[p] -= 1
                if dividing == 1:
                    break
            assert dividing == 1
        else:
            for p, e in other.factors.iteritems():
                assert p in quotient.factors and quotient.factors[p] >= e
                if quotient.factors[p] == e:
                    del quotient.factors[p]
                else:
                    quotient.factors[p] -= e
        quotient.integer //= divisor
        return quotient

    # maybe this is what you want, isn't it?
    __floordiv__ = exact_division

    def divisors(self):
        """
        Return all divisors.
        """
        l = [1]
        for p, e in self.factors.iteritems():
            for j in range(1, e + 1):
                l += [k*pow(p, j) for k in l if k % p]
        l.sort()
        return l

    def proper_divisors(self):
        """
        Return the proper divisors (divisors of n excluding 1 and n).
        """
        return self.divisors()[1:-1]

    def prime_divisors(self):
        """
        Return the list of primes that divides the number.
        """
        return self.factors.keys()

    def square_part(self, asfactored=False):
        """
        Return the largest integer whose square divides the number.

        If an optional argument asfactored is true, then the result is
        also a FactoredInteger object. (default is False)
        """
        result = FactoredInteger(1, {})
        for d, e in self.factors.iteritems():
            if e >= 2:
                result *= FactoredInteger(d ** (e >> 1), {d:e>>1})
        if asfactored:
            return result
        else:
            return result.integer

    def squarefree_part(self, asfactored=False):
        """
        Return the largest divisor of the number which is squarefree.

        If an optional argument asfactored is true, then the result is
        also a FactoredInteger object. (default is False)
        """
        result = FactoredInteger(1, {})
        for d, e in self.factors.iteritems():
            if e & 1:
                result *= FactoredInteger(d, {d:1})
        if asfactored:
            return result
        else:
            return result.integer


# for backward compatibility
allDivisors = lambda n: FactoredInteger(n).divisors()
primeDivisors = lambda n: FactoredInteger(n).prime_divisors()
squarePart = lambda n: FactoredInteger(n).square_part()
