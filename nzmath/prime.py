"""
A module for generating primes and testing primality.
"""

import logging
import warnings
import nzmath.arith1 as arith1
import nzmath.gcd as gcd
import nzmath.bigrandom as bigrandom
import nzmath.bigrange as bigrange
import nzmath.poly.array as array_poly
from nzmath.config import GRH
from nzmath.plugins import MATHMODULE as math

_log = logging.getLogger('nzmath.prime')
_log.setLevel(logging.DEBUG)


PRIMES_LE_31 = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31)
PRIMONIAL_31 = arith1.product(PRIMES_LE_31)


def trialDivision(n, bound=0):
    """
    Trial division primality test for an odd natural number.
    Optional second argument is a search bound of primes.
    If the bound is given and less than the sqaure root of n
    and True is returned, it only means there is no prime factor
    less than the bound.
    """

    if bound:
        m = min(bound, arith1.floorsqrt(n))
    else:
        m = arith1.floorsqrt(n)
    for p in bigrange.range(3, m+1, 2):
        if not (n % p):
            return False
    return True


def spsp(n, base, s=None, t=None):
    """
    Strong Pseudo-Prime test.  Optional third and fourth argument
    s and t are the numbers such that n-1 = 2**s * t and t is odd.
    """
    if s is None or t is None:
        s, t = arith1.vp(n-1, 2)
    z = pow(base, t, n)
    if z != 1 and z != n-1:
        j = 0
        while j < s:
            j += 1
            z = pow(z, 2, n)
            if z == n-1:
                break
        else:
            return False
    return True


def miller(n):
    """
    Miller's primality test.

    This test is valid under GRH.
    """
    s, t = arith1.vp(n - 1, 2)
    # The O-constant 2 by E.Bach
    ## ln(n) = log_2(n)/log_2(e) = ln(2) * log_2(n)
    ## 2 * ln(2) ** 2 <= 0.960906027836404
    bound = min(n - 1, 960906027836404 * (arith1.log(n) + 1) ** 2 // 10 ** 15 + 1)
    _log.info("bound: %d" % bound)
    for b in range(2, bound):
        if not spsp(n, b, s, t):
            return False
    return True

def millerRabin(n, times=20):
    """
    Miller-Rabin pseudo-primality test.  Optional second argument
    times (default to 20) is the number of repetition.  The error
    probability is at most 4**(-times).
    """
    s, t = arith1.vp(n - 1, 2)
    for _ in range(times):
        b = bigrandom.randrange(2, n-1)
        if not spsp(n, b, s, t):
            return False
    return True


def bigprimeq(z):
    """
    Giving up rigorous proof of primality, return True for a probable
    prime.
    """
    if int(z) != z:
        raise ValueError("non-integer for primeq()")
    elif z <= 1:
        return False
    elif gcd.gcd(z, PRIMONIAL_31) > 1:
        return (z in PRIMES_LE_31)
    return millerRabin(z)


def Lucas_chain(n, f, g, x_0, x_1):
    """
    Given an integer n, two functions f and g, and initial value (x_0, x_1),
    compute (x_n, x_{n+1}), where the sequence {x_i} is defined as:
      x_{2i} = f(x_i)
      x_{2i+1} = g(x_i, x_{i+1})
    """
    binary = arith1.expand(n, 2)
    u = x_0
    v = x_1
    while binary:
        if 1 == binary.pop():
            u, v = g(u, v), f(v)
        else:
            u, v = f(u), g(u, v)
    return u, v


def _lucas_test_sequence(n, a, b):
    """
    Return x_0, x_1, x_m, x_{m+1} of Lucas sequence of parameter a, b,
    where m = (n - (a**2 - 4*b / n)) >> 1.
    """
    d = a**2 - 4*b
    if (d >= 0 and arith1.floorsqrt(d) ** 2 == d) \
    or not(gcd.coprime(n, 2*a*b*d)):
        raise ValueError("Choose another parameters.")

    x_0 = 2
    inv_b = arith1.inverse(b, n)
    x_1 = ((a ** 2) * inv_b - 2) % n

    # Chain functions
    def even_step(u):
        """
        'double' u.
        """
        return (u**2 - x_0) % n

    def odd_step(u, v):
        """
        'add' u and v.
        """
        return (u*v - x_1) % n

    m = (n - arith1.legendre(d, n)) >> 1
    x_m, x_mplus1 = Lucas_chain(m, even_step, odd_step, x_0, x_1)

    return x_0, x_1, x_m, x_mplus1


def lpsp(n, a, b):
    """
    Lucas test.
    Return True if n is a Lucas pseudoprime of parameters a, b,
    i.e. with respect to x**2-a*x+b.
    """
    x_0, x_1, x_m, x_mplus1 = _lucas_test_sequence(n, a, b)

    return (x_1 * x_m - x_0 * x_mplus1) % n == 0


def fpsp(n, a, b):
    """
    Frobenius test.
    Return True if n is a Frobenius pseudoprime of parameters a, b,
    i.e. with respect to x**2-a*x+b.
    """
    x_0, x_1, x_m, x_mplus1 = _lucas_test_sequence(n, a, b)

    if (x_1 * x_m - x_0 * x_mplus1) % n == 0:
        euler_pow = pow(b, (n-1) >> 1, n)
        return (euler_pow * x_m) % n == 2
    else:
        return False


def by_primitive_root(n, divisors):
    """
    Return True iff n is prime.

    The method proves the primality of n by existence of a primitive
    root.  It requires a sequence of all prime divisors of n - 1.

    Lehmer,D.H., Tests for primality by the converse of Fermat's
    theorem, Bull.Amer.Math.Soc, Vol.33, pp.327-340, 1927.
    """
    m_order = n - 1
    primes = tuple(p for p in divisors if p > 2)
    for g in bigrange.range(2, n):
        jacobi = arith1.legendre(g, n)
        if jacobi == 0:
            return False
        elif jacobi == 1:
            continue
        if pow(g, m_order, n) != 1:
            return False
        if all(pow(g, m_order // p, n) != 1 for p in primes):
            return True
    return True # only when n=2, flow reaches this line


def full_euler(n, divisors):
    """
    Return True iff n is prime.

    The method proves the primality of n by the equality
    phi(n) = n - 1, where phi denotes the Euler totient.
    It requires a sequence of all prime divisors of n - 1.

    Brillhart,J. & Selfridge,J.L., Some Factorizations of $2^n\pm 1$
    and Related Results, Math.Comp., Vol.21, 87-96, 1967
    """
    m_order = n - 1
    primes = set(divisors)
    for g in bigrange.range(2, n):
        if pow(g, m_order, n) != 1:
            return False
        if 2 in primes:
            jacobi = arith1.legendre(g, n)
            if jacobi == 0:
                return False
            elif jacobi == -1:
                primes.remove(2)
        satisfied = [p for p in primes if p > 2 and pow(g, m_order // p, n) != 1]
        if satisfied:
            primes.difference_update(satisfied)
        if not primes:
            return True
    return True # only when n=2, flow reaches this line


def prime(s):
    """
    prime(n) returns the n-th prime number.
    """
    if s != int(s):
        raise ValueError("non-integer for prime()")
    elif s <= 0:
        raise ValueError("non-positive-integer for prime()")
    i = 1
    for p in generator():
        if i == s:
            return p
        i += 1
    # The following line should not be reached:
    raise ValueError("Too big number %d for prime(i)." % s)


def generator():
    """
    Generate primes from 2 to infinity.
    """
    yield 2
    yield 3
    yield 5
    coprimeTo30 = (7, 11, 13, 17, 19, 23, 29, 31)
    times30 = 0
    while True:
        for i in coprimeTo30:
            if primeq(i + times30):
                yield i + times30
        times30 += 30


def generator_eratosthenes(n):
    """
    Generate primes up to n (inclusive) using Eratosthenes sieve.
    """
    if n < 2:
        raise StopIteration

    yield 2
    if n <= 2:
        return

    yield 3

    # make list for sieve
    sieve_len_max = (n+1) >> 1
    sieve = [True, False, True]
    sieve_len = 3
    k = 5
    i = 2
    while sieve_len < sieve_len_max:
        if sieve[i]:
            yield k
            sieve_len *= k
            if sieve_len_max < sieve_len:
                sieve_len //= k
                # adjust sieve list length
                sieve *= sieve_len_max // sieve_len
                sieve += sieve[:(sieve_len_max - len(sieve))]
                sieve_len = sieve_len_max
            else:
                sieve = sieve * k
            for j in range(i, sieve_len, k):
                sieve[j] = False
        k += 2
        i += 1

    # sieve
    limit = arith1.floorsqrt(n)
    while k <= limit:
        if sieve[i]:
            yield k
            j = (k ** 2 - 1) >> 1
            while j < sieve_len_max:
                sieve[j] = False
                j += k
        k += 2
        i += 1

    # output result
    limit = (n - 1) >> 1
    while i <= limit:
        if sieve[i]:
            yield 2 * i + 1
        i += 1


def nextPrime(n):
    """
    Return the smallest prime bigger than the given integer.
    """
    if n <= 1:
        return 2
    if n == 2:
        return 3
    n += (1 + (n & 1)) # make n be odd.
    while not primeq(n):
        n += 2
    return n


def randPrime(n):
    """
    Return a random n-digits prime
    """
    if n <= 0:
        raise ValueError("input number must be natural number")

    if n == 1:
        return bigrandom.map_choice(lambda i: (2, 3, 5, 7)[i], 4)

    p = bigrandom.randrange(10**(n-1)+1, 10**n, 2)
    while not primeq(p):
        p += 2
    if p < 10**n:
        return p

    # in very rare case or n is too small case,
    # search continues from the lower bound.
    p = 10**(n-1) + 1
    while not primeq(p):
        p += 2
    return p


def smallSpsp(n, s=None, t=None):
    """
    4 spsp tests are sufficient to determine whether an integer less
    than 10**12 is prime or not.  Optional third and fourth argument
    s and t are the numbers such that n - 1 = 2**s * t and t is odd.
    """
    if s is None or t is None:
        s, t = arith1.vp(n - 1, 2)
    for p in (2, 13, 23, 1662803):
        if not spsp(n, p, s, t):
            return False
    return True


def primeq(n):
    """
    A convenient function for primatilty test. It uses one of
    trialDivision, smallSpsp or apr depending on the size of n.
    """
    if int(n) != n:
        raise ValueError("non-integer for primeq()")
    if n <= 1:
        return False

    if gcd.gcd(n, PRIMONIAL_31) > 1:
        return (n in PRIMES_LE_31)
    if n < 2000000:
        return trialDivision(n)
    if not smallSpsp(n):
        return False
    if n < 10 ** 12:
        return True
    if not GRH:
        return apr(n)
    else:
        return miller(n)


def primonial(p):
    """
    Return 2*3*...*p for given prime p.
    """
    return arith1.product(generator_eratosthenes(p))


# defs for APR algorithm

def _isprime(n, pdivisors=None):
    """
    Return True iff n is prime.

    The check is done without APR.
    Assume that n is very small (less than 10**12) or
    prime factorization of n - 1 is known (prime divisors are passed to
    the optional argument pdivisors as a sequence).
    """
    if gcd.gcd(n, PRIMONIAL_31) > 1:
        return (n in PRIMES_LE_31)
    elif n < 10 ** 12:
        # 1369 == 37**2
        # 1662803 is the only prime base in smallSpsp which has not checked
        return n < 1369 or n == 1662803 or smallSpsp(n)
    else:
        return full_euler(n, pdivisors)

def properDivisors(n):
    """
    Return proper divisors of n (divisors of n excluding 1 and n).

    It is only useful for a product of small primes.
    One can use FactoredInteger.proper_divisors() as well.

    DEPRECATION: this function will be removed in version 1.3.
    """
    warnings.warn(DeprecationWarning("properDivisor is deprecated"))

    if n in (1, 2, 3, 5, 7, 11, 13, 17, 19, 23):
        return []
    else:
        l = [1]
        for p, e in _factor(n):
            for j in range(1, e + 1):
                l += [k*pow(p, j) for k in l if k % p != 0]
        l.remove(1)
        l.remove(n)
        l.sort()
        return l

def _factor(n, bound=0):
    """
    Trial division factorization for a natural number.
    Optional second argument is a search bound of primes.
    If the bound is given and less than the sqaure root of n,
    result is not proved to be a prime factorization.
    """
    factors = []
    if not (n & 1):
        v2, n = arith1.vp(n, 2)
        factors.append((2, v2))
    m = _calc_bound(n, bound)
    p = 3
    while p <= m:
        if not (n % p):
            v, n = arith1.vp(n, p)
            factors.append((p, v))
            m = _calc_bound(n, bound)
        p += 2
    if n > 1:
        factors.append((n, 1))
    return factors

def _calc_bound(n, bound=0):
    if bound:
        m = min((bound, arith1.floorsqrt(n)))
    else:
        m = arith1.floorsqrt(n)
    return m

def primitive_root(p):
    """
    Return a primitive root of p.
    """
    pd = FactoredInteger(p - 1).proper_divisors()
    for i in bigrange.range(2, p):
        for d in pd:
            if pow(i, (p - 1)//d, p) == 1:
                break
        else:
            return i


class Zeta(object):
    """
    Represent linear combinations of roots of unity.
    """
    def __init__(self, size, pos=None, val=1):
        self.size = size
        self.z = [0] * self.size
        if pos is not None:
            self.z[pos % self.size] = val

    def __add__(self, other):
        if self.size == other.size:
            m = self.size
            zr_a = Zeta(m)
            for i in range(m):
                zr_a.z[i] = self.z[i] + other.z[i]
            return zr_a
        else:
            m = gcd.lcm(self.size, other.size)
            return self.promote(m) + other.promote(m)

    def __mul__(self, other):
        if not isinstance(other, Zeta):
            zr_m = Zeta(self.size)
            zr_m.z = [x*other for x in self.z]
            return zr_m
        elif self.size == other.size:
            zr_m = Zeta(self.size)
            other = +other
            for k in range(other.size):
                if not other.z[k]:
                    continue
                elif other.z[k] == 1:
                    zr_m = zr_m + (self<<k)
                else:
                    zr_m = zr_m + (self<<k)*other.z[k]
            return zr_m
        else:
            m = gcd.lcm(self.size, other.size)
            return self.promote(m)*other.promote(m)

    __rmul__ = __mul__

    def __lshift__(self, offset):
        """
        The name is shift but the meaning of function is rotation.
        """
        new = Zeta(self.size)
        new.z = self.z[-offset:] + self.z[:-offset]
        return new

    def __pow__(self, e, mod=0):
        if mod:
            return self._pow_mod(e, mod)

        r = Zeta(self.size, 0)
        if e == 0:
            return r
        z = +self
        while True:
            if e & 1:
                r = z*r
                if e == 1:
                    return r
            e //= 2
            z = z._square()

    def _pow_mod(self, e, mod):
        r = Zeta(self.size, 0)
        if e == 0:
            return r
        z = self % mod
        while True:
            if e & 1:
                r = z * r % mod
                if e == 1:
                    return r
            e //= 2
            z = z._square_mod(mod)

    def _square(self):
        return self * self

    def _square_mod(self, mod):
        return self * self % mod

    def __pos__(self):
        m = self.size
        z_p = Zeta(m)
        if m & 1 == 0:
            mp = m >> 1
            for i in range(mp):
                if self.z[i] > self.z[i+mp]:
                    z_p.z[i] = self.z[i] - self.z[i+mp]
                else:
                    z_p.z[i+mp] = self.z[i+mp] - self.z[i]
        else:
            p = 3
            while m % p:
                p += 2
            mp = m // p
            for i in range(mp):
                mini = self.z[i]
                for j in range(mp + i, m, mp):
                    if mini > self.z[j]:
                        mini = self.z[j]
                for j in range(i, m, mp):
                    z_p.z[j] = self.z[j] - mini
        return z_p

    def __mod__(self, mod):
        """
        Return a new Zeta instance modulo 'mod'.

        All entries are reduced to the least absolute residues.
        """
        new = Zeta(self.size)
        half = (mod - 1) >> 1 # -1 to make 1 (mod 2) be 1, not -1.
        new.z = [(x + half) % mod - half for x in self.z]
        return new

    def __setitem__(self, position, value):
        self.z[position % self.size] = value

    def __getitem__(self, position):
        return self.z[position % self.size]

    def promote(self, size):
        if size == self.size:
            return +self
        new = Zeta(size)
        r = size // self.size
        for i in range(self.size):
            new.z[i*r] = self.z[i]
        return new

    def __len__(self):
        return self.size

    def __eq__(self, other):
        for i in range(self.size):
            if self.z[i] != other.z[i]:
                return False
        return True

    def __hash__(self):
        return sum([hash(self.z[i]) for i in range(1, self.size)])

    def weight(self):
        """
        Return the number of nonzero entries.
        """
        return len(filter(None, self.z))

    def mass(self):
        """
        Return the sum of all entries.
        """
        return sum(self.z)


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
        self.integer = int(integer)
        if factors is None:
            self.factors = dict(_factor(self.integer))
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
        try:
            return self.integer == other.integer
        except AttributeError:
            return self.integer == int(other)

    def __hash__(self):
        return hash(self.integer)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __le__(self, other):
        try:
            return self.integer <= other.integer
        except AttributeError:
            return self.integer <= int(other)

    def __lt__(self, other):
        try:
            return self.integer < other.integer
        except AttributeError:
            return self.integer < int(other)

    def __gt__(self, other):
        try:
            return self.integer > other.integer
        except AttributeError:
            return self.integer > int(other)

    def __ge__(self, other):
        try:
            return self.integer >= other.integer
        except AttributeError:
            return self.integer >= int(other)

    def __long__(self):
        return int(self.integer)

    __int__ = __long__

    def __mod__(self, other):
        # maybe you want self.is_divisible_by(other)
        try:
            if other.integer in self.factors:
                return 0
            return self.integer % other.integer
        except AttributeError:
            if int(other) in self.factors:
                return 0
            return self.integer % int(other)

    def __rfloordiv__(self, other):
        # assume other is an integer.
        return other // self.integer

    def copy(self):
        return self.__class__(self.integer, self.factors.copy())

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
        divs = [FactoredInteger(1)]
        for p, e in self.factors.iteritems():
            q = FactoredInteger(1)
            pcoprimes = list(divs)
            for j in range(1, e + 1):
                q *= FactoredInteger(p, {p:1})
                divs += [k * q for k in pcoprimes]
        return divs

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


class TestPrime(object):
    primes = PRIMES_LE_31
    primecache = set(primes)

    def __init__(self, t=12):
        if isinstance(t, (int, long)):
            self.t = FactoredInteger(t)
        else:
            assert isinstance(t, FactoredInteger)
            self.t = t
        powerof2 = self.t.factors[2] + 2
        self.et = FactoredInteger(2 ** powerof2, {2:powerof2})
        for d in self.t.divisors():
            # d is an instance of FactoredInteger
            p = d.integer + 1
            if p & 1 and (p in self.primecache or _isprime(p, d.factors)):
                self.et = self.et * FactoredInteger(p, {p:1})
                if p in self.t.factors:
                    e = self.t.factors[p]
                    self.et = self.et * FactoredInteger(p**e, {p:e})
                self.primecache.add(p)

    def next(self):
        eu = []
        for p in self.primes:
            if p in self.t.factors:
                eu.append((p - 1) * p**(self.t.factors[p] - 1))
            else:
                eu.append(p - 1)
                break
        p = self.primes[eu.index(min(eu))]
        return self.__class__(self.t * FactoredInteger(p, {p:1}))


class Status:
    """
    status collector for apr.
    """
    def __init__(self):
        self.d = {}

    def yet(self, key):
        """
        set key's status be 'yet'.
        """
        self.d[key] = 0

    def done(self, key):
        """
        set key's status be 'done'.
        """
        self.d[key] = 1

    def yet_keys(self):
        """
        Return keys whose stati are 'yet'.
        """
        return [k for k in self.d if not self.d[k]]

    def isDone(self, key):
        """
        Return whether key's status is 'done' or not.
        """
        return self.d[key]

    def subodd(self, p, q, n, J):
        """
        Return the sub result for odd key 'p'.
        If it is True, the status of 'p' is flipped to 'done'.
        """
        s = J.get(1, p, q)
        Jpq = J.get(1, p, q)
        m = s.size
        for x in range(2, m):
            if x % p == 0:
                continue
            sx = Zeta(m)
            i = x
            j = 1
            while i > 0:
                sx[j] = Jpq[i]
                i = (i + x) % m
                j += 1
            sx[0] = Jpq[0]
            sx = pow(sx, x, n)
            s = s*sx % n
        s = pow(s, n//m, n)
        r = n % m
        t = 1
        for x in range(1, m):
            if x % p == 0:
                continue
            c = (r*x) // m
            if c:
                tx = Zeta(m)
                i = x
                j = 1
                while i > 0:
                    tx[j] = Jpq[i]
                    i = (i + x) % m
                    j += 1
                tx[0] = Jpq[0]
                tx = pow(tx, c, n)
                t = t*tx % n
        s = +(t*s % n)
        if s.weight() == 1 and s.mass() == 1:
            for i in range(1, m):
                if gcd.gcd(m, s.z.index(1)) == 1:
                    self.done(p)
                return True
        return False

    def sub8(self, q, k, n, J):
        s = J.get(3, q)
        J3 = J.get(3, q)
        m = len(s)
        sx_z = {1:s}
        x = 3
        step = 2
        while m > x:
            z_4b = Zeta(m)
            i = x
            j = 1
            while i != 0:
                z_4b[j] = J3[i]
                i = (i + x) % m
                j += 1
            z_4b[0] = J3[0]
            sx_z[x] = z_4b
            s = pow(sx_z[x], x, n) * s
            step = 8 - step
            x += step

        s = pow(s, n//m, n)

        r = n % m
        step = 2
        x = 3
        while m > x:
            c = r*x
            if c > m:
                s = pow(sx_z[x], c//m, n) * s
            step = 8 - step
            x += step
        r = r & 7
        if r == 5 or r == 7:
            s = J.get(2, q).promote(m) * s
        s = +(s % n)

        if s.weight() == 1 and s.mass() == 1:
            if gcd.gcd(m, s.z.index(1)) == 1 and pow(q, (n-1) >> 1, n) == n-1:
                self.done(2)
            return True
        elif s.weight() == 1 and s.mass() == n-1:
            if gcd.gcd(m, s.z.index(n-1)) == 1 and pow(q, (n-1) >> 1, n) == n-1:
                self.done(2)
            return True
        return False

    def sub4(self, q, n, J):
        j2 = J.get(1, 2, q)**2
        s = q*j2 % n
        s = pow(s, n >> 2, n)
        if n & 3 == 3:
            s = s*j2 % n
        s = +(s % n)
        if s.weight() == 1 and s.mass() == 1:
            i = s.z.index(1)
            if (i == 1 or i == 3) and pow(q, (n-1) >> 1, n) == n-1:
                self.done(2)
            return True
        return False

    def sub2(self, q, n):
        s = pow(n-q, (n-1) >> 1, n)
        if s == n-1:
            if n & 3 == 1:
                self.done(2)
        elif s != 1:
            return False
        return True

    def subrest(self, p, n, et, J, ub=200):
        if p == 2:
            q = 5
            while q < 2*ub + 5:
                q += 2
                if not _isprime(q) or et % q == 0:
                    continue
                if n % q == 0:
                    _log.info("%s divides %s.\n" % (q, n))
                    return False
                k = arith1.vp(q-1, 2)[0]
                if k == 1:
                    if n & 3 == 1 and not self.sub2(q, n):
                        return False
                elif k == 2:
                    if not self.sub4(q, n, J):
                        return False
                else:
                    if not self.sub8(q, k, n, J):
                        return False
                if self.isDone(p):
                    return True
            else:
                raise ImplementLimit("limit")
        else:
            step = p*2
            q = 1
            while q < step*ub + 1:
                q += step
                if not _isprime(q) or et % q == 0:
                    continue
                if n % q == 0:
                    _log.info("%s divides %s.\n" % (q, n))
                    return False
                if not self.subodd(p, q, n, J):
                    return False
                if self.isDone(p):
                    return True
            else:
                raise ImplementLimit("limit")


class JacobiSum:
    """
    A repository of the Jacobi sums.
    """
    def __init__(self):
        self.shelve = {}

    def get(self, group, p, q=None):
        if q:
            assert group == 1
            if (group, p, q) not in self.shelve:
                self.make(q)
            return self.shelve[group, p, q]
        else:
            assert group == 2 or group == 3
            if (group, p) not in self.shelve:
                self.make(p)
            return self.shelve[group, p]

    def make(self, q):
        fx = self.makefx(q)
        qpred = q - 1
        qt = _factor(qpred)
        qt2 = [k for (p, k) in qt if p == 2][0]
        k, pk = qt2, 2**qt2
        if k >= 2:
            J2q = Zeta(pk, 1 + fx[1])
            for j in range(2, qpred):
                J2q[j + fx[j]] = J2q[j + fx[j]] + 1
            self.shelve[1, 2, q] = +J2q
            if k >= 3:
                J2 = Zeta(8, 3 + fx[1])
                J3 = Zeta(pk, 2 + fx[1])
                for j in range(2, qpred):
                    J2[j*3 + fx[j]] = J2[j*3 + fx[j]] + 1
                    J3[j*2 + fx[j]] = J3[j*2 + fx[j]] + 1
                self.shelve[3, q] = +(self.shelve[1, 2, q]*J3)
                self.shelve[2, q] = +(J2*J2)
        else:
            self.shelve[1, 2, q] = 1
        for (p, k) in qt:
            if p == 2:
                continue
            pk = p**k
            Jpq = Zeta(pk, 1 + fx[1])
            for j in range(2, qpred):
                Jpq[j + fx[j]] = Jpq[j + fx[j]] + 1
            self.shelve[1, p, q] = +Jpq
        del fx

    @staticmethod
    def makefx(q):
        """
        Return a dict called 'fx'.
        The dict stores the information that fx[i] == j iff
          g**i + g**j = 1 mod q
        for g, a primitive root of the prime q.
        """
        g = primitive_root(q)
        qpred = q - 1
        qd2 = qpred >> 1
        g_mf = [0, g]
        for _ in range(2, qpred):
            g_mf.append((g_mf[-1]*g) % q)
        fx = {}
        for i in range(1, qd2):
            if i in fx:
                continue
            # search j s.t. g**j + g**i = 1 mod q
            j = g_mf.index(q + 1 - g_mf[i])
            fx[i] = j
            fx[j] = i
            fx[qpred - i] = (j - i + qd2) % qpred
            fx[fx[qpred - i]] = qpred - i
            fx[qpred - j] = (i - j + qd2) % qpred
            fx[fx[qpred - j]] = qpred - j
        del g_mf
        return fx


def apr(n):
    """
    apr is the main function for Adleman-Pomerance-Rumery primality test.
    Assuming n has no prime factors less than 32.
    Assuming n is spsp for several bases.
    """
    L = Status()

    rb = arith1.floorsqrt(n) + 1
    el = TestPrime()
    while el.et <= rb:
        el = el.next()

    plist = el.t.factors.keys()
    plist.remove(2)
    L.yet(2)
    for p in plist:
        if pow(n, p-1, p*p) != 1:
            L.done(p)
        else:
            L.yet(p)
    qlist = el.et.factors.keys()
    qlist.remove(2)
    J = JacobiSum()
    for q in qlist:
        for p in plist:
            if (q-1) % p != 0:
                continue
            if not L.subodd(p, q, n, J):
                return False
        k = arith1.vp(q-1, 2)[0]
        if k == 1:
            if not L.sub2(q, n):
                return False
        elif k == 2:
            if not L.sub4(q, n, J):
                return False
        else:
            if not L.sub8(q, k, n, J):
                return False
    for p in L.yet_keys():
        if not L.subrest(p, n, el.et, J):
            return False
    r = int(n)
    for _ in bigrange.range(1, el.t.integer):
        r = (r*n) % el.et.integer
        if n % r == 0 and r != 1 and r != n:
            _log.info("%s divides %s.\n" %(r, n))
            return False
    return True


class ImplementLimit (Exception):
    """
    Exception throwed when the execution hits the implementation limit.
    """
    pass


# AKS
def aks(n):
    """
    AKS ( Cyclotomic Congruence Test ) primality test for a natural number.
    
    Input: a natural number n ( n > 1 ).
    Output: n is prime => return True
            n is not prime => return False.
    """
    import nzmath.multiplicative as multiplicative

    if arith1.powerDetection(n)[1] != 1: #Power Detection
        return False

    lg = math.log(n, 2)
    k = int(lg * lg)

    start = 3
    while 1:
        d = gcd.gcd(start, n)
        if 1 < d < n:
            return False
        x = n % start
        N = x
        for i in range(1, k + 1):
            if N == 1:
                break
            N = (N * x) % start
        if i == k:
            r = start
            break
        start += 1
    d = gcd.gcd(r, n)
    if 1 < d < n:
        return False
    if n <= r:
        return True

    e = multiplicative.euler(r) #Cyclotomic Conguence
    e = math.sqrt(e)
    e = int(e * lg)
    for b in range(1, e + 1):
        f = array_poly.ArrayPolyMod([b, 1], n)
        total = array_poly.ArrayPolyMod([1], n)
        count = n
        while count > 0:
            if count & 1:
                total = total * f
                total = _aks_mod(total, r)
            f = f.power()
            f = _aks_mod(f, r)
            count = count >> 1
        total_poly = total.coefficients_to_dict()
        if total_poly != {0:b, n % r:1}:
            return False
    return True

def _aks_mod(polynomial, r):
    """
    This function is used in aks.
    polynomial modulo (x^r - 1)
    """
    total = polynomial.coefficients[:r]
    aks_mod = polynomial.coefficients[r:]
    while len(aks_mod) - 1 >= r:
        for i in range(r):
            total[i] += aks_mod[i]
        aks_mod = aks_mod[r:]
    for i in range(len(aks_mod)):
        total[i] += aks_mod[i]
    return array_poly.ArrayPolyMod(total, polynomial.mod)
