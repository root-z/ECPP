"""
Miscellaneous arithmetic functions
"""

import itertools
import random
import bisect
import nzmath.gcd as gcd
from nzmath.plugins import MATHMODULE as math


def floorsqrt(a):
    """
    Return the floor of square root of the given integer.
    """
    if a < (1 << 59):
        return int(math.sqrt(a))
    else:
        # Newton method
        x = pow(10, (log(a, 10) >> 1) + 1) # compute initial value
        while True:
            x_new = (x + a//x) >> 1
            if x <= x_new:
                return x
            x = x_new

def floorpowerroot(n, k, return_power = False):
    """
    Return the floor of k-th power root of the given integer n.
    Use Newton method.
    """
    if k == 1:
        if return_power:
            return n, n
        else:
            return n
    elif k == 2:
        rslt = floorsqrt(n)
        if return_power:
            return rslt, rslt ** 2
        else:
            return rslt

    if n <= 0:
        if 0 == n:
            if return_power:
                return 0, 0
            else:
                return 0
        elif 0 == k & 1:
            raise ValueError("%d has no real %d-th root." % (n, k))
        else:
            sign = -1
            n = -n
    else:
        sign = 1

    # compute initial values
    exp, rem = divmod(log(n) + 1, k)
    if 0 != rem: # not ceiling
        exp += 1
    q = n >> (exp * (k - 1))
    x = 1 << exp

    # iteration by tangent approximation
    while True:
        x += (q - x) // k
        z = x ** (k - 1)
        q = n // z
        if x <= q:
            break

    if sign < 0:
        x = -x

    if return_power:
        return x, x * z
    else:
        return x

def powerDetection(n, largest_exp=False):
    """
    param positive integer n
    param boolean largest_exp
    return integer x, k s.t. n == x ** k
           (2 <= k if exist else x, k == n, 1)
           if largest_exp is true then return largest k
    """
    from nzmath.prime import generator_eratosthenes as generator_eratosthenes

    ge = generator_eratosthenes(log(n, 2))
    for exp in ge:
        power_root, power = floorpowerroot(n, exp, True)
        if power == n:
            if largest_exp:
                x, k = powerDetection(power_root, True)
                return x, k * exp
            else:
                return power_root, exp

    return n, 1

def legendre(a, m):
    """
    This function returns the Legendre symbol (a/m).
    If m is an odd composite then this is the Jacobi symbol.
    """
    a = a % m
    symbol = 1
    while a != 0:
        while a & 1 == 0:
            a >>= 1
            if m & 7 == 3 or m & 7 == 5:
                symbol = -symbol
        a, m = m, a
        if a & 3 == 3 and m & 3 == 3:
            symbol = -symbol
        a = a % m
    if m == 1:
        return symbol
    return 0

def modsqrt(n, p, e=1):
    """
    This function returns one of the square roots of n for mod p**e.
    p must be an odd prime.
    e must be a positive integer.
    If 1 < e then n must be relatively prime with p.
    """
    if 1 < e:
        x = modsqrt(n, p)
        if 0 == x:
            raise ValueError("if 1 < e then n must be relatively prime with p")
        ppower = p
        z = inverse(x << 1, p)
        for i in range(e - 1):
            x += (n - x ** 2) // ppower * z % p * ppower
            ppower *= p
        return x
    
    symbol = legendre(n, p)
    if symbol == 1:
        pmod8 = p & 7
        if pmod8 != 1:
            n %= p
            if pmod8 == 3 or pmod8 == 7:
                x = pow(n, (p >> 2) + 1, p)
            else: # p & 7==5
                x = pow(n, (p >> 3) + 1, p)
                c = pow(x, 2, p)
                if c != n:
                    x = (x * pow(2, p >> 2, p)) % p
        else: #p & 7==1
            d = 2
            while legendre(d, p) != -1:
                d = random.randrange(3, p)
            s, t = vp(p-1, 2)
            A = pow(n, t, p)
            D = pow(d, t, p)
            m = 0
            for i in range(1, s):
                if pow(A*(D**m), 1 << (s-1-i), p) == (p-1):
                    m += 1 << i
            x = (pow(n, (t+1) >> 1, p) * pow(D, m >> 1, p)) % p
        return x
    elif symbol == 0:
        return 0
    else:
        raise ValueError("There is no solution")

def expand(n, m):
    """
    This function returns m-adic expansion of n.
    n and m should satisfy 0 <= n, 2 <= m.
    """
    k = []
    while n >= m:
        k.append(n % m)
        n //= m
    k.append(n)
    return k

def inverse(x, n):
    """
    This function returns inverse of x for modulo n.
    """
    x = x % n
    y = gcd.extgcd(n, x)
    if y[2] == 1:
        if y[1] < 0:
            r = n + y[1]
            return r
        else:
            return y[1]
    raise ZeroDivisionError("There is no inverse for %d modulo %d." % (x, n))

def CRT(nlist):
    """
    This function is Chinese Remainder Theorem using Algorithm 2.1.7
    of C.Pomerance and R.Crandall's book.

    For example:
    >>> CRT([(1,2),(2,3),(3,5)])
    23
    """
    r = len(nlist)
    if r == 1 :
        return nlist [ 0 ] [ 0 ]

    product = []
    prodinv = []
    m = 1
    for i in range(1, r):
        m = m*nlist[i-1][1]
        c = inverse(m, nlist[i][1])
        product.append(m)
        prodinv.append(c)

    M = product[r-2]*nlist[r-1][1]
    n = nlist[0][0]
    for i in range(1, r):
        u = ((nlist[i][0]-n)*prodinv[i-1]) % nlist[i][1]
        n += u*product[i-1]
    return n % M

def AGM(a, b):
    """
    Arithmetic-Geometric Mean.
    """
    x = (a+b) * 0.5
    y = math.sqrt(a*b)
    while abs(x-y) > y*1e-15:
        x, y = (x+y) * 0.5, math.sqrt(x*y)
    return x

#def _BhaskaraBrouncker(n):
#    """
#
#    _BhaskaraBrouncker returns the minimum tuple (p,q) such that:
#        p ** 2 - n * q ** 2 = 1 or -1,
#    for positive integer n, which is not a square.
#
#    A good approximation for square root of n is given by the ratio
#    p/q; the error is at most 1/2*q**2.
#
#    """
#    floorOfSqrt = floorsqrt(n)
#    a = floorOfSqrt
#    b0, b1 = 0, floorOfSqrt
#    c0, c1 = 1, n - floorOfSqrt * floorOfSqrt
#    p0, p1 = 1, floorOfSqrt
#    q0, q1 = 0, 1
#    while c1 != 1:
#        a = (floorOfSqrt + b1)//c1
#        b0, b1 = b1, a * c1 - b1
#        c0, c1 = c1, c0 + a * (b0 - b1)
#        p0, p1 = p1, p0 + a * p1
#        q0, q1 = q1, q0 + a * q1
#    return (p1, q1)

def vp(n, p, k=0):
    """
    Return p-adic valuation and indivisible part of given integer.

    For example:
    >>> vp(100, 2)
    (2, 25)

    That means, 100 is 2 times divisible by 2, and the factor 25 of
    100 is indivisible by 2.

    The optional argument k will be added to the valuation.
    """
    q = p
    while not (n % q):
        k += 1
        q *= p

    return (k, n // (q // p))

class _Issquare:
    """
    A class for testing whether a number is square or not.
    The function issquare is an instance of the class, indeed.
    """
    q11 = [1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0]
    q63 = [1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1,
           0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0,
           1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]
    q64 = [1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,
           0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
           0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
    q65 = [1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
           0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0,
           0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1]

    def __call__(self, c):
        """
        Test whether a given number is a square number or not.  If
        the number is a square number, the function returns its square
        root.  Otherwise zero is returned.
        """
        t = c & 63
        if not self.q64[t]:
            return 0
        r = c % 45045  # 45045 = 63 * 65 * 11
        if not self.q63[r % 63]:
            return 0
        if not self.q65[r % 65]:
            return 0
        if not self.q11[r % 11]:
            return 0
        t = floorsqrt(c)
        if t * t == c:
            return t
        else:
            return 0

# test whether a given number is a square number or not.
issquare = _Issquare()
"""test whether a given number is a square number or not."""


# choosing 4096 to cover at least 1000 digits
POWERS_LIST_LENGTH = 4096
POWERS_OF_2 = [2**i for i in range(POWERS_LIST_LENGTH)]

def log(n, base=2):
    """
    Return the integer part of logarithm of the given natural number
    'n' to the 'base'.  The default value for 'base' is 2.
    """
    if n < base:
        return 0
    if base == 10:
        return _log10(n)
    if base == 2 and n < POWERS_OF_2[-1]:
        # a shortcut
        return bisect.bisect_right(POWERS_OF_2, n) - 1
    fit = base
    result = 1
    stock = [(result, fit)]
    while fit < n:
        next = fit ** 2
        if next <= n:
            fit = next
            result += result
            stock.append((result, fit))
        else:
            break
    else: # just fit
        return result
    stock.reverse()
    for index, power in stock:
        prefit = fit * power
        if prefit == n:
            result += index
            break
        elif prefit < n:
            fit = prefit
            result += index
    return result

def _log10(n):
    return len(str(n))-1

def product(iterable, init=None):
    """
    product(iterable) is a product of all elements in iterable.  If
    init is given, the multiplication starts with init instead of the
    first element in iterable. If the iterable is empty, then init or
    1 will be returned.

    If iterable is an iterator, it will be exhausted.
    """
    iterator = iter(iterable)
    if init is None:
        try:
            result = iterator.next()
        except StopIteration:
            return 1 # empty product
    else:
        result = init
    for item in iterator:
        result *= item
    return result
