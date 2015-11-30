"""
Elliptic Curve Method Factorization

It's using Montgomery's parametrization.
"""

from __future__ import division
import logging
import random
import nzmath.arith1 as _arith1
import nzmath.gcd as _gcd
import nzmath.prime as _prime

_log = logging.getLogger('nzmath.factor.ecm')


# curve types
S = SUYAMA = 1
B = BERNSTEIN = 2
A1 = ASUNCION1 = 5
A2 = ASUNCION2 = 6
A3 = ASUNCION3 = 8
A4 = ASUNCION4 = 9
A5 = ASUNCION5 = 10
N3 = NAKAMURA3 = 13

class Curve (object):
    """
    Elliptic curves for factorization.
    """

    _CURVE_TYPES = (S, B, A1, A2, A3, A4, A5, N3)
    
    def __init__(self, c):
        """
        Initialize a Curve object with Montgomery's parameter c.
        """
        self.c = c
        self.c2 = (c + 2)//4

    @classmethod
    def get_random_curve_with_point(cls, curve_type, n, bounds):
        """
        Return the curve with parameter C and a point Q on the curve,
        according to the curve_type, factorization target n and the
        bounds for stages.

        curve_type should be one of the module constants corresponding
        to parameters:
          S: Suyama's parameter selection strategy
          B: Bernstein's [2:1], [16,18,4,2]
          A1: Asuncion's [2:1], [4,14,1,1]
          A2: Asuncion's [2:1], [16,174,4,41]
          A3: Asuncion's [3:1], [9,48,1,2]
          A4: Asuncion's [3:1], [9,39,1,1]
          A5: Asuncion's [4:1], [16,84,1,1]
          N3: Nakamura's [2:1], [28,22,7,3]
          This is a class method.
        """
        bound = bounds.first
        if curve_type not in cls._CURVE_TYPES:
            raise ValueError("Input curve_type is wrong.")
        if curve_type == SUYAMA:
            t = n
            while _gcd.gcd(t, n) != 1:
                sigma = random.randrange(6, bound + 1)
                u, v = (sigma**2 - 5) % n, (4*sigma) % n
                t = 4*(u**3)*v
            d = _arith1.inverse(t, n)
            curve = cls((((u - v)**3 * (3*u + v)) * d - 2) % n)
            start_point = Point(pow(u, 3, n), pow(v, 3, n))
        elif curve_type == BERNSTEIN:
            d = random.randrange(1, bound + 1)
            start_point = Point(2, 1)
            curve = cls((4*d + 2) % n)
        elif curve_type == ASUNCION1:
            d = random.randrange(1, bound + 1)
            start_point = Point(2, 1)
            curve = cls((d + 1) % n)
        elif curve_type == ASUNCION2:
            d = random.randrange(1, bound + 1)
            start_point = Point(2, 1)
            curve = cls((4*d + 41) % n)
        elif curve_type == ASUNCION3:
            d = random.randrange(1, bound + 1)
            start_point = Point(3, 1)
            curve = cls((d + 2) % n)
        elif curve_type == ASUNCION4:
            d = random.randrange(1, bound + 1)
            start_point = Point(3, 1)
            curve = cls((d + 1) % n)
        elif curve_type == ASUNCION5:
            d = random.randrange(1, bound + 1)
            start_point = Point(4, 1)
            curve = cls((d + 1) % n)
        elif curve_type == NAKAMURA3:
            d = random.randrange(1, bound + 1)
            start_point = Point(2, 1)
            curve = cls((7*d + 3) % n)
        return curve, start_point


class Point (tuple):
    """
    x and z projective coordinates of elliptic curves.
    """
    def __new__(cls, x, z):
        """
        Create a new Point object.
        """
        return super(Point, cls).__new__(Point, (x, z))

    def _get_x(self):
        "getter for x"
        return self[0]

    def _get_z(self):
        "getter for z"
        return self[1]

    x = property(_get_x, None, None, "x")
    z = property(_get_z, None, None, "z")

POINT_AT_INFINITY = Point(0, 0)


class Bounds (object):
    """
    Bounds for ECM.

    public attributes:
      first: bound for first stage
      second: bound for second stage
    """
    def __init__(self, n):
        """
        Create a Bounds object from target number n.
        """
        if _arith1.log(n, 10) < 20:
            self.first = 1000
        elif _arith1.log(n, 10) < 25:
            self.first = 10000
        else:
            self.first = 100000
        self.second = self.first * 50

    def increment(self, scale=10):
        """
        Multiply both bounds by optional argument scale (default to 10).
        """
        self.first *= scale
        self.second *= scale

    def __getitem__(self, key):
        """
        Return first by key 1, and second key 2.
        Otherwise KeyError is raised.
        """
        if key == 1:
            return self.first
        if key == 2:
            return self.second
        raise KeyError("Key must be 1 or 2.")

    def __hash__(self):
        """
        Bounds objects are mutable and thus unhashable.
        """
        raise TypeError("Bounds objects are unhashable")

    def __str__(self):
        """
        (first, second)
        """
        return "(%d, %d)" % (self.first, self.second)


def ecm(n, curve_type=A1, incs=3, trials=20, **options):
    """
    Find a factor of n with Elliptic Curve Method.
    An unsuccessful factorization returns 1.

    There are a few optional arguments.

    By 'curve_type', the function choose a family of curves.
    Please use a module constant to specify the curve_type.

    The second optional argument 'incs' specifies the number of times
    for changing bounds.  The function repeats factorization trials
    several times changing curves with a fixed bounds.

    The last optional argument 'trials' can control how quickly move
    on to the next higher bounds.
    """
    # verbosity
    verbose = options.get('verbose', False)
    if verbose:
        _log.setLevel(logging.DEBUG)
        _log.debug("verbose")
    else:
        _log.setLevel(logging.NOTSET)

    # trivial checks
    if _prime.primeq(n):
        _log.info("%d is prime!" % n)
        return n
    if _gcd.gcd(n, 6) != 1: 
        _log.info("%d is not coprime to 6!" % n)
        if n % 2 == 0:
            return 2
        if n % 3 == 0:
            return 3

    # main loop
    bounds = Bounds(n)
    for inc in range(incs):
        _log.info("bounds increment %d times" % inc)
        _log.debug("Bounds B1, B2 = %s" % (bounds,))
        for trial in range(trials):
            _log.info("Trial #: %d" % (trial + 1))
            curve, point = Curve.get_random_curve_with_point(curve_type, n, bounds)
            _log.debug("Curve param: %d" % curve.c)
            g = stage1(n, bounds, curve, point)
            if 1 < g < n:
                return g
            g = stage2(n, bounds, curve, point)
            if 1 < g < n:
                return g
        bounds.increment()
    _log.info("ECM 2 step test failed!")

    if not verbose:
        _log.setLevel(logging.DEBUG)
    return 1

def stage1(n, bounds, C, Q):
    """
    ECM stage 1 for factoring n.
    The upper bound for primes to be tested is bounds.first.
    It uses curve C and starting point Q.
    """
    for p in _prime.generator_eratosthenes(bounds.first):
        q = p
        while q < bounds.first:
            Q = mul(Q, p, C, n)
            q *= p
    g = _gcd.gcd(Q.z, n)
    _log.debug("Stage 1: %d" % g)
    return g

def stage2(n, bounds, C, Q):
    """
    ECM stage 2 for factoring n.
    The upper bounds for primes to be tested are stored in bounds.
    It uses curve C and starting point Q.
    """
    d1 = bounds.first + random.randrange(1, 16)
    d2 = 0
    while _gcd.gcd(d1, d2) != 1:
        d2 = random.randrange(2, d1//5 + 1) # We want to keep d2 small
    for i in range(1, bounds.second//d1):
        if _gcd.gcd(i, d2) != 1:
            continue
        for j in range(1, d1//2):
            if _gcd.gcd(j, d1) == 1:
                Q = mul(Q, i*d1 + j*d2, C, n)
                if i*d1 - j*d2 > bounds.first:
                    Q = mul(Q, i*d1 - j*d2, C, n)
        g = _gcd.gcd(Q.z, n)
        if 1 < g < n:
            _log.debug("Stage 2: %d" % g)
            return g
    _log.debug("Stage 2: %d" % g)
    return 1

def mul(Q, x, C, n):
    """
    Return x*Q on the curve C mod n.

    m*Q and (m+1)*Q are being tracked in the main loop.
    """
    if x == 0:
        return POINT_AT_INFINITY
    if x == 1:
        return Q
    if x == 2:
        return double(Q, C, n)
    minor = Q
    major = double(Q, C, n)
    binary = _arith1.expand(x, 2)
    lastbit, binary = binary[0], binary[1:]
    while binary:
        if binary.pop() == 1:
            minor = add(major, minor, Q, n)
            major = double(major, C, n)
        else:
            major = add(minor, major, Q, n)
            minor = double(minor, C, n)
    if lastbit:
        return add(minor, major, Q, n)
    return double(minor, C, n)

# Montgomery's Arithmetic

def double(Q, C, n):
    """
    Return the doubled Q on the curve C mod n.
    """
    u = (Q.x + Q.z)**2
    v = (Q.x - Q.z)**2
    t = C.c2 * (u - v) + v
    return Point(u * v % n, (u - v) * t % n)

def add(P, Q, M, n):
    """
    Return the sum of P and Q mod n with a help of auxiliary point M.
    M must be P - Q. The curve is not explicitly necessary.
    """
    u = (P.x + P.z) * (Q.x - Q.z)
    v = (P.x - P.z) * (Q.x + Q.z)
    w = (u + v)**2
    t = (u - v)**2
    return Point(M.z * w % n, M.x * t % n)
