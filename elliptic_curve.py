from nzmath.arith1 import modsqrt
from nzmath import arith1
import random
from jacobi import jacobi


class EllipticCurve(object):
    '''
    Elliptic Curve class. Allows for operations on points.
    Uses the short equation.
    y^2 = x^3 + a * x + b

    The curve is defined on Fp.

    A point on the curve is a pair (2-tuple) (x, y),
    with the special identity point a single value 0.
    '''
    def __init__(self, a, b, p):
        self.a = a
        self.b = b
        self.p = p

    '''
    def y(self, x):
        """
        Given x compute y such that (x, y) is on E(Fp)
        Args:
            x: x-coordinate

        Returns:
            y-coordinate
        """
        y_square = (x ** 3 + self.a * x + self.b) % self.p
        return modsqrt(y_square, self.p)
    '''

    def add(self, P, Q):
        """
        Elliptic Curve Addition
        Args:
            P: Point on curve
            Q: Point on curve

        Returns:
            R = P + Q
        """
        if P == 0:
            return Q
        if Q == 0:
            return P
        if P[0] == Q[0] and P[1] == -Q[1]:
            return 0
        if P == Q:
            l = (3 * P[0]**2 + self.a) / (2*P[1])
        else:
            l = (P[1] - Q[1]) / (P[0] - Q[0])

        x = l - P[0] - Q[0]
        y = l * (P[0] - x) - P[1]
        return x, y

    def sub(self, P, Q):
        """

        Args:
            P:
            Q:

        Returns:
            R = P - Q
        """
        if Q == 0:
            return P
        negative_Q = (Q[0], -Q[1])
        return self.add(P, negative_Q)

    def mul(self, k, P):
        """

        Args:
            k:
            P:

        Returns:

        """
        if k < 0:
            P = self.sub(0, P)

        l = arith1.expand(k, 2)
        Q = 0
        for j in range(len(l) -1, -1, -1):
            Q = self.add(Q, Q)
            if l[j] == 1:
                Q = self.add(Q, P)
        return Q

    def random_point(self):
        """

        Returns:
            A random point on the curve.
        """
        x = random.randrange(self.p)
        y_square = x**3 + self.a * x + self.b
        while jacobi(y, self.p) == -1:
            x = random.randrange(self.p)
            y_square = x**3 + self.a * x + self.b
        y = modsqrt(y_square, self.p)
        return x, y