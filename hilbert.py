'''
Compute Hilbert Class Polynomials
'''

from nzmath.arith1 import floorsqrt
import math


def solve(d):
    # initialize
    t = 1
    b = d % 2
    r = floorsqrt((-d)/3)
    h = 0
    red = {}

    # outer loop
    while (b <= r) :
        m = (b*b - d) / 4
        a = 0
        while a * a <= m:
            if m % a != 0:
                continue
            c = m/a
            if b > a:
                continue
            # optional polynomial setup
            tau = (-b + 1j * math.sqrt(-d)) / (2*a)

def delta(q):
    return q