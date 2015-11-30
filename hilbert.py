'''
Compute Hilbert Class Polynomials
'''

from nzmath.arith1 import floorsqrt

def solve(d):
    #initialize
    t = 1
    b = d % 2
    r = floorsqrt((-d)/3)
    h = 0
    red = {}
