"""
mpmath  -- mpmath plug-in
"""

import mpmath

# mpmath does not provide fmod/pow
def _fmod(x, y):
    """
    returns x - n * y, where n is the quotient of x / y, rounded
    towards zero to an integer.
    """
    fquot = mpmath.mpf(x) / y
    if fquot < 0:
        n = -mpmath.floor(-fquot)
    else:
        n = mpmath.floor(fquot)
    return x - n * y

mpmath.fmod = _fmod
mpmath.pow = mpmath.power

MATHMODULE = mpmath
CMATHMODULE = mpmath
FLOATTYPE = mpmath.mpf
COMPLEXTYPE = mpmath.mpc
PRECISION_CHANGEABLE = True

def CHECK_REAL_OR_COMPLEX(testee):
    """
    Return 1 if testee is a real number, 0 if complex.
    If testee is not a complex number, raise an exception.
    """
    try:
        COMPLEXTYPE(testee)
    except TypeError:
	# not a complex number
        raise
    try:
        FLOATTYPE(testee)
        return 1
    except TypeError:
        return 0

def SETPRECISION(prec):
    """
    Set precision.

    prec is the number of bits.
    """
    mpmath.mp.prec = prec
