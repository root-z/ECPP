"""
default -- default Python types and modules
"""

import math
import cmath

MATHMODULE = math
CMATHMODULE = cmath
FLOATTYPE = float
COMPLEXTYPE = complex
PRECISION_CHANGEABLE = False

def CHECK_REAL_OR_COMPLEX(testee):
    """
    Return 1 if testee is a real number, 0 if complex.
    If testee is not a complex number, raise an exception.
    """
    try:
        comp = COMPLEXTYPE(testee)
        return 0 == comp.imag
    except TypeError:
        raise

def SETPRECISION(prec):
    """
    Set precision.

    prec is the number of bits.
    """
    raise NotImplementedError("no precision setting")
    
