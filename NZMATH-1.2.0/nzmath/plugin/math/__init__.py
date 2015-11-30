"""
math plug-in modules.

plug-in for float, complex data types and math/cmath functions.

Choose the implementation among provided plug-ins through config.

API
===

All math plug-ins provides the following constants: MATHMODULE,
CMATHMODULE, FLOATTYPE and COMPLEXTYPE.

MATHMODULE
----------

A module provides at least the same functions in 'math' standard
library of Python 2.5.  All functions accept FLOATTYPE (see below).

CMATHMODULE
-----------

A module provides at least the same functions in 'cmath' standard
library of Python 2.5.  All functions accept COMPLEXTYPE (see below).

FLOATTYPE
---------

A data type represents an approximation of real number.

COMPLEXTYPE
-----------

A data type represents an approximation of complex number.  Typically,
real and imaginary (imag) parts are represented with FLOATTYPE.

CHECK_REAL_OR_COMPLEX
---------------------

A function return 1 if testee is a real number, 0 if complex.  If
testee is not a complex number, raise an exception.
"""
