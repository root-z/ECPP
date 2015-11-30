"""
plugins -- plug-in mechanism

Some function, such as floating point arithmetic, could have several
choices of modules. The 'plugins' module provides the mechanism to
choose one among these choices.  One can *plug-in* one module or
another for the function.

Usage:
The choice among plug-in modules is made through nzmath.config.
See the plug-in section of it for detail.

One can import all plugged-in constants as:
  from nzmath.plugins import *
but if he/she would like to select among features, it is possible:
  from nzmath.plugins import MATHMODULE as math
or something.

For Developers:
You should provide wrapper modules for plug-ins to ensure they have
the common interface.  Those plug-in wrapper modules should go into
the 'plugin' sub-package.
"""
from nzmath.config import PLUGIN_MATH

# math plug-ins
MATH_PLUGIN_CHOICE = ('mpmath', None)

if PLUGIN_MATH == MATH_PLUGIN_CHOICE[0]:
    from nzmath.plugin.math._mpmath import *
else:
    from nzmath.plugin.math.default import *


_MATH = ['MATHMODULE', 'CMATHMODULE', 'FLOATTYPE', 'COMPLEXTYPE',
         'CHECK_REAL_OR_COMPLEX',
         'PRECISION_CHANGEABLE', 'SETPRECISION']

__all__ = _MATH
