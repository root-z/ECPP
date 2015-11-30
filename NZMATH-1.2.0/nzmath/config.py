"""
config --- NZMATH configurations
"""

import os
import sys
import warnings

WINDOWS_PLATFORMS = ('win32', 'win64', 'cli', )

def nzmath_info():
	print("NZMATH [Version 1.2.0]")

# ----------------
# Default Settings
# ----------------
#
# Dependencies
# ============
#
# Some third party / platform dependent modules are possibly used,
# and they are configurable.
#
# mpmath
# ------
#
# mpmath is a package providing multiprecision math.
# http://code.google.com/p/mpmath
# This package is used in ecpp module.
#
# If you have mpmath installed, set as the following:
#HAVE_MPMATH = True
#CHECK_MPMATH = False
# Or, if you don't have mpmath installed, set as the following:
#HAVE_MPMATH = False
#CHECK_MPMATH = False
# The default values mean "I don't know; check it later":
HAVE_MPMATH = None
CHECK_MPMATH = True

#
# sqlite3
# -------
#
# sqlite3 is the default database module for Python, but it need to be
# enabled at the build time.
#
# If you have sqlite3 available, set as the following:
#HAVE_SQLITE3 = True
#CHECK_SQLITE3 = False
# Or, if you don't have sqlite3, set as the following:
#HAVE_SQLITE3 = False
#CHECK_SQLITE3 = False
# The default values mean "I don't know; check it later":
HAVE_SQLITE3 = None
CHECK_SQLITE3 = True


#
# net availability
# ----------------
#
# Some functions will connect to the Net.
# Desktop machines are usually connected to the Net, but notebooks
# may have connectivity only occasionally.
#
# If sometimes you have net connectivity, set as the following:
HAVE_NET = True
CHECK_NET = False
# Or, if your machine is not connected, set as the following:
#HAVE_NET = False
#CHECK_NET = False
# Or you can let NZMATH check connectivity every time, set as the following:
#HAVE_NET = None
#CHECK_NET = True

#
# psyco
# -----
#
# psyco is a package providing JIT compiler for Python.
# http://psyco.sourceforge.net/
# The package is usable only on x86 CPUs.
# This package is used in factor/mpqs module.
#
# If you have psyco installed, and would like to use it fully, set:
#HAVE_PSYCO = True
#CHECK_PSYCO = False
# If you would like to disable it no matter whether it is available:
#HAVE_PSYCO = False
#CHECK_PSYCO = False
# The default values mean "I don't know, check it and use it if possible":
HAVE_PSYCO = None
CHECK_PSYCO = True

#
# plug-ins
# ========
#
# math
# ----
# Python standard float/comple types and math/cmath modules only
# provide fixed precision (double precision), but sometimes
# multiprecision floating point is needed.
# If you have mpmath installed, set HAVE_MPMATH True and use:
#PLUGIN_MATH = 'mpmath'
#CHECK_PLUGIN_MATH = False
# Otherwise, use only Python float/complex as default (but use mpmath
# if possible):
PLUGIN_MATH = None
CHECK_PLUGIN_MATH = True

#
# Assumptions
# ===========
#
# Some conjectures are useful for assuring the validity of a fast
# algorithm.
#
# All assumptions are default to False, but you can set them True if
# you believe them.
#
# GRH
# ---
#
# Generalized Riemann Hypothesis.
# For example, primality test is O((log n)**2) if GRH is true
# while O((log n)**6) or something without it.
#
# If you believe GRH as the truth:
#GRH = True
# The default is, of course, conservatively doubting it:
GRH = False

#
# Files
# =====
#
# data directory
# --------------
#
# The directory where nzmath (static) data files are stored.
# The default will be os.path.join(sys.prefix, 'share', 'nzmath')
# or os.path.join(sys.prefix, 'Data', 'nzmath') on Windows.
#
# If your *nix computer installs NZMATH as a system program:
#DATADIR = '/usr/share/nzmath'
#CHECK_DATADIR = False
# If it is an optional program:
#DATADIR = '/usr/local/share/nzmath'
#CHECK_DATADIR = False
# Windows users may have more aggressive path:
#DATADIR = r'C:\Python25\Data'
#CHECK_DATADIR = False
# The default values mean "I don't know; check it later":
DATADIR = None
CHECK_DATADIR = True


# -------------------
# User Configurations
# -------------------

confdir = os.environ.get('NZMATHCONFDIR', None)
if confdir is None:
    if sys.platform in WINDOWS_PLATFORMS:
        # "C:\Documents and Settings\%USERNAME%\Application Data\nzmath"
        # APPDIR = "C:\Documents and Settings\%USERNAME%\Application Data"
        appdir = os.environ.get('APPDATA', None)
        # USERPROFILE = "C:\Documents and Settings\%USERNAME%"
        profdir = os.environ.get('USERPROFILE', None)
        if appdir is not None:
            t_confdir = os.path.join(appdir, 'nzmath')
        elif profdir is not None:
            t_confdir = os.path.join(profdir, 'Application Data', 'nzmath')
    else:
        # "~/.nzmath.d/"
        homedir = os.environ.get('HOME', None)
        if homedir is not None:
            t_confdir = os.path.join(homedir, '.nzmath.d')

    # nzmathconf.py exists?
    if os.path.exists(os.path.join(t_confdir, 'nzmathconf.py')):
        confdir = t_confdir

if confdir is not None and os.path.exists(confdir):
    sys.path.insert(0, confdir)
else:
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config'))

try:
    # overwrite the default settings with user's nzmathconf
    from nzmathconf import *
except ImportError:
    warnings.warn("nzmathconf.py not found")


# ------
# Checks 
# ------
#
# mpmath
# ------
#

def check_mpmath():
    """
    Check if mpmath is importable or not
    """
    try:
        import mpmath
        return True
    except ImportError:
        return False

if CHECK_MPMATH:
    HAVE_MPMATH = check_mpmath()

#
# sqlite3
# -------
#

def check_sqlite3():
    """
    Check if sqlite3 is importable or not.
    pysqlite2 may be a substitution.
    """
    try:
        try:
            import sqlite3
            return True
        except ImportError:
            import pysqlite2.dbapi2 as sqlite3
            return True
    except ImportError:
        return False

if CHECK_SQLITE3:
    HAVE_SQLITE3 = check_sqlite3()
    

#
# net availability
# ----------------
#
def check_net():
    """
    Check the net connection by http call.
    """
    import urllib2
    try:
        urllib2.urlopen('http://sourceforge.net/projects/nzmath/')
        return True
    except urllib2.HTTPError:
        # the problem is on server side, thus connected anyway
        return True
    except urllib2.URLError:
        # no dns, thus no connection
        return False
    except Exception:
        # I don't know the reason, but something wrong
        return False

if CHECK_NET:
    HAVE_NET = check_net()

#
# psyco availability
# ------------------
#
def check_psyco():
    """
    Check if psyco is importable or not.
    """
    try:
        import psyco
        return True
    except ImportError:
        return False

if CHECK_PSYCO:
    HAVE_PSYCO = check_psyco()
        
#
# math plug-in
#
def check_plugin_math():
    """
    Return 'mpmath' if HAVE_MPMATH, None otherwise.
    """
    if HAVE_MPMATH:
        return 'mpmath'
    else:
        return None

if CHECK_PLUGIN_MATH:
    PLUGIN_MATH = check_plugin_math()


#
# data directory
# --------------
#
def default_datadir():
    candidates = []
    
    # developing environment
    candidates.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'dist', 'data'))

    if DATADIR is not None:
        candidates.append(DATADIR)
    if sys.platform in WINDOWS_PLATFORMS:
        candidates.append(os.path.join(sys.prefix, 'Data', 'nzmath'))
    else:
        candidates.append(os.path.join(sys.prefix, 'share', 'nzmath'))
    # default install script will install here ...hmm
    candidates.append(os.path.join(sys.prefix, 'nzmath'))

    for canddir in candidates:
        if os.path.exists(canddir):
            return canddir
    return None

if CHECK_DATADIR:
    DATADIR = default_datadir()
    if DATADIR is None:
        warnings.warn('no datadir found')


# Declare exported variables.

__all__ = ['HAVE_MPMATH', 'HAVE_SQLITE3', 'HAVE_NET', 'HAVE_PSYCO',
           'PLUGIN_MATH', 'GRH', 'DATADIR']
