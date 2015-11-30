# -------------
# NZMATH config
# -------------
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
# If you have net connectivity now, set as the following:
#HAVE_NET = True
#CHECK_NET = False
# Or, if your machine is not connected, set as the following:
#HAVE_NET = False
#CHECK_NET = False
# The default values mean "I don't know; check it later":
HAVE_NET = None
CHECK_NET = True

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
