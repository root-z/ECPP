"""
Final program which combines Atkin-Morain ECPP with Miller-Rabin Primality test.
"""

from miller_rabin import miller_rabin
from ecpp import atkin_morain
import sys


def prime(n):
    """

    Args:
        n: Number to be tested

    Returns:
        certificate if the number is prime, False otherwise.
    """
    if not n > 0:
        raise ValueError("input must be greater than 0")
    repeat = 50
    if miller_rabin(n, 50):
        cert = atkin_morain(n)
        if not cert:
            raise InconsistencyError("ECPP")
        else:
            return cert
    else:
        return False


class InconsistencyError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self)

if __name__=='__main__':
    n = int(sys.argv[1])
    try:
        cert = prime(n)
        if not cert:
            print str(n) + " is composite."
        else:
            print str(n) + " is prime and the certificate is " + str(cert)
    except ValueError:
        print "Input must be greater than 0."
    except InconsistencyError:
        print "It's a probable prime yet no certificate was found. Please try again."
