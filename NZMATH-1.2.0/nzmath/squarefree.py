"""
Squarefreeness tests.

Definition:
  n: squarefree <=> there is no p whose square divides n.

Examples:
  - 0 is non-squarefree because any square of prime can divide 0.
  - 1 is squarefree because there is no prime dividing 1.
  - 2, 3, 5, and any other primes are squarefree.
  - 4, 8, 9, 12, 16 are non-squarefree composites.
  - 6, 10, 14, 15, 21 are squarefree composites.
"""

import math
import nzmath.arith1 as arith1
import nzmath.bigrange as bigrange
import nzmath.prime as prime
import nzmath.rational as rational
import nzmath.factor.methods as factor_methods
import nzmath.factor.misc as factor_misc


class Undetermined (Exception):
    """
    Undetermined state of calculation.
    """


def lenstra(n):
    """
    If return value is True, n is squarefree.  Otherwise, the
    squarefreeness is still unknown and Undetermined is raised.

    The condition is so strong that it seems n is a prime or a
    Carmichael number.

    pre-condition: n & 1
    reference: H.W.Lenstra 1973 ---
    """
    n = int(n) # see sf bug #1826712
    predn = n - 1
    bound = int(math.log(n)**2 + 1)
    for i in range(2, bound):
        if pow(i, predn, n) != 1:
            raise Undetermined("Lenstra's method can't determine squarefreeness")
    return True


def trial_division(n):
    """
    Test whether n is squarefree or not.

    The method is a kind of trial division.
    """
    try:
        return trivial_test(n)
    except Undetermined:
        pass

    for p in prime.generator():
        if not (n % (p*p)):
            # found a square factor
            return False
        elif not (n % p):
            # found a non-square factor
            n //= p
            try:
                return trivial_test(n)
            except Undetermined:
                pass
        if p*p*p > n:
            break
    # At the end of the loop:
    #   n doesn't have any factor less than its cubic root.
    #   n is not a prime nor a perfect square number.
    # The factor must be two primes p and q such that p < sqrt(n) < q.
    return True


def trivial_test(n):
    """
    Test whether n is squarefree or not.

    This method do anything but factorization.
    """
    if n == 1 or n == 2:
        return True
    if arith1.issquare(n):
        return False
    if n & 1:
        return lenstra(n)
    elif not (n & 3):
        return False
    raise Undetermined("trivial test can't determine squarefreeness")


def viafactor(n):
    """
    Test whether n is squarefree or not.

    It is obvious that if one knows the prime factorization of the number,
    he/she can tell whether the number is squarefree or not.
    """
    for p, e in factor_misc.FactoredInteger(n):
        if e >= 2:
            return False
    return True


# ternary logic versions
#
# The third logical value means "uncertain" or "proof unknown".
# We designate None to this value.  It lets "unknown" status
# be, at least, not true.
# There are nothing corresponding to boolean logic operators.
# They are out of scope of this module.

def lenstra_ternary(n):
    """
    Test the squarefreeness of n.
    The return value is one of the ternary logical constants.
    If return value is TRUE, n is squarefree.  Otherwise, the
    squarefreeness is still unknown and UNKNOWN is returned.

    The condition is so strong that it seems n is a prime or a
    Carmichael number.

    pre-condition: n & 1
    reference: H.W.Lenstra 1973 ---
    """
    n = int(n) # see sf bug #1826712
    predn = n - 1
    bound = int(math.log(n)**2 + 1)
    for i in range(2, bound):
        if pow(i, predn, n) != 1:
            return None
    return True


def trivial_test_ternary(n):
    """
    Test the squarefreeness of n.
    The return value is one of the ternary logical constants.

    The method uses a series of trivial tests.
    """
    if n == 1 or n == 2:
        return True
    if arith1.issquare(n):
        return False
    if n & 1:
        return lenstra_ternary(n)
    elif not (n & 3):
        return False
    return None


def trial_division_ternary(n):
    """
    Test the squarefreeness of n.
    The return value is one of the True or False, not None.

    The method is a kind of trial division.
    """
    result = trivial_test_ternary(n)
    if result is not None:
        return result

    for p in prime.generator():
        if not (n % (p*p)):
            # found a square factor
            return False
        elif not (n % p):
            # found a non-square factor
            n //= p
            result = trivial_test_ternary(n)
            if result is not None:
                return result
        if p*p*p > n:
            break
    # At the end of the loop:
    #   n doesn't have any factor less than its cubic root.
    #   n is not a prime nor a perfect square number.
    # The factor must be two primes p and q such that p < sqrt(n) < q.
    return True


# Just for symmetry, viafactor_ternary is defined as an alias of viafactor.
viafactor_ternary = viafactor


class SquarefreeDecompositionMethod (factor_methods.TrialDivision):
    """
    Decomposition of an integer into square part and squarefree part.
    """
    def __init__(self):
        factor_methods.TrialDivision.__init__(self)
        self.primeseq = None # initialized later

    def generate(self, target, **options):
        """
        Generate squarefree factors of the target number with their
        valuations.  The method may terminate with yielding (1, 1)
        to indicate the factorization is incomplete.

        If a keyword option 'strict' is False (default to True),
        factorization will stop after the first square factor no
        matter whether it is squarefree or not.
        """
        strict = options.get('strict', True)
        options['n'] = target
        primeseq = self._parse_seq(options)
        for p in primeseq:
            if not (target % p):
                e, target = arith1.vp(target, p)
                yield p, e
                if target == 1:
                    break
                elif e > 1 and not strict:
                    yield 1, 1
                    break
                elif trivial_test_ternary(target):
                    # the factor remained is squarefree.
                    yield target, 1
                    break
                q, e = factor_misc.primePowerTest(target)
                if e:
                    yield q, e
                    break
                sqrt = arith1.issquare(target)
                if sqrt:
                    if strict:
                        for q, e in self.factor(sqrt, iterator=primeseq):
                            yield q, 2 * e
                    else:
                        yield sqrt, 2
                    break
            if p ** 3 > target:
                # there are no more square factors of target,
                # thus target is squarefree
                yield target, 1
                break
        else:
            # primeseq is exhausted but target has not been proven prime
            yield 1, 1

    def issquarefree(self, n):
        """
        Return True if n is squarefree, False otherwise.

        The method uses partial factorization into squarefree parts,
        if such partial factorization is possible.  In other cases,
        It completely factor n by trial division.
        """
        if trivial_test_ternary(n):
            return True
        for s, e in self.generate(n, strict=False):
            if e > 1:
                return False
        return True


viadecomposition = SquarefreeDecompositionMethod().issquarefree
