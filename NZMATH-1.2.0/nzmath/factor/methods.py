"""
factoring methods.
"""

import nzmath.arith1 as arith1
import nzmath.bigrange as bigrange
import nzmath.prime as prime
import nzmath.factor.util as util
import nzmath.factor.find as find
from nzmath.factor.mpqs import mpqsfind
from nzmath.factor.ecm import ecm as ecmfind


class DefaultMethod (util.FactoringMethod):
    """
    A factor method used as the default.

    It tries the trial division method first, then the p-1 method,
    and finally calls the MPQS.
    """
    def __init__(self):
        util.FactoringMethod.__init__(self)

    def factor(self, number, **options):
        """
        Factor the given positive integer.
        The returned value is in the form of [(p1, e1), ..., (pn, en)].
        """
        if not self._validate_input_number(number):
            return []

        # backup
        original_return_type = options.get('return_type', '')

        # trial division first
        trial = TrialDivision()
        trial.verbose = self.verbose
        if number < 1000000:
            return trial.factor(number, **options)
        options['return_type'] = 'tracker'
        options['eratosthenes'] = options.get('eratosthenes', 100000)
        tracker = trial.factor(number, **options)

        # try p-1 method
        pmom = PMinusOneMethod()
        pmom.verbose = self.verbose
        tracker = pmom.continue_factor(tracker, **options)

        # finally mpqs
        options['return_type'] = original_return_type
        mpqs = MPQSMethod()
        mpqs.verbose = self.verbose
        result = mpqs.continue_factor(tracker, **options)
        if options['return_type'] == 'list' and options.get('need_sort', False):
            result.sort()
        return result


class TrialDivision (util.FactoringMethod):
    """
    Class for trial division method.
    """

    def __init__(self):
        util.FactoringMethod.__init__(self)

    def factor(self, number, **options):
        """
        Factor the given integer by trial division.

        options for the trial sequence can be either:
        1) 'start' and 'stop' as range parameters.
        2) 'iterator' as an iterator of primes.
        3) 'eratosthenes' as an upper bound to make prime sequence by sieve.
        If none of the options above are given, prime factor is
        searched from 2 to the square root of the given integer.

        an option 'return_type' is for the returned type, whose value can be:
        1) 'list' for default type [(p1, e1), ..., (pn, en)].
        2) 'tracker' for alternative type FactoringInteger.
        """
        # to be used in _parse_seq
        options['n'] = number

        return util.FactoringMethod.factor(self, number, **options)

    def continue_factor(self, tracker, **options):
        """
        Continue factoring and return the result of factorization.

        The argument 'tracker' should be an instance of FactoringInteger.
        The returned type is FactoringInteger.

        options is the same for factor, but the default value for
        'return_type' is 'tracker'.
        """
        return util.FactoringMethod.continue_factor(self, tracker, **options)

    def _parse_seq(self, options):
        """
        Parse 'options' to define trial sequaence.
        """
        if 'start' in options and 'stop' in options:
            if 'step' in options:
                trials = bigrange.range(options['start'], options['stop'], options['step'])
            else:
                trials = bigrange.range(options['start'], options['stop'])
        elif 'iterator' in options:
            trials = options['iterator']
        elif 'eratosthenes' in options:
            trials = prime.generator_eratosthenes(options['eratosthenes'])
        elif options['n'] < 1000000:
            trials = prime.generator_eratosthenes(arith1.floorsqrt(options['n']))
        else:
            trials = prime.generator()
        return trials

    def generate(self, target, **options):
        """
        Generate prime factors of the target number with their
        valuations.  The method may terminate with yielding (1, 1)
        to indicate the factorization is incomplete.
        """
        primeseq = self._parse_seq(options)
        for p in primeseq:
            if not (target % p):
                e, target = arith1.vp(target // p, p, 1)
                yield p, e
            if p ** 2 > target:
                # there are no more factors of target, thus target is a prime
                yield target, 1
                break
        else:
            # primeseq is exhausted but target has not been proven prime
            yield 1, 1


class PMinusOneMethod (util.FactoringMethod):
    """
    Class for p-1 method.
    """

    def __init__(self):
        util.FactoringMethod.__init__(self)

    def find(self, target, **options):
        """
        Find a factor from the target number.
        """
        return find.pmom(target, verbose=self.verbose)


class RhoMethod (util.FactoringMethod):
    """
    Class for Pollard's rho method.
    """

    def __init__(self):
        util.FactoringMethod.__init__(self)

    def find(self, target, **options):
        """
        Find a factor from the target number.
        """
        return find.rhomethod(target, verbose=self.verbose)


class MPQSMethod (util.FactoringMethod):
    """
    Class for Multi-Polynomial Quadratic Sieve method.
    """

    def __init__(self):
        util.FactoringMethod.__init__(self)

    def find(self, target, **options):
        """
        Find a factor from the target number.
        """
        limited_options = {}
        if 's' in options:
            limited_options['s'] = options['s']
        if 'f' in options:
            limited_options['f'] = options['f']
        if 'm' in options:
            limited_options['m'] = options['m']
        if 'verbose' in options:
            limited_options['verbose'] = options['verbose']
        return mpqsfind(target, **limited_options)


class EllipticCurveMethod (util.FactoringMethod):
    """
    Class for Elliptic Curve Method
    """
    def __init__(self):
        util.FactoringMethod.__init__(self)

    def find(self, target, **options):
        """
        Find a factor from the target number.
        """
        return ecmfind(target, **options)


def trialDivision(n, **options):
    """
    Factor the given integer by trial division.

    options for the trial sequence can be either:
    1) 'start' and 'stop' as range parameters.
    2) 'iterator' as an iterator of primes.
    3) 'eratosthenes' as an upper bound to make prime sequence by sieve.
    If none of the options above are given, prime factor is
    searched from 2 to the square root of the given integer.
    """
    method = TrialDivision()
    options['return_type'] = 'list'
    return method.factor(n, **options)

def pmom(n, **options):
    """
    Factor the given integer by Pollard's p-1 method.
    """
    method = PMinusOneMethod()
    options['return_type'] = 'list'
    options['need_sort'] = True
    return method.factor(n, **options)

def rhomethod(n, **options):
    """
    Factor the given integer by rho method.
    """
    method = RhoMethod()
    options['return_type'] = 'list'
    options['need_sort'] = True
    return method.factor(n, **options)

def mpqs(n, **options):
    """
    Factor the given integer by multi-polynomial quadratic sieve method.
    """
    method = MPQSMethod()
    options['return_type'] = 'list'
    options['need_sort'] = True
    return method.factor(n, **options)

def ecm(n, **options):
    """
    Factor the given integer by elliptic curve method.
    """
    method = EllipticCurveMethod()
    options['return_type'] = 'list'
    options['need_sort'] = True
    return method.factor(n, **options)

def factor(n, method='default', **options):
    """
    Factor the given integer.

    By default, use several methods internally.
    Optional argument 'method' can be:
      'ecm': use EllipticCurveMethod
      'mpqs': use MPQSMethod.
      'pmom': use PMinusOneMethod.
      'rhomethod': use RhoMethod.
      'trialDivision': use TrialDivision.
    (In fact, its initial letter suffice to specify.)
    """
    method_table = {'d': DefaultMethod,
                    'e': EllipticCurveMethod,
                    'm': MPQSMethod,
                    'p': PMinusOneMethod,
                    'r': RhoMethod,
                    't': TrialDivision}
    try:
        chosen_method = method_table[method[0]]()
    except KeyError:
        chosen_method = DefaultMethod()
    options['return_type'] = 'list'
    options['need_sort'] = True
    return chosen_method.factor(n, **options)
