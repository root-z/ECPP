"""
Term Order for polynomials.
"""

import re
import nzmath.compatibility
import nzmath.ring as ring


_INTERFACE_MSG = "%s is interface"

class TermOrderInterface (object):
    """
    (abstract term order)

    A term order is primalily a function, which determines precedence
    between two terms (or monomials).  By the precedence, all terms
    are ordered.

    More precisely in terms of Python, a term order accepts two tuples
    of integers, each of which represents power indices of the term,
    and returns 0, 1 or -1 just like cmp built-in function.

    A TermOrder object provides not only the precedence function, but
    also a method to format a string for a polynomial, to tell degree,
    leading coefficients, etc.
    """

    _PLUS_MINUS = re.compile(r"(^-|\+ -)")

    def __init__(self, comparator):
        """
        'comparator' accepts two tuples of integers, each of which
        represents power indices of the term, and returns 0, 1 or -1
        just like cmp built-in function.
        """
        if type(self) is TermOrderInterface:
            raise NotImplementedError(_INTERFACE_MSG % self.__class__.__name__)
        self.comparator = comparator

    def cmp(self, left, right):
        """
        Compare two indices left and right and determine precedence by
        self.comparator.
        """
        return self.comparator(left, right)

    def bisect(self, array, elem, lo=0, hi=None):
        """
        Return the index where to insert item `elem' in a list `array',
        assuming array is sorted with the order.  (In the following,
        a refers `array' and x refers `elem')

        The return value i is such that all e in a[:i] have e <= x, and
        all e in a[i:] have e > x.  So if x already appears in the list,
        a.insert(x) will insert just after the rightmost x already there.

        Optional args lo (default 0) and hi (default len(a)) bound the
        slice of a to be searched.

        This function is based on the bisect.bisect_right of the Python
        standard library.
        """
        if hi is None:
            hi = len(array)
        while lo < hi:
            mid = (lo + hi) >> 1
            if self.cmp(elem, array[mid]) < 0:
                hi = mid
            else: lo = mid + 1
        return lo

    def format(self, polynom, **kwds):
        """
        Return the formatted string of the polynomial.
        """
        return self.comparator(left, right)

    def bisect(self, array, elem, lo=0, hi=None):
        """
        Return the index where to insert item `elem' in a list `array',
        assuming array is sorted with the order.  (In the following,
        a refers `array' and x refers `elem')

        The return value i is such that all e in a[:i] have e <= x, and
        all e in a[i:] have e > x.  So if x already appears in the list,
        a.insert(x) will insert just after the rightmost x already there.

        Optional args lo (default 0) and hi (default len(a)) bound the
        slice of a to be searched.

        This function is based on the bisect.bisect_right of the Python
        standard library.
        """
        if hi is None:
            hi = len(array)
        while lo < hi:
            mid = (lo + hi) >> 1
            if self.cmp(elem, array[mid]) < 0:
                hi = mid
            else: lo = mid + 1
        return lo

    def leading_coefficient(self, polynom):
        """
        Return the leading coefficient of polynomial 'polynom' with
        respect to the term order.
        """
        raise NotImplementedError(_INTERFACE_MSG % self.__class__.__name__)

    def leading_term(self, polynom):
        """
        Return the leading term of polynomial 'polynom' as tuple of
        (degree index, coefficient) with respect to the term order.
        """
        raise NotImplementedError(_INTERFACE_MSG % self.__class__.__name__)


class UnivarTermOrder (TermOrderInterface):
    """
    term order for univariate polynomials.

    One thing special to univariate case is that powers are not tuples
    but integers.
    """
    def __init__(self, comparator):
        """
        UnivarTermOrder(comparator)

        'comparator' can be any callable that accepts two integers and
        returns 0, 1 or -1 just like cmp, i.e. if they are equal it
        returns 0, first one is greater 1, and otherwise -1.
        Theoretically acceptable comparator is only the cmp function.
        """
        TermOrderInterface.__init__(self, comparator)

    def format(self, polynom, varname="X", reverse=False):
        """
        Return the formatted string of the polynomial.

        - 'polynom' must be a univariate polynomial.
        - 'varname' can be set to the name of the variable (default to
          'X').
        - 'reverse' can be either True or False. If it's True, terms
          appear in reverse (descending) order.
        """
        degrees = _sort_with_cmp(
            [base for base in polynom.iterbases()],
            self.comparator)
        if reverse:
            degrees.reverse()
        str_terms = [("%s * %s ** %d" % (polynom[d], varname, d)) for d in degrees if polynom[d]]
        # constant
        if 0 in degrees and polynom[0]:
            const_term = str(polynom[0])
            if (hasattr(polynom, "getCoefficientRing") and
                polynom[0] == polynom.getCoefficientRing().one):
                const_term = "1"
            str_terms[str_terms.index("%s * %s ** 0" % (polynom[0], varname))] = const_term
        # degree 1
        if 1 in degrees and polynom[1]:
            str_terms[str_terms.index("%s * %s ** 1" % (polynom[1], varname))] = "%s * %s" % (polynom[1], varname)
        result = " + ".join(str_terms)
        # minus terms
        result = self._PLUS_MINUS.sub("- ", result)
        # coefficient is 1 (or -1)
        if hasattr(polynom, "getCoefficientRing"):
            one_times_x = re.compile(r"(^| )%s \* %s" % (polynom.getCoefficientRing().one, varname))
        else:
            one_times_x = re.compile(r"(^| )1 \* %s" % varname)
        result = one_times_x.sub(" " + varname, result)
        result = result.lstrip()
        return result

    def degree(self, polynom):
        """
        Return the degree of the polynomial 'polynom'.
        """
        if hasattr(polynom, "degree"):
            return polynom.degree()
        degree = -1
        for d in polynom.iterbases():
            if self.comparator(degree, d) < 0:
                degree = d
        return degree

    def leading_coefficient(self, polynom):
        """
        Return the leading coefficient of polynomial 'polynom' with
        respect to the term order.
        """
        if hasattr(polynom, 'leading_coefficient'):
            return polynom.leading_coefficient()
        degree, lc = -1, 0
        for d, c in polynom:
            if self.comparator(degree, d) < 0:
                degree, lc = d, c
        return lc

    def leading_term(self, polynom):
        """
        Return the leading term of polynomial 'polynom' as tuple of
        (degree, coefficient) with respect to the term order.
        """
        if hasattr(polynom, 'leading_term'):
            return polynom.leading_term()
        degree, lc = -1, 0
        for d, c in polynom:
            if self.comparator(degree, d) < 0:
                degree, lc = d, c
        return degree, lc

    def tail_degree(self, polynom):
        """
        Return the least degree among all terms of the polynomial
        'polynom'.

        This method is EXPERIMENTAL.
        """
        if hasattr(polynom, "tail_degree"):
            return polynom.tail_degree()
        degree = -1
        for d in polynom.iterbases():
            if degree == -1 or self.comparator(degree, d) > 0:
                degree = d
                if degree == 0:
                    break
        return degree


ascending_order = UnivarTermOrder(cmp)


class MultivarTermOrder (TermOrderInterface):
    """
    A class of term orders for multivariate polynomials.
    """
    def __init__(self, comparator):
        """
        'comparator' accepts two tuples of integers, each of which
        represents power indices of the term, and returns 0, 1 or -1
        just like cmp built-in function.
        """
        self.comparator = comparator

    def format(self, polynom, varnames=None, reverse=False, **kwds):
        """
        Return the formatted string of the polynomial.

        An additional keyword argument 'varnames' is required to name
        variables.
        """
        if varnames is None:
            raise TypeError("keyword argument 'varnames' is required")

        bases = _sort_with_cmp(polynom.bases(), self.comparator)
        if reverse:
            bases.reverse()

        result = " + ".join([self._format_term((base, polynom[base]), varnames) for base in bases if polynom[base]])
        # minus terms
        result = self._PLUS_MINUS.sub("- ", result)
        # coefficient is 1 (or -1)
        if hasattr(polynom, "getCoefficientRing"):
            one_times = re.compile(r"(^| )%s \* " % polynom.getCoefficientRing().one)
        else:
            one_times = re.compile(r"(^| )1 \* ")
        result = one_times.sub(" ", result)

        result = result.lstrip()
        if not result:
            result = "0"
        return result

    def _format_term(self, term, varnames):
        """
        Return formatted term string.

        'term' is a tuple of indices and coefficient.
        """
        if not term[1]:
            return ""
        if term[1] == ring.getRing(term[1]).one:
            powlist = []
        else:
            powlist = [str(term[1])]
        for v, d in zip(varnames, term[0]):
            if d > 1:
                powlist.append("%s ** %d" % (v, d))
            elif d == 1:
                powlist.append(v)

        if not powlist:
            # coefficient is plus/minus one and every variable has degree 0
            return str(term[1])
        return " * ".join(powlist)

    def leading_coefficient(self, polynom):
        """
        Return the leading coefficient of polynomial 'polynom' with
        respect to the term order.
        """
        if hasattr(polynom, 'leading_coefficient'):
            return polynom.leading_coefficient()
        return polynom[self._max(polynom.bases())]

    def leading_term(self, polynom):
        """
        Return the leading term of polynomial 'polynom' as tuple of
        (degree index, coefficient) with respect to the term order.
        """
        if hasattr(polynom, 'leading_term'):
            return polynom.leading_term()
        max_indices = self._max(polynom.bases())
        return max_indices, polynom[max_indices]

    def _max(self, indices_list):
        """
        Return the maximum indices with respect to the comparator.
        """
        if not indices_list:
            raise ValueError("max() arg is an empty sequence")
        it = iter(indices_list)
        maxi = it.next()
        for indices in it:
            if self.comparator(maxi, indices) < 0:
                maxi = indices
        return maxi


def weight_order(weight, tie_breaker=None):
    """
    Return a comparator of weight ordering.

    The weight ordering is defined for arguments x and y that x < y
    if w.x < w.y (dot products) or w.x == w.y and tie breaker tells
    x < y.

    The option 'tie_breaker' is another comparator that will be used
    if dot products with the weight vector leaves arguments tie.  If
    the option is None (default) and a tie breaker is indeed necessary
    to order given arguments, a TypeError is raised.
    """
    def _order(mono, mial):
        if mono == mial:
            return 0
        wx = sum(w * x for (w, x) in zip(weight, mono))
        wy = sum(w * y for (w, y) in zip(weight, mial))
        if wx != wy:
            return cmp(wx, wy)
        return tie_breaker(mono, mial)
    return _order


def _total_degree_lexicographic(left, right):
    """
    Total degree lexicographic (or graded lexicographic) term order :
      L < R iff
      (1) sum(li) < sum(ri) or
      (2) sum(li) = sum(ri) and
          there exists i s.t. l0 == r0, ..., l(i-1) == r(i-1), li < ri.
    """
    sum_left, sum_right = sum(left), sum(right)
    if sum_left != sum_right:
        return cmp(sum_left, sum_right)
    return cmp(left, right)

def _total_degree_reverse_lexicographic(left, right):
    """
    Total degree reverse lexicographic (or graded reverse
    lexicographic) term order :
      L < R iff
      (1) sum(li) < sum(ri) or
      (2) sum(li) = sum(ri) and
          there exists i s.t. ln == rn, ..., l(i+1) == r(i+1), li > ri.
    """
    sum_left, sum_right = sum(left), sum(right)
    if sum_left != sum_right:
        return cmp(sum_left, sum_right)
    return cmp(right[::-1], left[::-1])


lexicographic_order = MultivarTermOrder(cmp)
total_degree_lexicographic_order = MultivarTermOrder(_total_degree_lexicographic)
total_degree_reverse_lexicographic_order = MultivarTermOrder(_total_degree_reverse_lexicographic)


### for compatibility function
def _sort_with_cmp(seq, mycmp):
    """
    Return the sorted seq by using the comparator function mycmp.
    """
    
    try:
        return sorted(seq, cmp=mycmp)
    except TypeError: # cmp is no longer available
        class internal_cmp:
            def __init__(self, obj):
                self.obj = obj
            def __gt__(self, other):
                return mycmp(self.obj, other.obj) > 0
            def __ge__(self, other):
                return mycmp(self.obj, other.obj) >= 0
            def __lt__(self, other):
                return mycmp(self.obj, other.obj) < 0
            def __le__(self, other):
                return mycmp(self.obj, other.obj) <= 0
            def __cmp__(self, other):
                return mycmp(self.obj, other.obj)
        return sorted(seq, key=internal_cmp)

