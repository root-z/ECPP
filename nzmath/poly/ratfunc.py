"""
rational functions
"""

from __future__ import division
import nzmath.ring as ring
import nzmath.poly.ring as poly_ring


class RationalFunction(ring.QuotientFieldElement):
    """
    The class of rational functions.
    """
    def __init__(self, *arg, **kwd):
        if len(arg) == 1:
            if isinstance(arg[0], RationalFunction):
                numerator = arg[0].numerator
                denominator = arg[0].denominator
            else:
                numerator = arg[0]
                if "denominator" in kwd:
                    denominator = kwd["denominator"]
                else:
                    denominator = self.numerator.getRing().one
        elif len(arg) == 2:
            numerator = arg[0]
            denominator = arg[1]
        elif len(kwd) > 0:
            if "numerator" in kwd:
                numerator = kwd["numerator"]
            else:
                raise ValueError("numerator must be specified.")
            if "denominator" in kwd:
                denominator = kwd["denominator"]
            else:
                denominator = self.numerator.getRing().one
        else:
            raise ValueError("numerator must be specified.")
        ring.QuotientFieldElement.__init__(self, numerator, denominator)
        if self.numerator.number_of_variables == self.denominator.number_of_variables:
            self.number_of_variables = self.numerator.number_of_variables
        else:
            raise TypeError("numerator and denominator are inconsistent")

    def __eq__(self, other):
        """
        equality test
        """
        try:
            return ring.QuotientFieldElement.__eq__(self, other)
        except AttributeError:
            return NotImplemented

    def __hash__(self):
        try:
            return ring.QuotientFieldElement.__hash__(self)
        except AttributteError:
            return NotImplemented

    def __call__(self, *args):
        """
        evaluation

        The type of args depends on the type of polynomials
        representing numerator and denominator.
        """
        return self.numerator(*args) / self.denominator(*args)

    def __str__(self):
        """
        Return a simple string
        """
        return str(self.numerator) + " / " + str(self.denominator)

    def __repr__(self):
        """
        Return a string representation.
        """
        return "%s(%s, %s)" % (self.__class__.__name__, repr(self.numerator), repr(self.denominator))

    def getRing(self):
        """
        Return a ring to which the rational function belongs.
        """
        nring = self.numerator.getCoefficientRing()
        if not nring.isfield():
            nring = nring.getQuotientField()
        return poly_ring.RationalFunctionField(nring, self.number_of_variables)
