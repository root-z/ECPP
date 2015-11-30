from __future__ import division

import math
import nzmath.arith1 as arith1
import nzmath.equation as equation
import nzmath.gcd as gcd
import nzmath.lattice as lattice
import nzmath.vector as vector
import nzmath.matrix as matrix
import nzmath.factor.misc as misc
import nzmath.finitefield as finitefield
import nzmath.poly.uniutil as uniutil
import nzmath.poly.multiutil as multiutil
import nzmath.rational as rational
import nzmath.ring as ring
import nzmath.squarefree as squarefree
import nzmath.round2 as round2


class NumberField (ring.Field):
    """
    A class of number field.
    """
    def __init__(self, polynomial, precompute=False):
        """
        Initialize a number field with given polynomial coefficients
        (in ascending order).
        """
        ring.Field.__init__(self)
        self.polynomial = polynomial
        self.degree = len(polynomial) - 1
        if precompute:
            getConj()
            disc()
            signature()
            integer_ring()

    def __repr__(self):
        return_str = '%s(%s)' % (self.__class__.__name__, self.polynomial)
        return return_str

    def __mul__(self, other):
        """
        Output composite field of self and other.
        """
        common_options = {"coeffring": rational.theIntegerRing,
                          "number_of_variables": 2}
        flist = [((d, 0), c) for (d, c) in enumerate(self.polynomial)]
        f = multiutil.polynomial(flist, **common_options)
        g = zpoly(other.polynomial)
        diff = multiutil.polynomial({(1, 0): 1, (0, 1):-1}, **common_options)
        compos = f.resultant(g(diff), 0)
        return NumberField([compos[i] for i in range(compos.degree() + 1)])

    def getConj(self):
        """
        Return (approximate) solutions of self.polynomial.
        We can discriminate the conjugate field of self by these values.
        """
        if not hasattr(self, "conj"):
            conj = equation.SimMethod(self.polynomial)
            self.conj = conj
        return self.conj

    def disc(self):
        """
        Compute the discriminant of self.polynomial.
        (The output is not field disc of self but disc of self.polynomial.)
        """
        """
        Obasis, disc = round2.round2(self)
        A = [algfield.MatAlgNumber(obasis[0][i]) for i in range(degree)]
        return disc(A) <--- real disc of self field.
        """

        if not hasattr(self, "poldisc"):
            degree = self.degree
            traces = []
            for i in range(degree):
                for j in range(degree):
                    s = self.basis(i)*self.basis(j)
                    traces.append(s.trace())
            M = matrix.RingSquareMatrix(degree, degree, traces)
            self.poldisc = M.determinant()
        return self.poldisc

    def integer_ring(self):
        """
        Return the integer ring of self (using round2 module)
        """
        if not hasattr(self, "integer_basis"):
            int_basis, self.discriminant = round2.round2(self.polynomial)
            self.integer_basis = matrix.createMatrix(
                self.degree, [vector.Vector(int_basis[j]) for j in range(
                len(int_basis))])
        return self.integer_basis

    def field_discriminant(self):
        """
        Return the (field) discriminant of self (using round2 module)
        """
        if not hasattr(self, "discriminant"):
            int_basis, self.discriminant = round2.round2(self.polynomial)
            self.integer_basis = matrix.createMatrix(
                self.degree, [vector.Vector(int_basis[j]) for j in range(
                len(int_basis))])
        return self.discriminant

    def basis(self, j):
        """
        Return the j-th basis of self.
        """
        base = [0] * self.degree
        base[j] = 1
        return BasicAlgNumber([base, 1], self.polynomial)

    def signature(self):
        """
        Using Strum's algorithm, compute the signature of self.
        Algorithm 4.1.11 in Cohen's Book
        """
        if not hasattr(self, "sign"):
            degree = self.degree
            #Step 1.
            if degree == 0:
                return (0, 0)
            # no check for degree 1?

            minpoly = zpoly(self.polynomial)
            d_minpoly = minpoly.differentiate()
            A = minpoly.primitive_part()
            B = d_minpoly.primitive_part()
            g = 1
            h = 1
            pos_at_inf = A.leading_coefficient() > 0
            pos_at_neg = pos_at_inf == (degree & 1)
            r_1 = 1

            #Step 2.
            while True:
                deg = A.degree() - B.degree()
                residue = A.pseudo_mod(B)
                if not residue:
                    raise ValueError("not squarefree")
                if B.leading_coefficient() > 0 or deg & 1:
                    residue = - residue
                #Step 3.
                degree_res = residue.degree()
                pos_at_inf_of_res = residue.leading_coefficient() > 0

                if pos_at_inf_of_res != pos_at_inf:
                    pos_at_inf = not pos_at_inf
                    r_1 -= 1

                if pos_at_inf_of_res != (pos_at_neg == (degree_res & 1 == 0)):
                    pos_at_neg = not pos_at_neg
                    r_1 += 1

                #Step 4.
                if degree_res == 0:
                    return (r_1, (degree - r_1) >> 1)

                A, B = B, residue.scalar_exact_division(g*(h**deg))
                g = abs(A.leading_coefficient())
                if deg == 1:
                    h = g
                elif deg > 1:
                    h = g**deg // h**(deg - 1)
        return self.sign

    def POLRED(self):
        """
        Given a polynomial f i.e. a field self, output some polynomials
        defining subfield of self, where self is a field defined by f.
        Algorithm 4.4.11 in Cohen's book.
        """
        n = self.degree
        appr = self.getConj()

        #Step 1.
        Basis = self.integer_ring()
        BaseList = []
        for i in range(n):
            AlgInt = MatAlgNumber(Basis[i].compo, self.polynomial)
            BaseList.append(AlgInt)

        #Step 2.
        traces = []
        if self.signature()[1] == 0:
            for i in range(n):
                for j in range(n):
                    s = BaseList[i]*BaseList[j]
                    traces.append(s.trace())
        else:
            sigma = self.getConj()
            f = []
            for i in range(n):
                f.append(zpoly(Basis[i]))
            for i in range(n):
                for j in range(n):
                    m = 0 + 0j
                    for k in range(n):
                        m += f[i](sigma[k])*f[j](sigma[k].conjugate())
                    traces.append(m.real)

        #Step 3.
        M = matrix.createMatrix(n, n, traces)
        S = matrix.unitMatrix(n)
        L = lattice.LLL(S, M)[0]


        #Step 4.
        Ch_Basis = []
        for i in range(n):
            base_cor = changetype(0, self.polynomial).ch_matrix()
            for v in range(n):
                base_cor += BaseList[v]*L.compo[v][i]
            Ch_Basis.append(base_cor)

        C = []
        a = Ch_Basis[0]
        for i in range(n):
            coeff = Ch_Basis[i].ch_basic().getCharPoly()
            C.append(zpoly(coeff))

        #Step 5.
        P = []
        for i in range(n):
            diff_C = C[i].differentiate()
            gcd_C = C[i].subresultant_gcd(diff_C)
            P.append(C[i].exact_division(gcd_C))

        return P

    def isIntBasis(self):
        """
        Determine whether self.basis is integral basis of self field.
        """
        D = self.disc()
        if squarefree.trial_division(abs(D)):
            return True
        else:
            if D & 3 == 0:
                if squarefree.trial_division(abs(D) >> 2) and (D >> 2) & 3 != 1:
                    return True
        # compare with real integer ring
        return abs(self.integer_ring().determinant()) == 1

    def isGaloisField(self, other = None):
        """
        Determine whether self/other is Galois field.
        """
        if self.signature[0] == 0 or self.signature[1] == 0:
            return "Can not determined"
        else:
            return False

    def isFieldElement(self, A):
        """
        Determine whether A is field element of self field or not.
        """
        poly = A.polynomial
        if poly == self.polynomial:
            return True
        else:
            if poly == self.POLRED():
                return True
            else:
                return False

    def getCharacteristic(self):
        """
        Return characteristic of the field (it is always zero).
        """
        return 0

    def createElement(self, *seed):
        """
        createElement returns an element of the field with seed.
        """
        if len(seed) == 1:
            if isinstance(seed[0], list):
                if isinstance(seed[0][0], list):
                    return BasicAlgNumber(seed[0], self.polynomial)
                else:
                    return MatAlgNumber(seed[0], self.polynomial)
            else:
                raise NotImplementedError
        else:
            raise NotImplementedError

    def issubring(self, other):
        """
        Report whether another ring contains the field as a subring.
        """
        if self is other or self.isSubField(other):
            return True
        raise NotImplementedError("don't know how to tell")

    def issuperring(self, other):
        """
        Report whether the field is a superring of another ring.
        """
        if self is other:
            return True
        elif other.issubring(rational.theRationalField):
            return True
        raise NotImplementedError("don't know how to tell")

    def __eq__(self, other):
        """
        Equality test.
        """
        if self.issubring(other) and self.issuperring(other):
            return True
        return False

    def __hash__(self):
        raise NotImplementedError

    # properties
    def _getOne(self):
        "getter for one"
        if self._one is None:
            self._one = BasicAlgNumber([[1] + [0] * (self.degree - 1), 1],
                              self.polynomial)
        return self._one

    one = property(_getOne, None, None, "multiplicative unit.")

    def _getZero(self):
        "getter for zero"
        if self._zero is None:
            self._zero = BasicAlgNumber([[0] + [0] * (self.degree - 1), 1],
                              self.polynomial)
        return self._zero

    zero = property(_getZero, None, None, "additive unit.")


class BasicAlgNumber(object):
    """
    The class for algebraic number.
    """
    def __init__(self, valuelist, polynomial, precompute=False):
        if len(polynomial) != len(valuelist[0])+1:
            raise ValueError
        self.value = valuelist
        self.coeff = valuelist[0]
        self.denom = valuelist[1]
        self.degree = len(polynomial) - 1
        self.polynomial = polynomial
        self.field = NumberField(self.polynomial)
        Gcd = gcd.gcd_of_list(self.coeff)
        GCD = gcd.gcd(Gcd[0], self.denom)
        if GCD != 1:
            self.coeff = [i//GCD for i in self.coeff]
            self.denom = self.denom//GCD
        if precompute:
            self.getConj()
            self.getApprox()
            self.getCharPoly()

    def __repr__(self):
        return_str = '%s(%s, %s)' % (self.__class__.__name__, [self.coeff, self.denom], self.polynomial)
        return return_str

    def __neg__(self):
        coeff = []
        for i in range(len(self.coeff)):
            coeff.append(-self.coeff[i])
        return BasicAlgNumber([coeff, self.denom], self.polynomial)

    def _int_to_algnumber(self, other):
        """
        other (integer) -> BasicAlgNumber (over self.field)
        """
        return BasicAlgNumber([[other] + [0] * (self.degree - 1), 1],
                              self.polynomial)

    def _rational_to_algnumber(self, other):
        """
        other (rational) -> BasicAlgNumber (over self.field)
        """
        return BasicAlgNumber([[other.numerator] + [0] * (self.degree - 1),
                               other.denominator], self.polynomial)

    def __add__(self, other):
        if isinstance(other, BasicAlgNumber):
            d = self.denom*other.denom
            coeff = []
            for i in range(len(self.coeff)):
                coeff.append(
                    other.denom*self.coeff[i] + self.denom*other.coeff[i])
            return BasicAlgNumber([coeff, d], self.polynomial)
        elif isinstance(other, (int, long)):
            return self + self._int_to_algnumber(other)
        elif isinstance(other, rational.Rational):
            return self + self._rational_to_algnumber(other)
        else:
            return NotImplemented

    __radd__ = __add__

    def __sub__(self, other):
        return self.__add__(-other)

    def __mul__(self, other):
        if isinstance(other, BasicAlgNumber):
            d = self.denom*other.denom
            f = zpoly(self.polynomial)
            g = zpoly(self.coeff)
            h = zpoly(other.coeff)
            j = (g * h).pseudo_mod(f)
            jcoeff = [j[i] for i in range(self.degree)]
            return BasicAlgNumber([jcoeff, d], self.polynomial)
        elif isinstance(other, (int, long)):
            Coeff = [i*other for i in self.coeff]
            return BasicAlgNumber([Coeff, self.denom], self.polynomial)
        elif isinstance(other, rational.Rational):
            GCD = gcd.gcd(other.numerator, self.denom)
            denom = self.denom * other.denominator
            mul = other.numerator
            if GCD != 1:
                denom //= GCD
                mul //= GCD
            Coeff = [self.coeff[j] * mul for j in range(self.degree)]
            return BasicAlgNumber([Coeff, denom], self.polynomial)
        else:
            return NotImplemented

    __rmul__ = __mul__

    def __pow__(self, exponent, mod=None):
        d = self.denom**exponent
        f = zpoly(self.polynomial)
        g = zpoly(self.coeff)
        if mod is None:
            j = pow(g, exponent, f)
        else:
            j = pow(g, exponent, f) % mod
            # what does mod exactly means?
        jcoeff = [j[i] for i in range(self.degree)]
        return BasicAlgNumber([jcoeff, d], self.polynomial)

    def inverse(self):
        f = qpoly(self.polynomial)
        g = qpoly(self.coeff)
        quot = f.extgcd(g)

        icoeff = [self.denom * quot[1][i] for i in range(self.degree)]
        #print icoeff
        denom_list = []
        for i in range(self.degree):
            if icoeff[i] == 0:
                denom_list.append(1)
            else:
                denom_list.append(icoeff[i].denominator)
        new_denom = 1
        for i in range(self.degree):
            new_denom = gcd.lcm(new_denom, denom_list[i])
        icoeff = [icoeff[i].numerator * new_denom // icoeff[i].denominator for i in range(self.degree)]
        return BasicAlgNumber([icoeff, new_denom], self.polynomial)

    def __truediv__(self, other):
        f = zpoly(self.polynomial)
        g = zpoly(self.coeff)
        t = BasicAlgNumber([other.coeff, other.denom], self.polynomial)
        k = t.inverse()
        h = zpoly(k.coeff)
        d = self.denom * k.denom
        j = (g * h).monic_mod(f)
        jcoeff = [j[i] for i in range(self.degree)]
        return BasicAlgNumber([jcoeff, d], self.polynomial)

    __div__ = __truediv__
    __floordiv__ = __truediv__

    def getConj(self):
        """
        Return (approximate) solutions of self.polynomial.
        We can discriminate the conjugate field of self by these values.
        """
        if not hasattr(self, "conj"):
            self.conj = NumberField(self.polynomial).getConj()
        return self.conj

    def getApprox(self):
        """
        Return the approximations of all conjugates of self.
        """
        if not hasattr(self, "approx"):
            APP = []
            for i in range(self.degree):
                Approx = 0
                for j in range(self.degree):
                    Approx += self.coeff[j]*(self.getConj()[i]**j)
                APP.append(Approx)
            self.approx = APP
        return self.approx

    def getCharPoly(self):
        """
        Return the characteristic polynomial of self
        by compute products of (x-self^{(i)}).
        """
        if not hasattr(self, "charpoly"):
            Conj = self.getApprox()
            P = uniutil.polynomial({0:-Conj[0], 1:1}, ring.getRing(Conj[0]))
            for i in range(1, self.degree):
                P *= uniutil.polynomial({0:-Conj[i], 1:1},
                    ring.getRing(Conj[i]))
            charcoeff = []
            for i in range(self.degree + 1):
                if hasattr(P[i], "real"):
                    charcoeff.append(int(math.floor(P[i].real + 0.5)))
                else:
                    charcoeff.append(int(math.floor(P[i] + 0.5)))
            self.charpoly = charcoeff
        return self.charpoly

    def getRing(self):
        """
        Return the algebraic number field contained self.
        """
        return NumberField(self.polynomial)

    def trace(self):
        """
        Compute the trace of self in K.
        """
        denom = self.denom
        n = len(self.polynomial) - 1

        tlist = [n]
        for k in range(1, n):
            s = 0
            for i in range(1, k):
                s += tlist[k - i] * self.polynomial[n - i]
            tlist.append(-k * self.polynomial[n - k] - s)

        t = 0
        for j in range(len(tlist)):
            t += tlist[j] * self.coeff[j]

        if denom == 1:
            return t
        elif t % denom == 0:
            return t // denom
        else:
            return rational.Rational(t, denom)

    def norm(self):
        """
        Compute the norm of self in K.
        """
        f = zpoly(self.polynomial)
        g = zpoly(self.coeff)
        R = f.resultant(g)

        if self.denom == 1:
            return int(R)
        else:
            denom = self.denom**self.degree
            if isinstance(R, int):
                if R % denom == 0:
                    return R // denom
                else:
                    return rational.Rational(R, denom)
            else:
                return R / denom

    def isAlgInteger(self):
        """
        Determine whether self is an algebraic integer or not.
        """
        Norm = self.norm()
        if isinstance(Norm, int):
            return True
        else:
            return False

    def ch_matrix(self):
        """
        Change style to MatAlgNumber.
        """
        list = []
        if self.denom == 1:
            list = self.coeff
        else:
            for i in range(self.degree):
                list.append(rational.Rational(self.coeff[i], self.denom))
        return MatAlgNumber(list, self.polynomial)

class MatAlgNumber(object):
    """
    The class for algebraic number represented by matrix.
    """
    def __init__(self, coefficient, polynomial):
        """
        """
        self.coeff = coefficient
        self.degree = len(coefficient)
        List = []
        for i in range(self.degree):
            stbasis = [0] * self.degree
            stbasis[i] = 1
            List.append(stbasis)
        List.append([- polynomial[i] for i in range(self.degree)])
        for l in range(self.degree - 2):
            basis1 = []
            basis = []
            for j in range(self.degree):
                if j == 0:
                    basis1.append(0)
                if j < self.degree - 1:
                    basis1.append(List[l + self.degree][j])
                elif j == self.degree - 1:
                    basis2 = [List[l + self.degree][j] * - polynomial[k] for k in range(self.degree)]
            for i in range(self.degree):
                basis.append(basis1[i] + basis2[i])
            List.append(basis)
        Matrix = []
        flag = 0
        for i in range(self.degree):
            basis3 = []
            for j in range(self.degree):
                basis3.append([self.coeff[j] * k for k in List[j+flag]])
            for l in range(self.degree):
                t = 0
                for m in range(self.degree):
                    t += basis3[m][l]
                Matrix.append(t)
            flag += 1

        self.matrix = matrix.createMatrix(self.degree, self.degree, Matrix)
        self.polynomial = polynomial
        self.field = NumberField(self.polynomial)

    def __repr__(self):
        return_str = '%s(%s, %s)' % (self.__class__.__name__, self.matrix.__repr__(), self.polynomial)
        return return_str

    def __neg__(self):
        mat = - self.matrix
        coeff = []
        for i in range(mat.row):
            coeff.append(mat[0][i])
        return MatAlgNumber(coeff, self.polynomial)

    def __add__(self, other):
        if isinstance(other, (int, long, rational.Rational)):
            other = MatAlgNumber(
                [other] + [0] * (self.degree -1), self.polynomial)
        elif not isinstance(other, MatAlgNumber):
            return NotImplemented
        mat = self.matrix + other.matrix
        coeff = []
        for i in range(mat.row):
            coeff.append(mat[i+1][1])
        return MatAlgNumber(coeff, self.polynomial)

    __radd__ = __add__

    def __sub__(self, other):
        if isinstance(other, (int, long, rational.Rational)):
            other = MatAlgNumber(
                [other] + [0] * (self.degree -1), self.polynomial)
        elif not isinstance(other, MatAlgNumber):
            return NotImplemented
        mat = self.matrix - other.matrix
        coeff = []
        for i in range(mat.row):
            coeff.append(mat[i+1][1])
        return MatAlgNumber(coeff, self.polynomial)

    def __mul__(self, other):
        if isinstance(other, MatAlgNumber):
            mat = self.matrix * other.matrix
        elif isinstance(other, (int, long, rational.Rational)):
            mat = other * self.matrix
        else:
            return NotImplemented
        coeff = []
        for i in range(mat.row):
            coeff.append(mat[i+1][1])
        return MatAlgNumber(coeff, self.polynomial)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if isinstance(other, (int, long, rational.Rational)):
            return self * (rational.Rational(other) ** -1)
        elif not isinstance(other, MatAlgNumber):
            return NotImplemented
        other_mat = other.matrix.copy()
        other_mat.toFieldMatrix()
        mat = other_mat.inverse(self.matrix)
        coeff = []
        for i in range(mat.row):
            coeff.append(mat[i+1][1])
        return MatAlgNumber(coeff, self.polynomial)

    __div__ = __truediv__
    __floordiv__ = __truediv__

    def __pow__(self, other):
        mat = self.matrix ** other
        coeff = []
        for i in range(mat.row):
            coeff.append(mat[i+1][1])
        return MatAlgNumber(coeff, self.polynomial)

    def inverse(self):
        (self.matrix).toFieldMatrix()
        inv = (self.matrix).inverse()
        coeff = []
        for i in range(inv.row):
            coeff.append(inv[i+1][1])
        return MatAlgNumber(coeff, self.polynomial)

    def norm(self):
        return (self.matrix).determinant()

    def trace(self):
        return (self.matrix).trace()

    def getRing(self):
        """
        Return the algebraic number field contained self.
        """
        return NumberField(self.polynomial)

    def ch_basic(self):
        denom = 1
        for i in range(self.degree):
            if not isinstance(self.coeff[i], int):
                denom *= gcd.lcm(denom, (self.coeff[i]).denominator)
        coeff = []
        for i in range(self.degree):
            if isinstance(self.coeff[i], int):
                coeff.append(self.coeff[i] * denom)
            else:
                coeff.append(int((self.coeff[i]).numerator * denom / (self.coeff[i]).denominator))
        return BasicAlgNumber([coeff, denom], self.polynomial)

def prime_decomp(p, polynomial):
    """
    Return prime decomposition of (p) over Q[x]/(polynomial).
    """
    field = NumberField(polynomial)
    field_disc = field.field_discriminant()
    if (field.disc() // field_disc) % p != 0:
        return _easy_prime_decomp(p, polynomial)
    else:
        return _main_prime_decomp(p, polynomial)

def _easy_prime_decomp(p, polynomial):
    """
    prime decomposition by factoring polynomial
    """
    # import nzmath.module as module
    f = fppoly(polynomial, p)
    d = f.degree()
    factor = f.factor()
    p_alg = BasicAlgNumber([[p] + [0] * (d - 1), 1], polynomial)
    if len(factor) == 1 and factor[0][1] == 1:
        return [(p_alg , 1 )]
    else:
        dlist = []
        for i in range(len(factor)):
            basis_list = []
            for j in range(d):
                if factor[i][0][j] == 0:
                    basis_list.append(0)
                else:
                    basis_list.append(factor[i][0][j].toInteger())
            dlist.append([BasicAlgNumber([basis_list, 1], polynomial),
                          factor[i][1]])
        return [( (p_alg, dlist_ele[0]), dlist_ele[1] ) for dlist_ele in dlist]

def _main_prime_decomp(p, polynomial):
    """
    main step of prime decomposition
    """
    # import nzmath.module as module
    raise NotImplementedError

def changetype(a, polynomial=[0, 1]):
    """
    Change 'a' to be an element of field K defined polynomial
    """
    if isinstance(a, (int, long)):
        coeff = [a] + [0] * (len(polynomial) - 2)
        return BasicAlgNumber([coeff, 1], polynomial)
    elif isinstance(a, rational.Rational):
        coeff = [a.numerator] + [0] * (len(polynomial) - 2)
        return BasicAlgNumber([coeff, a.denominator], polynomial)
    else:
        deg = len(a) - 1
        coeff = [0, 1] + [0] * (deg - 2)
        if a[-1] != 1:
            polynomial = [a[-2], 1]
            mul = a[-1]
            for j in range(deg - 2, -1, -1):
                polynomial.insert(0, a[j] * mul)
                mul *= a[-1]
            return BasicAlgNumber([coeff, a[-1]], polynomial)
        else:
            return BasicAlgNumber([coeff, 1], a)

def disc(A):
    """
    Compute the discriminant of a_i, where A=[a_1,...,a_n]
    """
    n = A[0].degree
    list = []
    for i in range(n):
        for j in range(n):
            s = A[i] * A[j]
            list.append(s.trace())
    M = matrix.createMatrix(n, n, list)
    return M.determinant()

def qpoly(coeffs):
    """
    Return a rational coefficient polynomial constructed from given
    coeffs.  The coeffs is a list of coefficients in ascending order.
    """
    terms = [(i, rational.Rational(c)) for (i, c) in enumerate(coeffs)]
    return uniutil.polynomial(terms, rational.theRationalField)

def zpoly(coeffs):
    """
    Return an integer coefficient polynomial constructed from given
    coeffs.  The coeffs is a list of coefficients in ascending order.
    """
    return uniutil.polynomial(enumerate(coeffs), rational.theIntegerRing)

def fppoly(coeffs, p):
    """
    Return a Z_p coefficient polynomial constructed from given
    coeffs.  The coeffs is a list of coefficients in ascending order.
    """
    return uniutil.polynomial(enumerate(coeffs), finitefield.FinitePrimeField(p))
