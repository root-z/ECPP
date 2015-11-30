"""
module, ideal etc. for number field
"""
from __future__ import division

from nzmath.algfield import BasicAlgNumber, MatAlgNumber
import nzmath.arith1 as arith1
import nzmath.gcd as gcd
import nzmath.vector as vector
import nzmath.matrix as matrix
import nzmath.rational as rational
import nzmath.ring as ring
import nzmath.prime as prime


class Submodule(matrix.RingMatrix):
    """
    Submodule is a class for submodules (typically, Z^n).
    we assume that coeff_ring is PID.
    (i.e. hermite normal form (HNF) for matrices exists.)
    """
    
    def __init__(self, row, column, compo=0, coeff_ring=0, ishnf=None):
        """
        Submodule(row, column [,components, coeff_ring, ishnf])
        """
        if isinstance(compo, bool):
            ishnf = compo
            compo = 0
        elif isinstance(coeff_ring, bool):
            ishnf = coeff_ring
            coeff_ring = 0
        self._initialize(row, column, compo, coeff_ring)
        self.ishnf = ishnf

    @classmethod
    def fromMatrix(cls, mat, ishnf=None):
        """
        A constructor class method, which creates Submodule from a
        Matrix instance.
        """
        compo = []
        for row in mat.compo:
            compo += row
        return cls(mat.row, mat.column, compo, mat.coeff_ring, ishnf)

    def getGenerators(self):
        """
        return (current) generator of self (i.e. essentially, self.compo)
        """
        return [self[j] for j in range(1, self.column+1)]

    def isSubmodule(self, other):
        """ 
        Check self is in other as submodule
        """
        try:
            return self.sumOfSubmodules(other).isEqual(other)
        except:
            return False

    def isEqual(self, other):
        """
        determine self == other as module
        """
        if self.ishnf:
            if other.ishnf:
                return self == other
            else:
                other_copy = Submodule.fromMatrix(other)
                other_copy.toHNF()
                return self == other_copy
        else:
            if other.ishnf:
                self_copy = Submodule.fromMatrix(self)
                self_copy.toHNF()
                return self_copy == other
            else:
                self_copy = Submodule.fromMatrix(self)
                self_copy.toHNF()
                other_copy = Submodule.fromMatrix(other)
                other_copy.toHNF()
                return self_copy == other_copy

    def isContains(self, other):
        """
        determine whether other is in self or not.
        if you want other to be represented as self, use represent_element.
        """
        self_copy = Submodule.fromMatrix(self)
        return not isinstance(self_copy.represent_element(other), bool)

    def toHNF(self):
        """
        Change matrix to HNF(hermite normal form)
        Note that HNF do not always give basis of self
        (i.e. HNF may be redundant)
        """
        if not self.ishnf:
            hnf = self.hermiteNormalForm(True)
            if hnf == None: #zero module
                hnf = matrix.zeroMatrix(self.row, 1, self.coeff_ring)
            self.compo = hnf.compo
            self.column = hnf.column
            self.ishnf = True

    def sumOfSubmodules(self, other):
        """
        Return module which is sum of self and other.
        """
        if self.row != other.row:
            raise matrix.MatrixSizeError()
        N = self.copy()
        N.extendColumn(other)
        module = Submodule.fromMatrix(N)
        module.toHNF()
        return module

    def intersectionOfSubmodules(self, other):
        """
        Return module which is intersection of self and other.
        """
        if self.row != other.row:
            raise matrix.MatrixSizeError()
        M1 = self.copy()
        M1.extendColumn(other)
        N = M1.kernelAsModule()
        if N == None: # zero kernel
            N = matrix.zeroMatrix(self.row, 1, self.coeff_ring)
        N1 = N.getBlock(1, 1, self.column, N.column) # N.column is dim(ker(M1))
        module = Submodule.fromMatrix(self * N1)
        module.toHNF()
        return module

    def represent_element(self, other):
        """
        represent other as a linear combination with hnf generators
        if other not in self, return False
        this method execute self.toHNF()
        """
        self.toHNF()
        return self._HNF_solve(other)

    def linear_combination(self, coeff):
        """
        return a vector corresponding a linear combination of (current) basis
        with given Z-coefficients.
        """
        sol = coeff[0] * self[1]
        for i in range(2, self.column + 1):
            sol += coeff[i - 1] * self[i]
        return sol

    def _HNF_solve(self, other):
        """
        return X over coeff_ring s.t. self * X = other,
        where self is HNF, other is vector
        """
        if self.row != len(other):
            raise vector.VectorSizeError()
        n = self.column
        zero = self.coeff_ring.zero
        X = vector.Vector([zero] * n)
        # search non-zero entry
        for i in range(1, self.row + 1)[::-1]:
            if self[i, n] != zero:
                break
            if other[i] != zero:
                return False
        else:
            return X
        # solve HNF system (only integral solution)
        for j in range(1, n + 1)[::-1]:
            sum_i = other[i]
            for k in range(j + 1, n + 1):
                sum_i -= X[k] * self[i, k]
            try:
                X[j] = ring.exact_division(sum_i, self[i, j])
            except ValueError: # division is not exact
                return False
            ii = i
            i -= 1
            for i in range(1, ii)[::-1]:
                if j != 1 and self[i, j - 1] != zero:
                    break
                sum_i = X[j] * self[i, j]
                for k in range(j + 1, n + 1):
                    sum_i += X[k] * self[i, k]
                if sum_i != other[i]:
                    return False
            if i == 0:
                break
        return X

class Module(object):
    """
    for computing module with HNF over a number field.
    A module is a finitely generated sub Z-module.
    (Note that we do not assume rank of a module is deg(number_field).)
    """
    def __init__(self, pair_mat_repr, number_field, base=None, ishnf=False):
        """
        pair_mat_repr: generators represented with respect to base_module
                       1: [[deg(number_field))-integral tuple/vector], denominator]
                       2: [integral matrix whose the number of row is deg(number_field), 
                              denominator]
                       3: rational matrix whose the number of row is deg(number_field)
        number_field:  an instance of NumberField
        base:          a base module represented with respect to Z[theta]
                       1: [deg(number_field)-rational tuple/vector], which represents basis
                       2: deg(number_field)-square non-singular rational matrix
        ishnf:         if hnf_repr is already HNF-form, set True
        """
        self.number_field = number_field
        size = self.number_field.degree
        if isinstance(base, bool):
            ishnf = base
            base = None
        if base == None: # power basis
            self.base = matrix.unitMatrix(size, rational.theIntegerRing)
            self.base.toFieldMatrix()
        elif isinstance(base, list): # list repr
            self.base = matrix.FieldSquareMatrix(size, base)
        else: # matrix repr
            if base.row != size:
                raise TypeError("The number of row of base is incorrect.")
            self.base = base.copy()
            self.base.toFieldMatrix()
        if not isinstance(pair_mat_repr, list): #rational matrix repr
            self.mat_repr, self.denominator = _toIntegerMatrix(
                pair_mat_repr)
        else:
            self.denominator = pair_mat_repr[1]
            if isinstance(pair_mat_repr[0], list): # list repr
                column = len(pair_mat_repr[0])
                self.mat_repr = Submodule(
                    size, column, pair_mat_repr[0], ishnf)
            else: # integral matrix repr
                self.mat_repr = Submodule.fromMatrix(
                    pair_mat_repr[0].copy(), ishnf)

    def _simplify(self):
        """
        simplify self.denominator
        """
        mat_repr, numer = _toIntegerMatrix(self.mat_repr, 2)
        denom_gcd = gcd.gcd(numer, self.denominator)
        if denom_gcd != 1:
            self.denominator = ring.exact_division(self.denominator, denom_gcd)
            self.mat_repr = Submodule.fromMatrix(
                ring.exact_division(numer, denom_gcd) * mat_repr)

    def toHNF(self):
        """
        transform mat_repr to hnf form
        """
        self.mat_repr.toHNF()

    def __repr__(self):
        """
        formal representation
        """
        return_str = '%s([%s, %s], %s, %s)' % ( self.__class__.__name__, 
          repr(self.mat_repr), repr(self.denominator), 
          repr(self.number_field), repr(self.base) )
        return return_str

    def __str__(self):
        """
        simple representation
        """
        return str(
            (self.mat_repr, self.denominator)) + "\n over \n" + str(
            (self.base, self.number_field))

    def __eq__(self, other):
        """
        Check whether self == other or not (as module)
        """
        #if not(other in self.number_field):
        #    return False
        changed_self = self.change_base_module(other.base)
        return (changed_self.denominator == other.denominator
            ) and changed_self.mat_repr.isEqual(other.mat_repr)

    def __ne__(self, other):
        return self != other

    def __contains__(self, other):
        """
        Check whether other is an element of self or not
        """
        #if not(other in self.number_field):
        #    return False
        if isinstance(other, (tuple, vector.Vector)): #vector repr
            other_int, other_denom = _toIntegerMatrix(vector.Vector(other))
        elif isinstance(other, list):
            if isinstance(other[0], (tuple, vector.Vector)):
                other_int = vector.Vector(other[0])
                other_denom = other[1]
            else:
                raise ValueError("input [vector, denominator]!")
        else: # field element
            if not isinstance(other, BasicAlgNumber):
                other_copy = other.ch_basic()
            else:
                other_copy = other
            other_base_repr = self.base.inverse(
                vector.Vector(other_copy.coeff))
            other_int, other_int_denom = _toIntegerMatrix(other_base_repr)
            other_denom = other_int_denom * other_copy.denom
        if isinstance(other_int, matrix.Matrix):
            other_int = other_int[1]
        try:
            numer = ring.exact_division(self.denominator, other_denom)
        except ValueError: # division is not exact
            return False
        if numer != 1:
            other_vect = numer * other_int
        else:
            other_vect = other_int
        return not isinstance(
            self.mat_repr.represent_element(other_vect), bool)

    def __add__(self, other):
        """
        return self + other (as module)
        """
        # if self.number_field != other.number_field:
        #     raise ValueError("based on different number field")
        new_self = self.change_base_module(other.base)
        new_denominator = gcd.lcm(new_self.denominator, other.denominator)
        if new_self.denominator != new_denominator:
            multiply_denom = ring.exact_division(new_denominator,
                new_self.denominator)
            new_self_mat_repr = multiply_denom * new_self.mat_repr
        else:
            new_self_mat_repr = new_self.mat_repr
        if other.denominator != new_denominator:
            multiply_denom = ring.exact_division(new_denominator,
                other.denominator)
            new_other_mat_repr = multiply_denom * other.mat_repr
        else:
            new_other_mat_repr = other.mat_repr
        new_mat_repr = Submodule.fromMatrix(
            new_self_mat_repr).sumOfSubmodules(new_other_mat_repr)
        new_module = self.__class__(
            [new_mat_repr, new_denominator], self.number_field, other.base)
        new_module._simplify()
        return new_module

    def __mul__(self, other):
        """
        return self * other (as ideal computation)
        """
        if isinstance(other, 
            (int, long, rational.Integer, rational.Rational)):
            return self._rational_mul(other)
        #try:
        #     use __contains__ function?
        #    if other in self.number_field:
        if isinstance(other, (BasicAlgNumber, MatAlgNumber)):
            return self._scalar_mul(other)
        #except: # other is a module
        return self._module_mul(other)

    __rmul__ = __mul__

    def _rational_mul(self, other):
        """
        return other * self, assuming other is a rational element
        """
        if rational.isIntegerObject(other):
            other_numerator = rational.Integer(other)
            other_denominator = rational.Integer(1)
        else:
            other_numerator = other.numerator
            other_denominator = other.denominator
        denom_gcd = gcd.gcd(self.denominator, other_numerator)
        if denom_gcd != 1:
            new_denominator = ring.exact_division(
                self.denominator, denom_gcd) * other_denominator
            multiply_num = other_numerator.exact_division(denom_gcd)
        else:
            new_denominator = self.denominator * other_denominator
            multiply_num = other_numerator
        new_module = self.__class__(
                         [self.mat_repr * multiply_num, new_denominator],
                         self.number_field, self.base, self.mat_repr.ishnf)
        new_module._simplify()
        return new_module

    def _scalar_mul(self, other):
        """
        return other * self, assuming other is a scalar
        (i.e. an element of self.number_field)
        """
        # use other is an element of higher or lower degree number field ?
        #if other.getRing() != self.number_field:
        #    raise NotImplementedError(
        #        "other is not a element of number field")
        if not isinstance(other, BasicAlgNumber):
            try:
                other = other.ch_basic()
            except:
                raise NotImplementedError(
                    "other is not an element of a number field")
        # represent other with respect to self.base
        other_repr, pseudo_other_denom = _toIntegerMatrix(
            self.base.inverse(vector.Vector(other.coeff)))
        other_repr = other_repr[1]
        # compute mul using self.base's multiplication
        base_mul = self._base_multiplication()
        n = self.number_field.degree
        mat_repr = []
        for k in range(1, self.mat_repr.column +1):
            mat_repr_ele = vector.Vector([0] * n)
            for i1 in range(1, n + 1):
                for i2 in range(1, n + 1):
                    mat_repr_ele += self.mat_repr[i1, k] * other_repr[
                        i2] * _symmetric_element(i1, i2, base_mul)
            mat_repr.append(mat_repr_ele)
        mat_repr, denom, numer = _toIntegerMatrix(
            matrix.FieldMatrix(n, len(mat_repr), mat_repr), 1)
        denom = self.denominator * pseudo_other_denom * other.denom * denom
        gcds = gcd.gcd(denom, numer)
        if gcds != 1:
            denom = ring.exact_division(denom, gcds)
            numer = ring.exact_division(numer, gcds)
        if numer != 1:
            mat_repr = numer * mat_repr
        else:
            mat_repr = mat_repr
        return self.__class__([mat_repr, denom], self.number_field, self.base)

    def _module_mul(self, other):
        """
        return self * other as the multiplication of modules
        """
        #if self.number_field != other.number_field:
        #    raise NotImplementedError
        if self.base != other.base:
            new_self = self.change_base_module(other.base)
        else:
            new_self = self.copy()
        base_mul = other._base_multiplication()
        n = self.number_field.degree
        new_mat_repr = []
        for k1 in range(1, new_self.mat_repr.column + 1):
            for k2 in range(1, other.mat_repr.column + 1):
                new_mat_repr_ele = vector.Vector([0] * n)
                for i1 in range(1, n + 1):
                    for i2 in range(1, n + 1):
                        new_mat_repr_ele += new_self.mat_repr[
                            i1, k1] * other.mat_repr[i2, k2
                            ] * _symmetric_element(i1, i2, base_mul)
                new_mat_repr.append(new_mat_repr_ele)
        new_mat_repr, denom, numer = _toIntegerMatrix(
            matrix.FieldMatrix(n, len(new_mat_repr), new_mat_repr), 1)
        denom = new_self.denominator * other.denominator * denom
        gcds = gcd.gcd(denom, numer)
        if gcds != 1:
            denom = ring.exact_division(denom, gcds)
            numer = ring.exact_division(numer, gcds)
        if numer != 1:
            mat_repr = numer * new_mat_repr
        else:
            mat_repr = new_mat_repr
        sol = self.__class__(
            [mat_repr, denom], self.number_field, other.base)
        sol.toHNF()
        return sol

    def __pow__(self, other):
        """
        self ** other (based on ideal multiplication)
        """
        if other <= 0:
            raise ValueError, "only support non-negative powering" 
        mul_part = self.copy()
        index = other
        while True:
            if index & 1:
                try:
                    sol *= mul_part
                except NameError:
                    sol = mul_part
            index >>= 1
            if not index:
                return sol
            else:
                mul_part *= mul_part    

    def _base_multiplication(self):
        """
        return [base[i] * base[j]] (as a numberfield element)
        this is a precomputation for computing _module_mul
        """
        if not hasattr(self, "_base_multiply"):
            base_mul = _base_multiplication(self.base, self.number_field)
            self._base_multiply = base_mul
        return self._base_multiply

    def copy(self):
        """
        create copy of self
        """
        new_pair_mat_repr = [[self.mat_repr.copy()[i] for i in range(1, 
            self.mat_repr.column + 1)], self.denominator]
        new_base_module = [self.base.copy()[i] for i in range(
            1, self.base.column + 1)]
        return self.__class__(new_pair_mat_repr, self.number_field, 
            new_base_module, True)

    def intersect(self, other):
        """
        return intersection of self and other
        """
        # if self.number_field != other.number_field:
        #     raise ValueError("based on different number field")
        if self.base != other.base:
            new_self = self.change_base_module(other.base)
        else:
            new_self = self.copy()
        new_denominator = gcd.lcm(new_self.denominator, other.denominator)
        if new_self.denominator != new_denominator:
            multiply_denom = ring.exact_division(
                new_denominator, new_self.denominator)
            new_self_mat_repr = multiply_denom * new_self.mat_repr
        else:
            new_self_mat_repr = new_self.mat_repr
        if other.denominator != new_denominator:
            multiply_denom = ring.exact_division(
                new_denominator, other.denominator)
            new_other_mat_repr = multiply_denom * other.mat_repr
        else:
            new_other_mat_repr = other.mat_repr
        new_mat_repr = Submodule.fromMatrix(
            new_self_mat_repr).intersectionOfSubmodules(new_other_mat_repr)
        new_module = self.__class__(
            [new_mat_repr, new_denominator], self.number_field, other.base)
        new_module._simplify()
        return new_module
    
    def issubmodule(self, other):
        """
        Check self is submodule of other.
        """
        if self.__class__.__name__ != other.__class__.__name__:
            return False
        if (self + other) == other: # as module
            return True
        else:
            return False

    def issupermodule(self, other):
        """
        Check other is supermodule of self.
        """
        if self.__class__.__name__ != other.__class__.__name__:
            return False
        else:
            return other.issubmodule(self)

    def represent_element(self, other):
        """
        represent other as a linear combination with generators of self
        if other not in self, return False
        Note that we do not assume self.mat_repr is HNF
        """
        #if other not in self.number_field:
        #    return False
        theta_repr = (self.base * self.mat_repr)
        theta_repr.toFieldMatrix()
        pseudo_vect_repr = theta_repr.inverseImage(
            vector.Vector(other.coeff))
        pseudo_vect_repr = pseudo_vect_repr[1]
        gcd_self_other = gcd.gcd(self.denominator, other.denom)
        multiply_numerator = self.denominator // gcd_self_other
        multiply_denominator = other.denom // gcd_self_other

        def getNumerator_and_Denominator(ele):
            """
            get the pair of numerator and denominator of ele
            """
            if rational.isIntegerObject(ele):
                return (ele, 1)
            else: # rational
                return (ele.numerator, ele.denominator)

        list_repr = []
        for j in range(1, len(pseudo_vect_repr) + 1):
            try:
                numer, denom = getNumerator_and_Denominator(
                    pseudo_vect_repr[j])
                list_repr_ele = ring.exact_division(
                    numer * multiply_numerator, denom * multiply_denominator)
                list_repr.append(list_repr_ele)
            except ValueError: # division is not exact
                return False
        return list_repr

    def change_base_module(self, other_base):
        """
        change base_module to other_base_module
        (i.e. represent with respect to base_module)
        """
        if self.base == other_base:
            return self
        n = self.number_field.degree
        if isinstance(other_base, list):
            other_base = matrix.FieldSquareMatrix(n,
                [vector.Vector(other_base[i]) for i in range(n)])
        else: # matrix repr
            other_base = other_base.copy()
        tr_base_mat, tr_base_denom = _toIntegerMatrix(
            other_base.inverse(self.base))
        denom = tr_base_denom * self.denominator
        mat_repr, numer = _toIntegerMatrix(
            tr_base_mat * self.mat_repr, 2)
        d = gcd.gcd(denom, numer)
        if d != 1:
            denom = ring.exact_division(denom, d)
            numer = ring.exact_division(numer, d)
        if numer != 1:
            mat_repr = numer * mat_repr
        return self.__class__([mat_repr, denom], self.number_field, other_base)

    def index(self):
        """
        return order of a residue group over base_module 
        (i.e. [base_module : module] or [module : base_module] ** (-1))
        """
        if not  self.mat_repr.ishnf:
            mat_repr = Submodule.fromMatrix(self.mat_repr.copy())
            mat_repr.toHNF()
        else:
            mat_repr = self.mat_repr
        n = self.base.column
        if mat_repr.column < n:
            raise ValueError("index is infinite")
        det = 1
        for i in range(1, n + 1):
            # for HNF, determinant = product of diagonal elements
            det *= mat_repr[i, i]
        if self.denominator != 1:
            return  rational.Rational(det, self.denominator ** n)
        else:
            return det

    def smallest_rational(self):
        """
        return the Z-generator of intersection of self and rational field
        """
        if not  self.mat_repr.ishnf:
            mat_repr = Submodule.fromMatrix(self.mat_repr.copy())
            mat_repr.toHNF()
        else:
            mat_repr = self.mat_repr
        theta_repr, pseudo_denom = _toIntegerMatrix(
            self.base * mat_repr)
        return rational.Rational(theta_repr.hermiteNormalForm()[1, 1],
            self.denominator * pseudo_denom)


class Ideal(Module):
    """
    for computing ideal with HNF (as Z-module).
    """
    def __init__(self, pair_mat_repr, number_field, base=None, ishnf=False):
        """
        Ideal is subclass of Module.
        Please refer to Module.__init__.__doc__
        """
        Module.__init__(self, pair_mat_repr, number_field, base, ishnf)

    issubideal = Module.issubmodule
    issuperideal = Module.issupermodule

    gcd = Module.__add__
    lcm = Module.intersect

    def __pow__(self, other):
        if other < 0:
            return Module.__pow__(self.inverse(), -other)
        elif other == 0:
            field = self.number_field
            n = field.degree
            int_basis = field.integer_ring()
            return Ideal([matrix.unitMatrix(n), 1], field, int_basis)
        else:
            return Module.__pow__(self, other)

    __pow__.__doc__ = Module.__pow__.__doc__   

    def inverse(self):
        """
        Return the inverse ideal of self
        """
        # Algorithm 4.8.21 in CCANT
        field = self.number_field
        # step 1 (Note that det(T)T^-1=(adjugate matrix of T))
        T = Ideal._precompute_for_different(field)
        d = int(T.determinant())
        inv_different = Ideal([T.adjugateMatrix(), 1], field,
                              field.integer_ring())
        # step 2
        inv_I = self * inv_different # already base is taken for integer_ring
        inv_I.toHNF()
        # step 3
        inv_mat, denom = _toIntegerMatrix(
            (inv_I.mat_repr.transpose() * T).inverse(), 0)
        numer = d * inv_I.denominator
        gcd_denom = gcd.gcd(denom, numer)
        if gcd_denom != 1:
            denom = ring.exact_division(denom, gcd_denom)
            numer = ring.exact_division(numer, gcd_denom)
        return Ideal([numer * inv_mat, denom], field, field.integer_ring())

    @classmethod
    def _precompute_for_different(cls, number_field):
        """
        Return T such that T^-1 represents HNF of inverse of different (codifferent)
        """
        field = number_field
        n = field.degree
        int_ring_mul = _base_multiplication(field.integer_ring(), field)
        T_lst = []
        for i in range(1, n + 1):
            each_T_lst = []
            for j in range(1, n + 1):
                each_T_lst.append(MatAlgNumber(
                    _symmetric_element(
                    i, j, int_ring_mul).compo, field.polynomial
                    ).trace())
            T_lst.append(each_T_lst)
        T = matrix.RingSquareMatrix(n, T_lst)
        return T        

    def twoElementRepresentation(self):
        # Exercise 31 + Proposition 4.7.8 in CCANT
        # (need prime ideal decompostion)
        raise NotImplementedError
        # return Ideal_with_generator([alpha, beta], self.number_field)

    def norm(self):
        """
        return the norm of self
        (Note that Norm(I)=[Z_K : I] for an ideal I and an integral ring Z_K)
        """
        # use method of returning integral ring for a number field
        return self.change_base_module(
            self.number_field.integer_ring()).index()

    def isIntegral(self):
        """
        determine whether self is integral ideal or not
        """
        mdl = self.change_base_module(self.number_field.integer_ring())
        return mdl.denominator == 1

    def isPrime(self):
        """
        determine whether self is prime ideal or not
        """
        if not self.isIntegral():
            return False
        size = self.number_field.degree
        nrm = self.norm()
        p = arith1.floorpowerroot(nrm, size)
        if p ** size != nrm or not prime.primeq(p):
            return False
        return None #undefined


class Ideal_with_generator(object):
    """
    ideal given as a generator
    (i.e. (a_1,...,a_n) = a_1 Z_K + ... + a_n Z_K)
    """
    def __init__(self, generator):
        """
        Ideal_with_generator(generator)
        generator: list of instances in BasicAlgNumber
                   (over same number_field)
        """
        self.generator = generator
        self.number_field = self.generator[0].field

    def __repr__(self):
        """
         formal representation
        """
        return '%s(%s)' % ( self.__class__.__name__, repr(self.generator) )

    def __str__(self):
        """
        simple representation
        """
        return str(self.generator)

    def __add__(self, other):
        """
        self + other as ideal
        other must be an instance of Ideal_with_generator
        """
        #if self.number_field != other.number_field:
        if self.number_field.polynomial != other.number_field.polynomial:
            raise NotImplementedError("use same number field")
        new_generator = self.generator[:]
        new_generator.extend(other.generator)
        return Ideal_with_generator(new_generator)
    
    def __mul__(self, other):
        """
        self * other as ideal
        other must be an instance of Ideal_with_generator
        """
        # if self.number_field != other.number_field:
        if self.number_field.polynomial != other.number_field.polynomial:
            raise NotImplementedError("use same number field")
        new_generator = [gen1 * gen2 
            for gen1 in self.generator for gen2 in other.generator]
        return Ideal_with_generator(new_generator)

    def __pow__(self, other):
        """
        self ** other (based on ideal multiplication)
        """
        if other <= 0:
            raise ValueError, "only support non-negative powering" 
        mul_part = self.copy()
        index = other
        while True:
            if index & 1:
                try:
                    sol *= mul_part
                except NameError:
                    sol = mul_part
            index >>= 1
            if not index:
                return sol
            else:
                mul_part *= mul_part

    def copy(self):
        """
        create copy of self
        """
        new_generator = []
        for gen in self.generator:
            new_generator.append(
                BasicAlgNumber([gen.coeff[:], gen.denom], gen.polynomial))
        return self.__class__(new_generator)

    def to_HNFRepresentation(self):
        """
        Transform self to the corresponding ideal as HNF representation
        """
        int_ring = self.number_field.integer_ring()
        n = self.number_field.degree
        polynomial = self.number_field.polynomial
        k = len(self.generator)
        # separate coeff and denom for generators and base
        gen_denom = reduce( gcd.lcm,
                            (self.generator[j].denom for j in range(k)) )
        gen_list = [gen_denom * self.generator[j] for j in range(k)]
        base_int, base_denom = _toIntegerMatrix(int_ring)
        
        base_alg = [BasicAlgNumber([base_int[j], 1], 
        polynomial) for j in range(1, n + 1)]
        repr_mat_list = []
        for gen in gen_list:
            for base in base_alg:
                new_gen = gen * base
                repr_mat_list.append(
                    vector.Vector([new_gen.coeff[j] for j in range(n)]))
        mat_repr = matrix.RingMatrix(n, n * k, repr_mat_list)
        new_mat_repr, numer = _toIntegerMatrix(mat_repr, 2)
        denom = gen_denom * base_denom
        denom_gcd = gcd.gcd(numer, denom)
        if denom_gcd != 1:
            denom = ring.exact_division(denom, denom_gcd)
            numer = ring.exact_division(numer, denom_gcd)
            if numer != 1:
                new_mat_repr = numer * new_mat_repr
        else:
            new_mat_repr = mat_repr
        return Ideal([new_mat_repr, denom], self.number_field)

    def twoElementRepresentation(self):
        # refer to [Poh-Zas] 
        # Algorithm 4.7.10 and Exercise 29 in CCANT
        """
        Reduce the number of generator to only two elements
        Warning: If self is not a prime ideal, this method is not efficient
        """
        k = len(self.generator)
        gen_denom = reduce( gcd.lcm, 
            (self.generator[j].denom for j in range(k)) )
        gen_list = [gen_denom * self.generator[j] for j in range(k)]
        int_I = Ideal_with_generator(gen_list)
        R = 1
        norm_I = int_I.norm()
        l_I = int_I.smallest_rational()
        if l_I.denominator > 1:
            raise ValueError, "input an integral ideal"
        else:
            l_I = l_I.numerator
        while True:
            lmda = [R for i in range(k)]
            while lmda[0] > 0:
                alpha = lmda[0] * gen_list[0]
                for i in range(1, k):
                    alpha += lmda[i] * gen_list[i]
                if gcd.gcd(
                  norm_I, ring.exact_division(
                  alpha.norm(), norm_I)
                  ) == 1 or gcd.gcd(
                  norm_I, ring.exact_division(
                  (alpha + l_I).norm(), norm_I)) == 1:
                    l_I_ori = BasicAlgNumber(
                        [[l_I] + [0] * (self.number_field.degree - 1),
                        gen_denom], self.number_field.polynomial)
                    alpha_ori = BasicAlgNumber([alpha.coeff, 
                        gen_denom], self.number_field.polynomial)
                    return Ideal_with_generator([l_I_ori, alpha_ori])
                for j in range(k)[::-1]:
                    if lmda[j] != -R:
                        lmda[j] -= 1
                        break
                    else:
                        lmda[j] = R
            R += 1

mono_method = ["index", "smallest_rational", "inverse", "norm"]
di_method = ["__eq__", "__ne__", "__contains__", "intersect", 
    "issubideal", "issuperideal"]

for new_method in mono_method:
    exec("def new_func(myself): " + 
        "return myself.to_HNFRepresentation()." + new_method + "()")
    exec("Ideal_with_generator." + new_method + " = new_func")
for new_method in di_method:
    exec("def new_func(myself, myother): " + 
        "return myself.to_HNFRepresentation()." + 
        new_method + "(myother.to_HNFRepresentation())")
    exec("Ideal_with_generator." + new_method + " = new_func")


######################################
#functions for internal computations
######################################
def _toIntegerMatrix(mat, option=0):
    """
    transform a (including integral) rational matrix to 
        some integral matrix as the following
    [option]
    0: return integral-matrix, denominator
       (mat = 1/denominator * integral-matrix)
    1: return integral-matrix, denominator, numerator
       (mat = numerator/denominator * reduced-int-matrix)
    2: return integral-matrix, numerator (assuming mat is integral)
       (mat = numerator * numerator-reduced-rational-matrix)
    """
    def getDenominator(ele):
        """
        get the denominator of ele
        """
        if rational.isIntegerObject(ele):
            return 1
        else: # rational
            return ele.denominator
    def getNumerator(ele):
        """
        get the numerator of ele
        """
        if rational.isIntegerObject(ele):
            return ele
        else: # rational
            return ele.numerator
    if isinstance(mat, vector.Vector):
        mat = mat.toMatrix(True)
    if option <= 1:
        denom = mat.reduce(
            lambda x, y: gcd.lcm(x, getDenominator(y)), 1)
        new_mat = mat.map(
            lambda x: getNumerator(x) * ring.exact_division(
            denom, getDenominator(x)))
        if option == 0:
            return Submodule.fromMatrix(new_mat), denom
    else:
        new_mat = mat
    numer = new_mat.reduce(
        lambda x, y: gcd.gcd(x, getNumerator(y)))
    if numer != 0:
        new_mat2 = new_mat.map(
            lambda x: ring.exact_division(getNumerator(x), numer))
    else:
        new_mat2 = new_mat
    if option == 1:
        return Submodule.fromMatrix(new_mat2), denom, numer
    else:
        return Submodule.fromMatrix(new_mat2), numer

def _symmetric_element(i, j, mat):
    """
    get (i, j)-element from lst ordered (1, 1), (1, 2), (2, 2), ...
    we assume that
    1: (j, i)-element is same as (i, j)-element (i.e. symmetric)
    2: index of mat starts 1 (i.e. for Matrix only)
    """
    if i <= j:
        return mat[int((j - 1) * j / 2 + i)]
    else:
        return mat[int((i - 1) * i / 2 + j)]

def _base_multiplication(base, number_field):
    """
    return [base[i] * base[j]] (as a numberfield element)
    this is a precomputation for computing _module_mul
    """
    n = number_field.degree
    polynomial = number_field.polynomial
    base_int, base_denom = _toIntegerMatrix(base)
    base_alg = [BasicAlgNumber([base_int[j], 1], 
        polynomial) for j in range(1, n + 1)]
    base_ij_list = []
    for j in range(n):
        for i in range(j + 1):
            base_ij = base_alg[i] * base_alg[j]
            base_ij_list.append(
                vector.Vector([base_ij.coeff[k] * 
                rational.Rational(1, base_denom ** 2) for k in range(n)]))
    base_ij = base.inverse(
        matrix.FieldMatrix(n, len(base_ij_list), base_ij_list))
    return base_ij
