from __future__ import division

import nzmath.ring as ring
import nzmath.vector as vector


class Matrix(object):
    """
    Matrix is a class for matrices.
    """

    def __init__(self, row, column, compo=0, coeff_ring=0):
        """
        Matrix(row, column [,components, coeff_ring])
        """
        self._initialize(row, column, compo, coeff_ring)
        self._selectMatrix()

    def _initialize(self, row, column, compo=0, coeff_ring=0):
        """
        initialize matrix.
        """
        if (isinstance(row, (int, long))
            and isinstance(column, (int, long))
            and row > 0
            and column > 0 ): # row and column check
            self.row = row
            self.column = column
            self.compo = []
            if isinstance(compo, (ring.Ring, int)):
                coeff_ring = compo
                compo = 0
            if compo == 0:
                zero = 0
                if coeff_ring != 0:
                    coeff_ring = _getRing(coeff_ring)
                    zero = coeff_ring.zero
                for i in range(self.row):
                    self.compo.append([zero] * self.column)
            else:
                if isinstance(compo[0], list):
                    if len(compo) != self.row:
                        raise ValueError("len(compo) " + 
                                         "is not match the matrix size")
                    for i in range(self.row):
                        if len(compo[i]) != self.column:
                            raise ValueError("len(compo[%d]) " % i + 
                                             "is not match the matrix size")
                        self.compo.append(compo[i][:])
                elif isinstance(compo[0], tuple):
                    if len(compo) != self.column:
                        raise ValueError("len(compo) " +
                                         "is not match the matrix size")
                    self.compo = [[] for i in range(self.row)]
                    for i in range(self.column):
                        if len(compo[i]) != self.row:
                            raise ValueError("len(compo[%d]) " % i + 
                                             "is not match the matrix size")
                        j = 0
                        for ele in compo[i]:
                            self.compo[j].append(ele)
                            j += 1
                elif isinstance(compo[0], vector.Vector):
                    if len(compo) != self.column:
                        raise ValueError("len(compo) " +
                                         "is not match the matrix size")
                    self.compo = [[] for i in range(self.row)]
                    for i in range(self.column):
                        if len(compo[i]) != self.row:
                            raise ValueError("len(compo[%d]) " % i +
                                             "is not match the matrix size")
                        j = 0
                        for ele in compo[i].compo:
                            self.compo[j].append(ele)
                            j += 1
                else:
                    if (len(compo) != self.row * self.column):
                        raise ValueError("len(compo) " +
                                         "is not match the matrix size")
                    for i in range(self.row):
                        self.compo.append(
                        compo[self.column * i : self.column * (i + 1)])
            if coeff_ring == 0:
                self.coeff_ring = ring.getRing(self.compo[0][0])
            else:
                self.coeff_ring = coeff_ring = _getRing(coeff_ring)
                if not(isinstance(ring.getRing(self.compo[0][0]), coeff_ring.__class__)):
                    compos = []
                    for i in range(self.row):
                        compos.append(map(coeff_ring.createElement, self.compo[i]))
                    self.compo = compos
        else:
            raise ValueError("invalid value for matrix size")

    def _selectMatrix(self):
        """
        Select Matrix class.
        """
        if self.coeff_ring.isfield():
            if self.row == self.column:
                self.__class__ = FieldSquareMatrix
            else:
                self.__class__ = FieldMatrix
        else:
            if self.row == self.column:
                self.__class__ = RingSquareMatrix
            else:
                self.__class__ = RingMatrix

    def __getitem__(self, index):
        """
        M[i,j] : Return (i,j)-component of M.
        M[i] <==> M.getColumn(i)
        """
        if isinstance(index, tuple):
            return self.compo[index[0] - 1][index[1] - 1]
        elif isinstance(index, (int, long)):
            return self.getColumn(index)
        else:
            raise IndexError("Matrix invalid index: %s" % index)

    def __setitem__(self, key, value):
        """
        M[i,j] = a   :   Substitute a for (i,j)-component of M.
        M[i] = V <==> M.setColumn(i, V)
        """
        if isinstance(key, tuple):
            self.compo[key[0] - 1][key[1] - 1] = value
        elif isinstance(key, (int, long)):
            self.setColumn(key, value)
        else:
            raise TypeError(self.__setitem__.__doc__)

    def __eq__(self, other):
        """
        Check self == other.
        self == 0 means whether self == zeromatrix or not.
        """
        if isinstance(other, Matrix):
            if (self.row != other.row) or (self.column != other.column):
                return False
            for i in range(self.row):
                for j in range(self.column):
                    if self.compo[i][j] != other.compo[i][j]:
                        return False
            return True
        elif isinstance(other, int) and other == 0: # zero matrix ?
            return not bool(self)
        else:
            return False

    def __hash__(self):
        val = 0
        for i in range(self.row):
            for j in range(self.column):
                val += hash(self.compo[i][j])
        return val

    def __ne__(self, other):
        """
        Check self != other.
        self != 0 means whether self != zeromatrix or not.
        """
        return not (self == other)

    def __nonzero__(self):
        """
        Check self != zeromatrix
        """
        for i in range(self.row):
            for j in range(self.column):
                if self.compo[i][j]:
                    return True
        return False

    def __contains__(self, item):
        """
        Check whether item == one of self's element.
        """
        for i in range(self.row):
            if item in self.compo[i]:
                return True
        return False

    def __repr__(self):
        return_str = ""
        for i in range(self.row):
            return_str += str(self.compo[i])
            if i + 1 != self.row:
                return_str += "+"
        return return_str

    def __str__(self):
        return_str = ""
        width = [1] * self.column      # width of each column
        for j in range(self.column):
            for i in range(self.row):
                if len(str(self.compo[i][j])) > width[j]:
                    width[j] = len(str(self.compo[i][j]))
        for i in range(self.row):
            for j in range(self.column):
                return_str += "%*s " % (width[j], str(self.compo[i][j]))
            return_str += "\n"
        return return_str[:-1]

    def __call__(self, arg):
        """
        Return matrix applied __call__ to each elements.
        """
        sol = []
        for i in range(self.row):
            for j in range(self.column):
                ele = self.compo[i][j]
                if callable(ele):
                    sol.append(ele(arg))
                else:
                    sol.append(ele)
        return createMatrix(self.row, self.column, sol)


    def map(self, function):
        """
        Return matrix applied function to all self elements.
        """
        compos = []
        for i in range(self.row):
            compos = compos + (map(function, self.compo[i]))
        return createMatrix(self.row, self.column, compos)

    def reduce(self, function, initializer=None):
        """
        Return value applied binomial function to all self elements.
        """
        if initializer:
            val = reduce(function, self.compo[0], initializer)
        else:
            val = reduce(function, self.compo[0])
        for i in range(1, self.row):
            val = reduce(function, self.compo[i], val)
        return val

    def copy(self):
        """
        Create a copy of the instance.
        """
        compos = []
        for i in range(self.row):
            for j in range(self.column):
                compos.append(self.compo[i][j])
        Mat = self.__class__(self.row, self.column, compos, self.coeff_ring)
        return Mat
    
    def set(self, compo):
        """
        set(compo) : Substitute list for components
        """
        if (len(compo) != self.row * self.column):
            raise ValueError("len(compo) is not match the matrix size")
        for i in range(self.row):
            self.compo[i] = compo[self.column*i : self.column*(i + 1)]

    def setRow(self, m, arg):
        """
        setRow(m, new_row) : new_row should be a list/Vector
        """
        if isinstance(arg, list):
            if (len(arg) != self.column):
                raise vector.VectorSizeError(
                    "len(compo) is not match the row size")
            self.compo[m - 1] = arg[:]
        elif isinstance(arg, vector.Vector):
            if (len(arg) != self.column):
                raise vector.VectorSizeError(
                    "len(compo) is not match the row size")
            self.compo[m - 1] = arg.compo[:]
        else:
            raise TypeError(self.setRow.__doc__)

    def setColumn(self, n, arg):
        """
        setColumn(n, new_column) : new_column should be a list/Vector
        """
        if isinstance(arg, list):
            if (len(arg) != self.row):
                raise ValueError("len(compo) is not match the column size")
            for i in range(self.row):
                self.compo[i][n - 1] = arg[i]
        elif isinstance(arg, vector.Vector):
            if (len(arg) != self.row):
                raise ValueError("len(compo) is not match the column size")
            for i in range(self.row):
                self.compo[i][n - 1] = arg.compo[i]
        else:
            raise TypeError(self.setColumn.__doc__)

    def getRow(self, i):
        """
        getRow(i) : Return i-th row in form of Matrix
        """
        return vector.Vector(self.compo[i - 1][:])

    def getColumn(self, j):
        """
        getColumn(j) : Return j-th column in form of Matrix
        """
        column = []
        for k in range(self.row):
            column.append(self.compo[k][j - 1])
        return vector.Vector(column)

    def swapRow(self, m1, m2):
        """
        swapRow(m1, m2) : Swap self's m1-th row and m2-th row.
        """
        tmp = self.compo[m1 - 1][:]
        self.compo[m1 - 1] = self.compo[m2 - 1][:]
        self.compo[m2 - 1] = tmp[:]

    def swapColumn(self, n1, n2):
        """
        swapColumn(n1, n2) : Swap self's n1-th column and n2-th column.
        """
        for k in range(self.row):
            tmp = self.compo[k][n1 - 1]
            self.compo[k][n1 - 1] = self.compo[k][n2 - 1]
            self.compo[k][n2 - 1] = tmp

    def insertRow(self, i, arg):
        """
        insertRow(i, new_row) : added new_row
        new_row can be a list or a Matrix
        """
        if isinstance(arg, list):
            if self.column != len(arg):
                raise vector.VectorSizeError()
            self.compo.insert(i - 1, arg)
            self.row += 1
        elif isinstance(arg, vector.Vector):
            if self.column != len(arg):
                raise vector.VectorSizeError()
            self.compo.insert(i - 1, arg.compo)
            self.row += 1
        elif isinstance(arg, Matrix):
            if self.column != arg.column:
                raise MatrixSizeError()
            self.compo += arg.compo
            self.row += arg.row
        else:
            raise TypeError()
        self._selectMatrix()

    def insertColumn(self, j, arg):
        """
        insertColumn(j, new_column) : added new_column
        new_column can be a list or a Matrix
        """
        if isinstance(arg, list):
            if self.row != len(arg):
                raise vector.VectorSizeError()
            for k in range(self.row):
                self.compo[k].insert(j-1, arg[k])
            self.column += 1
        elif isinstance(arg, vector.Vector):
            if self.row != len(arg):
                raise vector.VectorSizeError()
            for k in range(self.row):
                self.compo[k].insert(j-1, arg.compo[k])
            self.column += 1
        elif isinstance(arg, Matrix):
            if self.row != arg.row:
                raise MatrixSizeError()
            for k in range(arg.row):
                self.compo[k] = \
                self.compo[k][:j - 1] + arg.compo[k] + self.compo[k][j - 1:]
            self.column += arg.column
        else:
            raise TypeError()
        self._selectMatrix()

    def extendRow(self, arg):
        """
        extendRow(new_row) : join new_row in vertical way.
        """
        self.insertRow(self.row + 1, arg)

    def extendColumn(self, arg):
        """
        extendColumn(new_column) : join new_column in horizontal way.
        """
        self.insertColumn(self.column + 1, arg)

    def deleteRow(self, i):
        """
        deleteRow(i) : deleted i-th row
        """
        self.row -= 1
        del self.compo[i - 1]
        self._selectMatrix()

    def deleteColumn(self, j):
        """
        deleteColumn(j) : deleted j-th column
        """
        self.column -= 1
        for k in range(self.row):
            del self.compo[k][j - 1]
        self._selectMatrix()


    def transpose(self):
        """
        Return transposed matrix of self.
        """
        trans = []
        for j in range(1, self.column + 1):
            for i in range(1, self.row + 1):
                trans.append(self[i, j])
        return self.__class__(self.column, self.row, trans, self.coeff_ring)

    def getBlock(self, i, j, row, column=None):
        """
        Return a block whose size is row*column, (1,1)-element is self[i,j].
        """
        if column == None:
            column = row
        if i + row - 1 > self.row or j + column - 1 > self.column:
            raise MatrixSizeError()
        mat = []
        for k in range(i - 1, i + row - 1):
            mat.extend(self.compo[k][j - 1:j + column - 1])
        return createMatrix(row, column, mat, self.coeff_ring)

    def subMatrix(self, I, J=None):
        """
        Return submatrix whose element is self[i, j] for i in I and j in J.
        If I, J is not index(list or tuple) but integer,
         return submatrix which deleted I-th row and J-th column from self.
        """
        if J == None:
            J = I
        if isinstance(I, (int, long)) and isinstance(J, (int, long)):
            mat = self.copy()
            mat.deleteRow(I)
            mat.deleteColumn(J)
            return mat
        else:
            mat = []
            for i in I:
                for j in J:
                    mat.append(self[i, j])
            return createMatrix(len(I), len(J), mat, self.coeff_ring)

    def toMatrix(self, flag=True):
        """
        return the matrix (self)
        (this method is for compatibility to vector)
        """
        return self


class SquareMatrix(Matrix):
    """
    SquareMatrix is a class for square matrices.
    """

    def __init__(self, row, column=0, compo=0, coeff_ring=0):
        """
        SquareMatrix(row [, column ,components, coeff_ring])
        SquareMatrix must be row == column .
        """
        self._initialize(row, column, compo, coeff_ring)
        self._selectMatrix()

    def _initialize(self, row, column=0, compo=0, coeff_ring=0):
        """
        initialize matrix.
        """
        if isinstance(compo, (ring.Ring, int)):
            coeff_ring = compo
            compo = 0
        if isinstance(column, list):
            compo = column
            column = row
        elif isinstance(column, ring.Ring):
            coeff_ring = column
            column = row
        elif column == 0:
            column = row
        if row != column:
            raise ValueError("not square matrix")
        Matrix._initialize(self, row, column, compo, coeff_ring)

    def isUpperTriangularMatrix(self):
        """
        Check whether self is upper triangular matrix or not.
        """
        for j in range(1, self.column + 1):
            for i in range(j + 1, self.row + 1):
                if self[i, j]:
                    return False
        return True

    def isLowerTriangularMatrix(self):
        """
        Check whether self is lower triangular matrix or not.
        """
        for i in range(1, self.row + 1):
            for j in range(i + 1, self.column + 1):
                if self[i, j]:
                    return False
        return True

    def isDiagonalMatrix(self):
        """
        Check whether self is diagonal matrix or not.
        """
        return self.isUpperTriangularMatrix() and self.isLowerTriangularMatrix()

    def isScalarMatrix(self):
        """
        Check whether self is scalar matrix or not.
        Scalar matrix is matrix which is unit matrix times some scalar.
        """
        if not(self.isDiagonalMatrix()):
            return False
        chk = self[1, 1]
        for i in range(2, self.row + 1):
            if self[i, i] != chk:
                return False
        return True

    def isSymmetricMatrix(self):
        """
        Check whether self is symmetric matrix or not.
        """
        for i in range(1, self.row + 1):
            for j in range(i + 1, self.column + 1):
                if self[i, j] != self[j, i]:
                    return False
        return True


class RingMatrix(Matrix):
    """
    RingMatrix is a class for matrices whose elements are in ring.
    """

    def __init__(self, row, column, compo=0, coeff_ring=0):
        """
        RingMatrix(row, column [,components, coeff_ring])
        """
        self._initialize(row, column, compo, coeff_ring)
        self._selectMatrix()

    def _selectMatrix(self):
        """
        Select Matrix class.
        """
        if self.row == self.column:
            self.__class__ = RingSquareMatrix
        else:
            self.__class__ = RingMatrix

    def __add__(self, other):
        """
        Return matrix addition.
        """
        if (self.row != other.row) or (self.column != other.column):
            raise MatrixSizeError()
        sums = []
        for i in range(1, self.row + 1):
            for j in range(1, self.column + 1):
                sums.append(self[i, j] + other[i, j])
        return createMatrix(self.row, self.column, sums,
                  self.coeff_ring.getCommonSuperring(other.coeff_ring))

    def __sub__(self, other):
        """
        Return matrix subtraction.
        """
        if (self.row != other.row) or (self.column != other.column):
            raise MatrixSizeError()
        diff = []
        for i in range(1, self.row + 1):
            for j in range(1, self.column + 1):
                diff.append(self[i, j] - other[i, j])
        return createMatrix(self.row, self.column, diff,
                  self.coeff_ring.getCommonSuperring(other.coeff_ring))

    def __mul__(self, other):
        """
        multiplication with a Matrix or a scalar
        """
        if isinstance(other, Matrix):
            if self.column != other.row:
                raise MatrixSizeError()
            product = []
            for i in range(1, self.row + 1):
                for j in range(1, other.column + 1):
                    part_product = self[i, 1] * other[1, j]
                    for k in range(2, self.column + 1):
                        part_product = part_product + self[i, k] * other[k, j]
                    product.append(part_product)
            return createMatrix(self.row, other.column, product,
                      self.coeff_ring.getCommonSuperring(other.coeff_ring))
        elif isinstance(other, vector.Vector):
            if self.column != len(other):
                raise vector.VectorSizeError()
            product = []
            for i in range(1, self.row + 1):
                part_product = self[i, 1] * other[1]
                for j in range(2, self.column + 1):
                    part_product = part_product + self[i, j] * other[j]
                product.append(part_product)
            return vector.Vector(product)
        else: #scalar mul
            try:
                rings = self.coeff_ring.getCommonSuperring(ring.getRing(other))
            except:
                return NotImplemented
            product = []
            for i in range(1, self.row + 1):
                for j in range(1, self.column + 1):
                    product.append(self[i, j] * other)
            return createMatrix(self.row, self.column, product, rings)

    def __rmul__(self, other):
        if isinstance(other, Matrix):
            if other.column != self.row:
                raise MatrixSizeError()
            product = []
            for i in range(1, other.row + 1):
                for j in range(1, self.column + 1):
                    part_product = other[i, 1] * self[1, j]
                    for k in range(2, self.row + 1):
                        part_product = part_product + other[i, k] * self[k, j]
                    product.append(part_product)
            return createMatrix(other.row, self.column, product,
                     self.coeff_ring.getCommonSuperring(other.coeff_ring))
        elif isinstance(other, vector.Vector):
            if self.row != len(other):
                raise vector.VectorSizeError()
            product = []
            for j in range(1, self.column + 1):
                part_product = other[1] * self[1, j]
                for i in range(2, self.row + 1):
                    part_product = part_product + other[i] * self[i, j]
                product.append(part_product)
            return vector.Vector(product)
        else:
            try:
                rings = self.coeff_ring.getCommonSuperring(ring.getRing(other))
            except:
                return NotImplemented
            product = []
            for i in range(1, self.row + 1):
                for j in range(1, self.column + 1):
                    product.append(self[i, j] * other)
            return createMatrix(self.row, self.column, product, rings)

    def __mod__(self, other):
        """
        return self modulo other.
        """
        if not bool(other):
            raise ZeroDivisionError()
        mod = []
        for i in range(1, self.row + 1):
            for j in range(1, self.column + 1):
                mod.append(self[i, j] % other)
        return createMatrix(self.row, self.column, mod, self.coeff_ring)

    def __pos__(self):
        """
        return copy of self.
        """
        return self.copy()

    def __neg__(self):
        """
        return -self.
        """
        one = self.coeff_ring.one
        try:
            minus_one = -one
        except:
            minus_one = self.coeff_ring.zero - one
        return self.map(lambda ele: minus_one * ele)

    def getCoefficientRing(self):
        """
        Set and return coefficient ring.
        """
        if not hasattr(self, "_coeff_ring"):
            scalars = None
            for i in range(1, self.row + 1):
                for j in range(1, self.column + 1):
                    cring = ring.getRing(self[i, j])
                    if scalars is None or \
                    scalars != cring and scalars.issubring(cring):
                        scalars = cring
                    elif not scalars.issuperring(cring):
                        scalars = scalars.getCommonSuperring(cring)
            if scalars != self.coeff_ring or scalars.issubring(self.coeff_ring):
                scalars = self.coeff_ring
            elif not scalars.issuperring(self.coeff_ring):
                scalars = scalars.getCommonSuperring(self.coeff_ring)
            self._coeff_ring = self.coeff_ring = scalars
        return self._coeff_ring

    def toFieldMatrix(self):
        """RingMatrix -> FieldMatrix"""
        self.__class__ = FieldMatrix
        self.coeff_ring = self.coeff_ring.getQuotientField()

    def toSubspace(self, isbasis=None):
        """RingMatrix -> Subspace"""
        self.toFieldMatrix()
        self.toSubspace(isbasis)

    def hermiteNormalForm(self, non_zero=False):
        # Algorithm 2.4.5 in CCANT (fixed for parameter l)
        """
        Find the Hermite Normal Form for integer matrix.
        If non_zero is True, return non-zero columns for self's HNF
        """
        A = self.copy()
        rings = self.coeff_ring
        # step 1 [Initialize]
        j0 = k = self.column
        for i in range(1, self.row + 1)[::-1]:
            while 1:
                # step 2 [Check zero]
                for j in range(1, j0)[::-1]:
                    if bool(A[i, j]):
                        break
                else: # j==1
                    break
                # step 3 [Euclidean step]
                u, v, d = rings.extgcd(A[i, k], A[i, j])
                A_ik = ring.exact_division(A[i, k], d)
                A_ij = ring.exact_division(A[i, j], d)
                B = u * A[k] + v * A[j]
                A[j] = A_ik * A[j] - A_ij * A[k]
                A[k] = B
            # step4 [Final reductions]
            b = A[i, k]
            if b < 0:
                A[k] = -A[k]
                b = -b
            if not bool(b):
                k += 1
            else:
                for j in range(k + 1, self.column + 1):
                    q = A[i, j] // b
                    A[j] = A[j] - q * A[k]
            # step 5 [Finished?]
            k -= 1
            j0 = k
            if k == 0:
                break
            # go to step 2
        if non_zero:
            if k == self.column: # zero module
                return None
            W = createMatrix(
                self.row, self.column - k,
                [A[j + k] for j in range(1, self.column - k + 1)])
            return W
        else:
            return A

    HNF = hermiteNormalForm

    def _SimplifyHNF(self):
        """
        This method is a common process used in 
        extHNF() and kernelAsModule()
        """
        A = self.copy()
        U = unitMatrix(A.column, A.coeff_ring)
        rings = self.coeff_ring
        # step 1 [Initialize]
        j0 = k = self.column
        for i in range(1, self.row + 1)[::-1]:
            while 1:
                # step 2 [Check zero]
                for j in range(1, j0)[::-1]:
                    if bool(A[i, j]):
                        break
                else: # j==1
                    break
                # step 3 [Euclidean step]
                u, v, d = rings.extgcd(A[i, k], A[i, j])
                A_ik = ring.exact_division(A[i, k], d)
                A_ij = ring.exact_division(A[i, j], d)
                B = u * A[k] + v * A[j]
                A[j] = A_ik * A[j] - A_ij * A[k]
                A[k] = B
                B = u * U[k] + v * U[j]
                U[j] = A_ik * U[j] - A_ij * U[k]
                U[k] = B
            # step4 [Final reductions]
            b = A[i, k]
            if b < 0:
                A[k] = -A[k]
                U[k] = -U[k]
                b = -b
            if not bool(b):
                k += 1
            else:
                for j in range(k + 1, self.column + 1):
                    q = A[i, j] // b
                    A[j] = A[j] - q * A[k]
                    U[j] = U[j] - q * U[k]
            # step 5 [Finished?]
            k -= 1
            j0 = k
            if k == 0:
                break
            # go to step 2
        return (U, A, k)

    def exthermiteNormalForm(self, non_zero=False):
        # Modified Algorithm 2.4.5 in CCANT
        """
        Find the Hermite Normal Form M for integer matrix.
        Computing U which satisfied M=self*U.
        Return matrices tuple,(U, M).
        """
        U, A, k = self._SimplifyHNF()
        if non_zero:
            if k == self.column: # zero module
                return None
            new_A = createMatrix(
                self.row, self.column - k,
                [A[j + k] for j in range(1, self.column - k + 1)])
            new_U = createMatrix(
                self.column, self.column - k,
                [U[j + k] for j in range(1, self.column - k + 1)])
            return (new_U, new_A)
        else:
            return (U, A)

    extHNF = exthermiteNormalForm

    def kernelAsModule(self):
        """
        Compute kernel as Z-module.
        """
        U, A, k = self._SimplifyHNF()
        if k == 0:
            return None
        else:
            ker = createMatrix(
                self.column, k, [U[j] for j in range(1, k + 1)])
            return ker


class RingSquareMatrix(SquareMatrix, RingMatrix, ring.RingElement):
    """
    RingSquareMatrix is a class for square matrices whose elements are in ring.
    """

    def __init__(self, row, column=0, compo=0, coeff_ring=0):
        """
        RingSquareMatrix(row [, column ,components, coeff_ring])
        RingSquareMatrix must be row == column .
        """
        self._initialize(row, column, compo, coeff_ring)

    def __pow__(self, other):
        """
        powering self to integer.
        """
        n = +other
        if not isinstance(n, (int, long)):
            raise TypeError("index must be an integer")
        power = unitMatrix(self.row, self.coeff_ring)
        # check n
        if n == 0:
            return power
        if n > 0:
            z = self.copy()
        else:
            if hasattr(self, "inverse"):
                n = abs(n)
                z = self.inverse()
            else:
                raise NoInverse()
        while 1:
            if n & 1:
                power = power * z
            n //= 2
            if n == 0:
                return power
            z = z * z

    def toFieldMatrix(self):
        """RingSquareMatrix -> FieldSquareMatrix"""
        self.__class__ = FieldSquareMatrix
        self.coeff_ring = self.coeff_ring.getQuotientField()

    def getRing(self):
        """
        Return matrix ring of self.
        """
        return MatrixRing.getInstance(self.row, self.getCoefficientRing())

    def isOrthogonalMatrix(self):
        """
        Check whether self is orthogonal matrix or not.
        Orthogonal matrix satisfies M*M^T equals unit matrix.
        """
        return self * self.transpose() == unitMatrix(self.row, self.coeff_ring)

    def isAlternatingMatrix(self):
        """
        Check whether self is alternating matrix or not.
        Alternating (skew symmetric, or antisymmetric) matrix satisfies M=-M^T.
        """
        for i in range(1, self.row + 1):
            for j in range(i, self.column + 1):
                if self[i, j] != -self[j, i]:
                    return False
        return True

    isAntisymmetricMatrix = isAlternatingMatrix
    isSkewsymmetricMatrix = isAlternatingMatrix

    def isSingular(self):
        """
        Check determinant == 0 or not.
        """
        return not bool(self.determinant())

    def trace(self):
        """
        Return trace of self.
        """
        trace = self.coeff_ring.zero
        for i in range(1, self.row + 1):
            trace = trace + self[i, i]
        return trace

    def determinant(self): # Algorithm 2.2.6 of Cohen's book
        """
        Return determinant of self.
        """
        M = self.copy()
        n = self.row
        c = self.coeff_ring.one
        sign = True
        for k in range(1, n):
            p = M[k, k]
            if not bool(p): # p==0
                i = k + 1
                while not bool(M[i, k]):
                    if i == n:
                        return self.coeff_ring.zero
                    else:
                        i += 1
                for j in range(k, n + 1):
                    tmp = M[i, j]
                    M[i, j] = M[k, j]
                    M[k, j] = tmp
                sign = not(sign)
                p = M[k, k]
            for i in range(k + 1, n + 1):
                for j in range(k + 1, n + 1):
                    t = p * M[i, j] - M[i, k] * M[k, j]
                    M[i, j] = ring.exact_division(t, c)
            c = p
        if sign:
            return M[n, n]
        else:
            return -M[n, n]

    def cofactor(self, i, j):
        """
        Return (i, j)-cofactor of self.
        """
        cofactor = (self.subMatrix(i, j)).determinant()
        if (i + j) & 1:
            cofactor = -cofactor
        return cofactor

    def commutator(self, other):
        """
        Return commutator defined as follows:
        [self, other] = self * other - other * self .
        """
        return self * other - other * self

    def characteristicMatrix(self):
        """
        Return the characteristic matrix (i.e. xI-A) of self.
        """
        import nzmath.poly.uniutil as uniutil
        x = uniutil.polynomial({1:1}, self.coeff_ring)
        return x * unitMatrix(self.row, x.getRing()) - self

    def _characteristicPolyList(self): # Algorithm 2.2.7 of Cohen's book
        """
        for characteristicPolynomial, adjugateMatrix
        
        Assume self.row >= 2.
        """
        unit = unitMatrix(self.row, self.coeff_ring)
        coeff = [self.coeff_ring.one, -self.trace()]
        C = self + coeff[-1] * unit
        i = 2
        while i < self.row:
            C = self * C
            coeff.append(ring.exact_division(-C.trace(), i))
            C = C + coeff[-1] * unit
            i += 1
        coeff.append(ring.exact_division(-(self * C).trace(), i))
        coeff.reverse()
        return coeff, C

    def characteristicPolynomial(self):
        """
        characteristicPolynomial() -> Polynomial
        """
        import nzmath.poly.uniutil
        genPoly = nzmath.poly.uniutil.polynomial
        if self.row == 1:
            rings = self.coeff_ring
            return genPoly({0:-self.trace(), 1:rings.one}, rings)
        coeff = self._characteristicPolyList()[0]
        return genPoly(dict(enumerate(coeff)), self.coeff_ring)

    def adjugateMatrix(self):
        """
        Return adjugate(classical adjoint) matrix.
        """
        if self.row == 1:
            return unitMatrix(self.row, self.coeff_ring)
        C = self._characteristicPolyList()[1]
        if self.row & 1:
            return C
        else:
            return -C

    def cofactorMatrix(self):
        """
        Return cofactor matrix.
        """
        return self.adjugateMatrix().transpose()

    cofactors = cofactorMatrix

    def smithNormalForm(self):# Algorithm 2.4.14 of Cohen's book
        """
        Find the Smith Normal Form for square non-singular integral matrix.
        Return the list of diagonal elements.
        """
        M = self.copy()
        n = M.row
        R = M.determinant()
        rings = self.coeff_ring
        if not bool(R):
            raise ValueError("Don't input singular matrix")
        if R < 0:
            R = -R
        lst = []
        while n != 1:
            j = n
            c = 0
            while j != 1:
                j -= 1
                if M[n, j]:
                    u, v, d = rings.extgcd(M[n, j], M[n, n])
                    B = v * M.getColumn(n) + u * M.getColumn(j)
                    M_nn = ring.exact_division(M[n, n], d)
                    M_nj = ring.exact_division(M[n, j], d)
                    M.setColumn(j, ((M_nn * M.getColumn(j)
                                     - M_nj * M.getColumn(n)) % R))
                    M.setColumn(n, (B % R))
            j = n
            while j != 1:
                j -= 1
                if M[j, n]:
                    u, v, d = rings.extgcd(M[j, n], M[n, n])
                    B = v * M.getRow(n) + u * M.getRow(j)
                    M_nn = ring.exact_division(M[n, n], d)
                    M_jn = ring.exact_division(M[j, n], d)
                    M.setRow(j, ((M_nn * M.getRow(j)
                                  - M_jn * M.getRow(n)) % R))
                    M.setRow(n, (B % R))
                    c += 1
            if c <= 0:
                b = M[n, n]
                flag = False
                if not bool(b):
                    b = R
                for k in range(1, n):
                    for l in range(1, n):
                        if (M[k, l] % b):
                            M.setRow(n, M.getRow(n) + M.getRow(k))
                            flag = True
                if not flag:
                    dd = rings.gcd(M[n, n], R)
                    lst.append(dd)
                    R = ring.exact_division(R, dd)
                    n -= 1
        dd = rings.gcd(M[1, 1], R)
        lst.append(dd)
        lst.reverse()
        return lst

    SNF = smithNormalForm
    elementary_divisor = smithNormalForm

    def extsmithNormalForm(self):
        """
        Find the Smith Normal Form M for square matrix,
        Computing U,V which satisfied M=U*self*V.
        Return matrices tuple,(U,V,M).
        """
        M = self.copy()
        n = M.row
        U = unitMatrix(M.row, M.coeff_ring)
        V = unitMatrix(M.row, M.coeff_ring)
        rings = self.coeff_ring
        while n != 1:
            j = n
            c = 0
            while j != 1:
                j -= 1
                if M[n, j]:
                    u, v, d = rings.extgcd(M[n, j], M[n, n])
                    M_nn = ring.exact_division(M[n, n], d)
                    M_nj = ring.exact_division(M[n, j], d)
                    B = v * M.getColumn(n) + u * M.getColumn(j)
                    M.setColumn(j, (M_nn * M.getColumn(j)
                                    - M_nj * M.getColumn(n)))
                    M.setColumn(n, B)
                    B = v * V.getColumn(n) + u * V.getColumn(j)
                    V.setColumn(j, (M_nn * V.getColumn(j)
                                    - M_nj * V.getColumn(n)))
                    V.setColumn(n, B)
            j = n
            while j != 1:
                j -= 1
                if M[j, n]:
                    u, v, d = rings.extgcd(M[j, n], M[n, n])
                    M_nn = ring.exact_division(M[n, n], d)
                    M_jn = ring.exact_division(M[j, n], d)
                    B = v * M.getRow(n) + u * M.getRow(j)
                    M.setRow(j, (M_nn * M.getRow(j) - M_jn * M.getRow(n)))
                    M.setRow(n, B)
                    B = v * U.getRow(n) + u * U.getRow(j)
                    U.setRow(j, (M_nn * U.getRow(j) - M_jn * U.getRow(n)))
                    U.setRow(n, B)
                    c += 1
            if c <= 0:
                b = M[n, n]
                flag = False
                for k in range(1, n):
                    for l in range(1, n):
                        if (M[k, l] % b):
                            M.setRow(n, M.getRow(n) + M.getRow(k))
                            U.setRow(n, U.getRow(n) + U.getRow(k))
                            flag = True
                if not flag:
                    n -= 1
        for j in range(1, M.column + 1):
            if M[j, j] < 0:
                V[j] = -V[j]
                M[j, j] = -M[j, j]
        return (U, V, M)

    extSNF = extsmithNormalForm


class FieldMatrix(RingMatrix):
    """
    FieldMatrix is a class for matrices whose elements are in field.
    """

    def __init__(self, row, column, compo=0, coeff_ring=0):
        """
        FieldMatrix(row, column [,components, coeff_ring])
        """
        self._initialize(row, column, compo, coeff_ring)
        self._selectMatrix()

    def _initialize(self, row, column, compo=0, coeff_ring=0):
        """initialize matrix"""
        RingMatrix._initialize(self, row, column, compo, coeff_ring)
        if not self.coeff_ring.isfield():
            self.coeff_ring = self.coeff_ring.getQuotientField()

    def _selectMatrix(self):
        """
        Select Matrix class.
        """
        if self.__class__ != Subspace:
            if self.row == self.column:
                self.__class__ = FieldSquareMatrix
            else:
                self.__class__ = FieldMatrix

    def __truediv__(self, other):
        """
        division by a scalar.
        """
        return ring.inverse(other) * self

    __div__ = __truediv__ # backward compatibility?

    def toSubspace(self, isbasis=None):
        """FieldMatrix -> Subspace"""
        self.__class__ = Subspace
        self.isbasis = isbasis

    def _cohensSimplify(self):
        """
        _cohensSimplify is a common process used in image() and kernel()

        Return a tuple of modified matrix M, image data c and kernel data d.
        """
        M = self.copy()
        c = [0] * (M.row + 1)
        d = [-1] * (M.column + 1)
        for k in range(1, M.column + 1):
            for j in range(1, M.row + 1):
                if not c[j] and M[j, k]:
                    break
            else:           # not found j such that m(j, k)!=0 and c[j]==0
                d[k] = 0
                continue
            top = -ring.inverse(M[j, k])
            M[j, k] = -self.coeff_ring.one
            for s in range(k + 1, M.column + 1):
                M[j, s] = top * M[j, s]
            for i in range(1, M.row + 1):
                if i == j:
                    continue
                top = M[i, k]
                M[i, k] = self.coeff_ring.zero
                for s in range(k + 1, M.column + 1):
                    M[i, s] = M[i, s] + top * M[j, s]
            c[j] = k
            d[k] = j
        return (M, c, d)

    def kernel(self):       # Algorithm 2.3.1 of Cohen's book
        """
        Return a Matrix whose column vectors are one basis of self's kernel,
        or return None if self's kernel is 0.
        """
        tmp = self._cohensSimplify()
        M, d = tmp[0], tmp[2]
        basis = []
        for k in range(1, M.column + 1):
            if d[k]:
                continue
            vect = []
            for i in range(1, M.column + 1):
                if d[i] > 0:
                    vect.append(M[d[i], k])
                elif i == k:
                    vect.append(self.coeff_ring.one)
                else:
                    vect.append(self.coeff_ring.zero)
            basis.append(vect)
        dimension = len(basis)
        if dimension == 0:
            return None
        output = zeroMatrix(self.column, dimension, self.coeff_ring)
        for j in range(1, dimension + 1):
            output.setColumn(j, basis[j - 1])
        return output

    def image(self):        # Algorithm 2.3.2 of Cohen's book
        """
        Return a Matrix which column vectors are one basis of self's image,
        or return None if self's image is 0.
        """
        tmp = self._cohensSimplify()
        M, c = tmp[0], tmp[1]
        basis = []
        for j in range(1, M.row + 1):
            if c[j]:
                basis.append(self[c[j]])
        dimension = len(basis)
        if dimension == 0:
            return None
        output = zeroMatrix(self.row, dimension, self.coeff_ring)
        for j in range(1, dimension + 1):
            output.setColumn(j, basis[j - 1])
        output._selectMatrix()
        return output

    def rank(self):
        """
        Return rank of self.
        """
        img = self.image()
        if img:
            return img.column
        else:
            return 0

    def inverseImage(self, V):    # modified Algorithm 2.3.5 of Cohen's book
        """
        inverseImage(V) -> X
        
        such that self * X == V
        """
        if isinstance(V, vector.Vector):
            if self.row != len(V):
                raise vector.VectorSizeError()
            B = createMatrix(len(V), 1, V.compo)
        else:
            if self.row != V.row:
                raise MatrixSizeError()
            B = V.copy() # step 1
        M = self.copy()
        m = M.row
        n = M.column
        r = B.column
        X = zeroMatrix(n, r, self.coeff_ring)
        non_zero = []
        i = 1
        # step 2
        for j in range(1, n + 1):
            # step 3
            for k in range(i, m + 1):
                if M[k, j]:
                    break
            else:
                continue
            # step 4
            if k > i:
                for l in range(j, n + 1):
                    t = M[i, l]
                    M[i, l] = M[k, l]
                    M[k, l] = t
                B.swapRow(i, k)
            # step 5
            d = ring.inverse(M[i, j])
            for k in range(i + 1, m + 1):
                ck = d * M[k, j]
                for l in range(j + 1, n + 1):
                    M[k, l] = M[k, l] - ck * M[i, l]
                for l in range(r):
                    B[k, l] = B[k, l] - ck * B[i, l]
            non_zero.insert(0, j)
            i += 1
            if i > m:
                break
        # step 6
        i -= 1
        zero = self.coeff_ring.zero
        for j in non_zero:
            d = ring.inverse(M[i, j])
            for k in range(r):
                sums = zero
                for l in range(j + 1, n + 1):
                    sums = sums + M[i, l] * X[l, k]
                X[j, k] = (B[i, k] - sums) * d
            i -= 1
        # step 7
        i = len(non_zero) + 1
        for j in range(1, r + 1):
            for k in range(i, m + 1):
                if B[k, j]:
                    raise NoInverseImage()
        return X

    def solve(self, B):  # modified Algorithm 2.3.4 of Cohen's book
        """
        Return solution X for self * X = B (B is a vector).
        This function returns tuple (V, M) below.
          V: one solution as vector
          M: kernel of self as list of basis vectors.
        If you want only one solution, use 'inverseImage'.

        Warning: B should not be a matrix instead of a vector
        """
        M_1 = self.copy()
        M_1.insertColumn(self.column + 1, B.compo)
        V = M_1.kernel()
        ker = []
        flag = False
        if not V:
            raise NoInverseImage("no solution")
        n = V.row
        for j in range(1, V.column + 1):
            if not bool(V[n, j]): # self's kernel
                ker.append(vector.Vector([V[i, j] for i in range(1, n)]))
            elif not(flag):
                d = -ring.inverse(V[n, j])
                sol = vector.Vector([V[i, j] * d for i in range(1, n)])
                flag = True
        if not(flag):
            raise NoInverseImage("no solution")
        return sol, ker

    def columnEchelonForm(self):  # Algorithm 2.3.11 of Cohen's book
        """
        Return a Matrix in column echelon form whose image is equal to 
        the image of self.
        """
        M = self.copy()
        k = M.column
        for i in range(M.row, 0, -1):
            for j in range(k, 0, -1):
                if M[i, j]:
                    break
            else:
                continue
            d = ring.inverse(M[i, j])
            for l in range(1, i + 1):
                t = d * M[l, j]
                M[l, j] = M[l, k]
                M[l, k] = t
            for j in range(1, M.column + 1):
                if j == k:
                    continue
                for l in range(1, i + 1):
                    M[l, j] = M[l, j] - M[l, k] * M[i, j]
            k -= 1
        return M


class FieldSquareMatrix(RingSquareMatrix, FieldMatrix):
    """
    FieldSquareMatrix is a class for square matrices in field.
    """

    def __init__(self, row, column=0, compo=0, coeff_ring=0):
        """
        FieldSquareMatrix(row [, column, components, coeff_ring])
        FieldSquareMatrix must be row == column .
        """
        self._initialize(row, column, compo, coeff_ring)

    def triangulate(self):
        """
        Return triangulated matrix of self.
        """
        triangle = self.copy()
        flag = False # for calculation of determinant
        for i in range(1, triangle.row + 1):
            if not triangle[i, i]:
                for k in range(i + 1, triangle.row + 1):
                    if triangle[k, i]:
                        triangle.swapRow(i + 1, k + 1)
                        flag = not(flag)
                        break # break the second loop
                else:
                    # the below components are all 0. Back to the first loop
                    continue
            for k in range(i + 1, triangle.row + 1):
                inv_i_i = ring.inverse(triangle[i, i])
                ratio = triangle[k, i] * inv_i_i
                for l in range(i, triangle.column + 1):
                    triangle[k, l] = triangle[k, l] - triangle[i, l] * ratio
        if flag:
            for j in range(triangle.row, triangle.column + 1):
                triangle[triangle.row, j] = -triangle[triangle.row, j]
        return triangle

    def determinant(self):
        """
        Return determinant of self.
        """
        triangle = self.triangulate()
        det = self.coeff_ring.one
        for i in range(1, self.row + 1):
            det = det * triangle[i, i]
        return det

    def inverse(self, V=1): # modified Algorithm 2.2.2, 2.3.5 of Cohen's book
        """
        Return inverse matrix of self or self.inverse() * V.
        If inverse does not exist, raise NoInverse error.
        """
        if isinstance(V, vector.Vector):
            if self.row != len(V):
                raise vector.VectorSizeError()
            B = createMatrix(len(V), 1, V.compo)
        elif isinstance(V, Matrix):
            if self.row != V.row:
                raise MatrixSizeError()
            B = V.copy() # step 1
        else: # V==1
            B = unitMatrix(self.row, self.coeff_ring)
        M = self.copy()
        n = M.row
        r = B.column
        X = zeroMatrix(n, r, self.coeff_ring)
        # step 2
        for j in range(1, n + 1):
            # step 3
            for i in range(j, n + 1):
                if M[i, j]:
                    break
            else:
                raise NoInverse()
            # step 4
            if i > j:
                for l in range(j, n + 1):
                    t = M[i, l]
                    M[i, l] = M[j, l]
                    M[j, l] = t
                B.swapRow(i, j)
            # step 5
            d = ring.inverse(M[j, j])
            for k in range(j + 1, n + 1):
                ck = d * M[k, j]
                for l in range(j + 1, n + 1):
                    M[k, l] = M[k, l] - ck * M[j, l]
                for l in range(1, r + 1):
                    B[k, l] = B[k, l] - ck * B[j, l]
        # step 6
        for i in range(n, 0, -1):
            d = ring.inverse(M[i, i])
            for k in range(1, r + 1):
                sums = self.coeff_ring.zero
                for j in range(i + 1, n + 1):
                    sums = sums + M[i, j] * X[j, k]
                X[i, k] = (B[i, k] - sums) * d
        if r != 1:
            return X
        else:
            return X[1]

    def hessenbergForm(self):      # Algorithm 2.2.9 of Cohen's book
        """Return a Matrix in Hessenberg Form."""
        n = self.row
        zero = self.coeff_ring.zero
        # step 1
        H = self.copy()
        for m in range(2, H.row):
            # step 2
            for i in range(m + 1, n + 1):
                if H[i, m - 1]:
                    break
            else:
                continue
            tinv = ring.inverse(H[i, m - 1])
            if i > m:
                for j in range(m - 1, n + 1):
                    tmp = H[i, j]
                    H[i, j] = H[m, j]
                    H[m, j] = tmp
                H.swapColumn(i, m)
            # step 3
            for i in range(m + 1, n + 1):
                if H[i, m - 1]:
                    u = H[i, m - 1] * tinv
                    for j in range(m, n + 1):
                        H[i, j] = H[i, j] - u * H[m, j]
                    H[i, m - 1] = zero
                    H.setColumn(m, H[m] + u * H[i])
        return H

    def LUDecomposition(self):
        """
        LUDecomposition() -> (L, U)
        
        L and U are matrices such that
            self == L * U
            L : lower triangular matrix
            U : upper triangular matrix
        """

        n = self.row
        L = unitMatrix(n, self.coeff_ring)
        U = self.copy()
        # initialize L and U
        for i in range(1, n + 1):
            for j in range(i + 1, n + 1):
                L[j, i] = U[j, i] * ring.inverse(U[i, i])
                for k in range(i, n + 1):
                    U[j, k] = U[j, k] - U[i, k] * L[j, i]
        return (L, U)


class MatrixRing (ring.Ring):
    """
    MatrixRing is a class for matrix rings.
    """

    _instances = {}

    def __init__(self, size, scalars):
        """
        MatrixRing(size, scalars)
        
        size: size of matrices (positive integer)
        scalars: ring of scalars
        """
        ring.Ring.__init__(self)
        self.size = size
        self.scalars = scalars

    def __eq__(self, other):
        """
        self == other
        """
        return (self.__class__ == other.__class__ and
                self.size == other.size and
                self.scalars == other.scalars)

    def __hash__(self):
        return self.scalars ** self.size

    def __repr__(self):
        return "MatrixRing(%d, %s)" % (self.size, self.scalars)

    def __str__(self):
        return "M_%d(%s)" % (self.size, str(self.scalars))

    @classmethod
    def getInstance(cls, size, scalars):
        """
        Return the cached instance of the specified matrix ring.  If
        the specified ring is not cached, it is created, cached and
        returned.
        
        The method is a class method.
        """
        if (size, scalars) not in cls._instances:
            anInstance = MatrixRing(size, scalars)
            cls._instances[size, scalars] = anInstance
        return cls._instances[size, scalars]

    def unitMatrix(self):
        """
        Return the unit matrix.
        """
        return self.one.copy()

    def _getOne(self):
        """
        getter for one (unit matrix)
        """
        if self._one is None:
            self._one = unitMatrix(self.size, self.scalars)
        return self._one

    one = property(_getOne, None, None, "multiplicative unit")

    def zeroMatrix(self):
        """
        Return the zero matrix.
        """
        return self.zero.copy()

    def _getZero(self):
        """
        Return zero matrix.
        """
        if self._zero is None:
            self._zero = zeroMatrix(self.size, self.scalars)
        return self._zero

    zero = property(_getZero, None, None, "additive unit")

    def createElement(self, compo):
        """
        Return a newly created matrix from 'compo'.

        'compo' must be a list of n*n components in the scalar ring,
        where n = self.size.
        """
        return createMatrix(self.size, compo, self.scalars)

    def getCharacteristic(self):
        """
        Return the characteristic of the ring.
        """
        return self.scalars.getCharacteristic()

    def issubring(self, other):
        """
        Report whether another ring contains the ring as a subring.
        """
        if other is self:
            return True
        if not isinstance(other, MatrixRing):
            return False
        return self.size == other.size and self.scalars.issubring(other.scalars)

    def issuperring(self, other):
        """
        Report whether the ring is a superring of another ring.
        """
        if other is self:
            return True
        if not isinstance(other, MatrixRing):
            return False
        return self.size == other.size and \
        self.scalars.issuperring(other.scalars)

    def getCommonSuperring(self, other):
        """
        Return common super ring of self and another ring.
        """
        if not isinstance(other, MatrixRing) or self.size != other.size:
            raise TypeError("no common super ring")
        return MatrixRing.getInstance(self.size, 
        self.scalars.getCommonSuperring(other.scalars))


class Subspace(FieldMatrix):
    """
    Subspace is a class for subspaces.
    """

    def __init__(self, row, column, compo=0, coeff_ring=0, isbasis=None):
        """
        Subspace(row, column [,components, coeff_ring, isbasis])
        """
        if isinstance(compo, bool):
            isbasis = compo
            compo = 0
        elif isinstance(coeff_ring, bool):
            isbasis = coeff_ring
            coeff_ring = 0
        self._initialize(row, column, compo, coeff_ring)
        self.isbasis = isbasis

    @classmethod
    def fromMatrix(cls, mat, isbasis=None):
        """
        A constructor class method, which creates Subspace from a
        Matrix instance.
        """
        compo = []
        for row in mat.compo:
            compo += row
        return cls(mat.row, mat.column, compo, mat.coeff_ring, isbasis)

    def toFieldMatrix(self):
        """
        Subspace -> Field(Square)Matrix
        """
        if self.row == self.column:
            self.__class__ = FieldSquareMatrix
        else:
            self.__class__ = FieldMatrix

    def isSubspace(self, other):
        """ 
        Check self is in other as subspace
        """
        try:
            other.inverseImage(self)
            return True
        except:
            return False

    def toBasis(self):
        """
        Change matrix to basis.
        """
        if not self.isbasis:
            basis = self.image()
            if not basis: # zero space
                basis = zeroMatrix(self.row, 1, self.coeff_ring)
            self.compo = basis.compo
            self.column = basis.column
            self.isbasis = True

    def supplementBasis(self):     # Modified Algorithm 2.3.6 of Cohen's book
        """
        Return a basis of full space, which including self's column vectors.
        """
        self.toBasis()
        if self.row < self.column:
            raise MatrixSizeError()
        n = self.row
        k = self.column
        M = self.copy()
        pnt = range(1, self.row + 1)
        for s in range(1, k + 1):
            for t in range(s, n + 1):
                if M[t, s]:
                    break
            else: # zero space
                return unitMatrix(self.row, self.coeff_ring)
            d = ring.inverse(M[t, s])
            M[t, s] = 1
            if t != s:
                pnt[t - 1] = pnt[s - 1]
            for j in range(s + 1, k + 1):
                if t != s:
                    tmp = M[s, j]
                    M[s, j] = M[t, j]
                    M[t, j] = tmp
                M[s, j] *= d
                for i in range(1, n + 1):
                    if i != s and i != t:
                        M[i, j] = M[i, j] - M[i, s] * M[s, j]
        B = self.copy()
        one = self.coeff_ring.one
        zeros = [self.coeff_ring.zero] * n
        for i in pnt[k: ]:
            e_i = zeros
            e_i[i - 1] = one
            B.extendColumn(e_i)
        return Subspace.fromMatrix(B, True)

    def sumOfSubspaces(self, other): # Algorithm 2.3.8 of Cohen's book
        """
        Return space which is sum of self and other.
        """
        if self.row != other.row:
            raise MatrixSizeError()
        N = self.copy()
        N.extendColumn(other)
        return Subspace.fromMatrix(N.image(), True)

    def intersectionOfSubspaces(self, other): # Algorithm 2.3.9 of Cohen's book
        """
        Return space which is intersection of self and other.
        """
        if self.row != other.row:
            raise MatrixSizeError()
        M1 = self.copy()
        M1.extendColumn(other)
        N = M1.kernel()
        if not N:
            zeroMatrix(self.row, 1, self.coeff_ring)
        N1 = N.getBlock(1, 1, self.column, N.column) # N.column is dim(ker(M1))
        return Subspace.fromMatrix((self * N1).image(), True)


# --------------------------------------------------------------------
#       the belows are not class methods
# --------------------------------------------------------------------

def createMatrix(row, column=0, compo=0, coeff_ring=0):
    """
    generate new Matrix or SquareMatrix class.
    """
    if isinstance(compo, (ring.Ring, int)):
        coeff_ring = compo
        compo = 0
    if isinstance(column, list):
        compo = column
        column = row
    elif isinstance(column, ring.Ring):
        coeff_ring = column
        column = row
    if coeff_ring != 0 and isinstance(coeff_ring, int):
        coeff_ring = _getRing(coeff_ring)
    if compo == 0:
        return zeroMatrix(row, column, coeff_ring)
    if coeff_ring == 0:
        if isinstance(compo[0], (list, tuple)):
            coeff_ring = ring.getRing(compo[0][0])
        elif isinstance(compo[0], vector.Vector):
            coeff_ring = ring.getRing(compo[0][1])
        else:
            coeff_ring = ring.getRing(compo[0])
    if coeff_ring.isfield():
        if row == column:
            return FieldSquareMatrix(row, compo, coeff_ring)
        else:
            return FieldMatrix(row, column, compo, coeff_ring)
    else:
        if row == column:
            return RingSquareMatrix(row, compo, coeff_ring)
        else:
            return RingMatrix(row, column, compo, coeff_ring)

def unitMatrix(size, coeff=1):
    """
    return unit matrix of size.
    coeff is subclass for ring.Ring or ring.Ring.one.
    """
    if isinstance(coeff, ring.Ring):
        one = coeff.one
        zero = coeff.zero
    else:
        one = coeff
        coeff = ring.getRing(one)
        zero = coeff.zero
    unit_matrix = [one]
    units = [zero] * size + [one]
    for i in range(size - 1):
        unit_matrix = unit_matrix + units
    return createMatrix(size, size, unit_matrix, coeff)

identityMatrix = unitMatrix

def zeroMatrix(row, column=None, coeff=0):
    """
    return zero matrix.
    coeff is subclass for ring.Ring or ring.Ring.zero.
    """
    if column == 0:
        coeff = 0
        column = row
    if not(isinstance(column, (int, long))):
        if column == None:
            column = row
        else:
            coeff = column
            column = row
    if isinstance(coeff, ring.Ring):
        zero = coeff.zero
    else:
        zero = coeff
        coeff = ring.getRing(coeff)
    zero_matrix = [zero] * (row * column)
    return createMatrix(row, column, zero_matrix, coeff)


#--------------------------------------------------------------------
#   define exceptions
#--------------------------------------------------------------------

class MatrixSizeError(Exception):
    """Invalid input error for matrix size."""
    pass

class VectorsNotIndependent(Exception):
    """Invalid input error because column vectors are linear dependent."""
    pass

class NoInverseImage(Exception):
    """Invalid input error because self do not have inverse image for input."""
    pass

class NoInverse(Exception):
    """Invalid input error because matrix is not invertible."""
    pass

########## for utility
def _getRing(coeff_ring):
    """
    If 'coeff_ring' represents characteristic, return F_p or Z_n instance.
    """
    if isinstance(coeff_ring, int):
        try:
            import nzmath.finitefield
            coeff_ring = nzmath.finitefield.FinitePrimeField(coeff_ring)
        except:
            import nzmath.intresidue
            coeff_ring = nzmath.intresidue.IntegerResidueClassRing(coeff_ring)
    return coeff_ring
