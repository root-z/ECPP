"""
Group Theorical module
"""
import math
import nzmath.rational as rational
import nzmath.factor.misc as factor_misc
import nzmath.vector as vector
import nzmath.matrix as matrix
import nzmath.compatibility


class Group:
    """
    This is a class for finite group.
    """

    def __init__(self, value, operation=-1):
        if isinstance(value, Group):
            self.entity = value.entity
            if operation == -1:
                self.operation = value.operation
            else:
                self.setOperation(operation)
        else:
            self.entity = value
            if operation == -1:
                self.setOperation(0)
            else:
                self.setOperation(operation)

    def __repr__(self):
        if hasattr(self.entity, "__repr__"):
            return self.entity.__repr__()
        else:
            return repr(self.entity.__class__.__name__)

    def __str__(self):
        if hasattr(self.entity, "__str__"):
            return self.entity.__str__()
        else:
            return str(self.entity.__class__.__name__)

    def __eq__(self, other):
        if self.entity == other.entity and self.operation == other.operation:
            return True
        else:
            return False

    def __ne__(self, other):
        return not (self == other)

    def setOperation(self, value):
        """
        Change group type for additive(0) or multiplicative(1).
        """
        if isinstance(value, int) :
            self.operation = (value & 1)
        else:
            raise TypeError("invalid input")

    def createElement(self, *value):
        """
        Create group element with value.
        Return GroupElement instance.
        """
        return GroupElement(self.entity.createElement(*value), self.operation)

    def identity(self):
        """
        Return identity element(unit).
        Return addtive 0 or multiplicative 1.
        """
        if hasattr(self.entity, "identity"):
            return GroupElement(self.entity.identity(), self.operation)
        else:
            if self.operation:
                return GroupElement(self.entity.one, 1)
            else:
                return GroupElement(self.entity.zero, 0)

    def grouporder(self):
        """
        Return group order(Cardinality).
        """
        order = 0
        if hasattr(self.entity, "grouporder"):
            order = self.entity.grouporder()
        elif hasattr(self.entity, "card"):
            order = card(self.entity)
        else:
            order = self.entity.m
        if self.operation and hasattr(self.entity, "zero"): # *-group
            order -= 1
        return order


class GroupElement:
    """
    This is a class for finite group element.
    """

    def __init__(self, value, operation=-1):
        self.operation = operation
        if isinstance(value, GroupElement):
            self.entity = value.entity
            self.set = value.set
            self._set_name = value._set_name
            if operation == -1:
                self.operation = value.operation
            else:
                self.setOperation(operation)
        else:
            self.entity = value
            if operation == -1:
                if self._type_check(1):
                    self.operation = 1
                if self._type_check(0):
                    self.operation = 0
                if self.operation == -1:
                    raise TypeError("This element isn't Group")
            self.set = self.getGroup()
            self.setOperation(self.operation) # mainly operation
            self._set_name = self.set.entity.__class__.__name__

    def __repr__(self):
        return self._set_name + ',' + repr(self.entity)

    def __str__(self):
        return self._set_name + ',' + str(self.entity)

    def __eq__(self, other):
        if self.entity == other.entity and self.operation == other.operation:
            return True
        else:
            return False

    def __ne__(self, other):
        return not (self == other)

    def _type_check(self, value):
        """
        Check group type is value or not.
        """
        a = self.entity
        if not (value & 1):
            if hasattr(a, "__add__") and hasattr(a, "__mul__"):
                return True
            else:
                return False
        else:
            if hasattr(a, "__mul__") and hasattr(a, "__pow__"):
                return True
            else:
                return False

    def setOperation(self, value):
        """
        Change group type for additive(0) or multiplicative(1).
        """
        if isinstance(value, int) and self._type_check(value):
            self.operation = (value & 1)
        else:
            raise TypeError("invalid input")
        self.set.setOperation(self.operation)

    def ope(self, other):
        """
        Group basic operation.
        """
        if  not self.operation:
            return GroupElement(self.entity + other.entity, self.operation)
        else:
            return GroupElement(self.entity * other.entity, self.operation)

    def ope2(self, other):
        """
        Group extended operation
        """
        if rational.isIntegerObject(other):
            if not self.operation:
                return GroupElement(self.entity * other, self.operation)
            else:
                return GroupElement(self.entity ** other, self.operation)
        else:
            raise TypeError("input integer")

    def inverse(self):
        """
        Return inverse element.
        """
        ele = self.entity
        cla = self.set
        operation = self.operation
        if hasattr(cla.entity, "zero") and (ele == cla.entity.zero):
            return self
        else:
            if not operation:
                return GroupElement(-ele, operation)
            elif hasattr(ele, "inverse"):
                return GroupElement(ele.inverse(), operation)
            else:
                return GroupElement(self.ope2(cla.order() - 1), operation)

    def order(self):
        """
        Compute order using grouporder factorization.
        """
        clas = self.set.entity
        if hasattr(clas, "zero") and self.entity == clas.zero:
            return 1
        ordfact = factor_misc.FactoredInteger(self.set.grouporder())
        identity = self.set.identity()
        k = 1
        for p, e in ordfact:
            b = self.ope2(int(ordfact) // (p ** e))
            while b != identity:
                k = k * p
                b = b.ope2(p)
        return k

    def t_order(self, v=2):
        """
        Compute order using Terr's Baby-step Giant-step algorithm.
        """
        if (v < 1) or not(rational.isIntegerObject(v)):
            raise TypeError("input integer v >= 1")
        e = self.set.identity()
        a = self.set.identity()
        R = [(e, 0)]
        for i in range(1, v + 1):
            a = self.ope(a)
            if a == e:
                return i
            else:
                R.append((a, i))
        j = 0
        b = a.ope2(2)
        t = 2 * v
        while True:
            for (c, k) in R:
                if b == c:
                    return (t - k)
            a = self.ope(a)
            j += 1
            R.append((a, j + v))
            b = a.ope(b)
            t += j + v

    def getGroup(self):
        """
        Return the group which element belongs.
        """
        if self._type_check(0) and self._type_check(1):
            if hasattr(self.entity, "getRing"):
                return Group(self.entity.getRing(), self.operation)
            else:
                return Group(self.entity, self.operation)
        else:
            if hasattr(self.entity, "getGroup"):
                return Group(self.entity.getGroup(), self.operation)
            else:
                return Group(self.entity, self.operation)


class GenerateGroup(Group):
    """
    This is a class for finite group with generator.
    """

    def __init__(self, value, operation=-1, entity=None):
        if isinstance(value, list):
            temp = value
        else:
            TypeError("invalid input")
        if entity:
            self.entity = entity
        else:
            self.entity = GroupElement(value[0]).set.entity
        self.generator = []
        for a in temp:
            self.generator.append(GroupElement(a))
        if operation == -1:
            self.operation = self.generator[0].operation
        else:
            self.setOperation(operation)

    def isGroupElement(self, other):
        """
        Check whether other is self group element or not.
        """
        pass

    def setOperation(self, value):
        """
        Change group type for additive(0) or multiplicative(1).
        """
        if isinstance(value, int) :
            self.operation = (value & 1)
        else:
            raise TypeError("invalid input")
        for a in self.generator:
            a.setOperation(value)


class AbelianGenerate(GenerateGroup):
    """
    This is a class for finite abelian group with genarator.
    """

    def relationLattice(self):
        """
        Return relation lattice basis as column vector matrix for generator.
        If B[j]=transpose(b[1,j],b[2,j],..,b[l,j]),
        it satisfies that product(generator[i]**b[i,j])=1 for each j.
        """
        if not hasattr(self, "relation"):
            l = len(self.generator)
            b = matrix.RingSquareMatrix(l)
            H1 = [(self.identity(), vector.Vector([0] * l))]
            H2 = list(H1)
            m = 1
            a_baby_s, giant_s = list(H1), list(H2)
            pro_I1, pro_I2, pro_diag = 1, 1, 1
            e_vec = []
            g_gen = []
            for i in range(1, l + 1):
                e_i = vector.Vector([0] * l)
                e_i[i] = 1
                e_vec.append(e_i)
                g_gen.append(self.generator[i - 1])
            for j in range(1, l + 1):
                e = 1
                baby_s = list(a_baby_s)
                baby_e = GroupElement(g_gen[j - 1])
                giant_e = GroupElement(g_gen[j - 1])
                flag = False
                while not flag:
                    for (g_g, v) in giant_s:
                        for (g_b, w) in baby_s:
                            if g_g.ope(giant_e) == g_b:
                                b[j] = v + w + (e * (e + 1) >> 1) * e_vec[j - 1]
                                flag = True
                                break
                        if flag:
                            break
                    if flag:
                        break
                    for (g_a, v_a) in a_baby_s:
                        baby_s.append((g_a.ope(baby_e), v_a - e * e_vec[j - 1]))
                    e += 1
                    baby_e = baby_e.ope(g_gen[j - 1])
                    giant_e = giant_e.ope(baby_e)
                if (j < l) and (b[j, j] > 1):
                    pro_diag *= b[j, j]
                    pro_diag_root = math.sqrt(pro_diag)
                    if (b[j, j] * pro_I1) <= pro_diag_root or j == 1:
                        temp = list(H1)
                        for (g, v) in temp:
                            g_j_inv = g_gen[j - 1].inverse()
                            H1_1 = GroupElement(g)
                            for x in range(1, b[j, j]):
                                H1_1 = H1_1.ope(g_j_inv)
                                H1.append((H1_1, v + x * e_vec[j - 1]))
                        pro_I1 *= b[j, j]
                    else:
                        if m > 1:
                            temp = list(H2)
                            for (g, v) in temp:
                                H2_1 = GroupElement(g)
                                for x in range(1, b[m, m]):
                                    H2_1 = H2_1.ope(g_gen[m - 1])
                                    H2.append((H2_1, v + x * e_vec[m - 1]))
                            pro_I2 *= b[m, m]
                        m = j
                    s = int(math.ceil(pro_diag_root / pro_I1))
                    if len(H2) > 1:
                        t = int(math.ceil(pro_diag_root / pro_I2))
                    else:
                        t = 1
                    a_baby_s, giant_s = list(H1), list(H2)
                    g_m_inv = g_gen[m - 1].inverse()
                    for (h1, v) in H1:
                        H1_1 = GroupElement(h1)
                        for r in range(1, s):
                            H1_1 = H1_1.ope(g_m_inv)
                            a_baby_s.append((H1_1, v + r * e_vec[m - 1]))
                    g_m_s = g_gen[m - 1].ope2(s)
                    for (h2, v) in H2:
                        H2_1 = GroupElement(h2)
                        for q in range(1, t):
                            H2_1 = H2_1.ope(g_m_s)
                            giant_s.append((H2_1, v + (q * s) * e_vec[m - 1]))
            self.relation = b
        return self.relation

    def computeStructure(self):
        """
        Compute Finite Abelian Group Structure.
        """
        B = self.relationLattice()
        U_d, V, M = B.extsmithNormalForm()
        det = M.determinant()
        U_d.toFieldMatrix()
        U = U_d.inverse()
        for i in range(U_d.row):
            U_d[i] = (U_d[i] % det)
        structure = []
        l = M.row
        for j in range(1, l):
            if M[j, j] != 1 or j == 1:
                g = self.identity()
                for i in range(1, l+1):
                    g = g.ope(self.generator[i-1].ope2(int(U[i, j])))
                structure.append((g, M[j, j]))
            else:
                break
        return structure, det
