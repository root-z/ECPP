import random
import nzmath.combinatorial as combinatorial
import nzmath.gcd as gcd
import nzmath.matrix as matrix
import nzmath.rational as rational

class Permute:

    """
    This is a class for 'normal' type element in the permutation group.
    Example, [2,3,1,5,4]
    This means [1 2 3 4 5]
               [2 3 1 5 4]
    (It is 1:1 onto mapping, 1->2, 2->3, 3->1, 4->5, 5->4)
    """

    def __init__(self, value, key=None, flag=False):
        """
        You can initialize with various one-to-one onto mapping.
        Example,
        Permute([2,3,4,5,1]) -> normal type
        Permute([3,4,2,1,0], 0) -> [4,5,3,2,1]-normal type(index start with 0)
        Permute(['b','c','d','e','a'], 1) -> [2,3,4,5,1]-normal type(identity is ascending order)
        Permute(['b','c','d','e','a'], -1) -> [4,3,2,1,5]-normal type(identity is descending order)
        Permute(['b','c','d','e','a'], ['b','a', 'c','d','e']) -> [1,3,4,5,2]-normal type(identity=key)
        Permute({'a':'b','b':'c','c':'d','d':'e','e':'a'}) -> [2,3,4,5,1]-normal type
        """
        if isinstance(value, dict):
            if key:
                raise TypeError("cannot convert Permute. I think `key` should be None.")
            data = value.values()
            key = value.keys()
        elif isinstance(value, (list, tuple)):
            data = list(value)
        else:
            raise TypeError("cannot convert Permute. `value` should be a list or a dict.")
        if key == 0:
            self.data = [i + 1 for i in data]
            self.key = range(len(data))
        elif key:
            if isinstance(key, (list, tuple)):
                self.key = list(key)
                if len(value) != len(key):
                    raise TypeError("cannot convert Permute. The length of `key` should be equal to that of `value`.")
            elif key == 1:
                p_key = list(data)
                p_key.sort()
                self.key = p_key
            elif key == -1:
                p_key = list(data)
                p_key.sort()
                p_key.reverse()
                self.key = p_key
            else:
                raise TypeError("cannot convert Permute. `key` should be a list.")
            key = self.key
            if flag:
                self.data = data
            else:
                self.data = [key.index(x) + 1 for x in data]
        else:
            self.data = data
            self.key = range(1, len(data) + 1)
        data = self.data
        p_data = range(len(data))
        for x in data:
            if not rational.isIntegerObject(x):
                raise TypeError("cannot convert Permute. `flag` should be False.")
            elif x <= 0 or x > len(data):
                raise TypeError("cannot convert Permute. The map should be onto.")
            elif p_data[x-1] == -1:
                raise ValueError("cannot convert Permute. The map should be one-to-one.")
            else:
                p_data[x-1] = -1

    def __getitem__(self, other):
        try:
            idx = self.key.index(other)
        except ValueError:
            raise ValueError("The indices must be elements of self.key.")
        return self.key[self.data[idx]  - 1]

    def __mul__(self, other):
        """
        Compute the multiplication, that is, the composite of two maps self \circ other.
        The map is the application of `self` to the result of `other`.
        """
        s_data = self.data
        o_data = other.data
        lth = len(s_data)
        if self.key != other.key or lth != len(o_data):
            raise TypeError("cannot multiply a Permute by a Permute which has a different type.")
        sol = [s_data[o_data[i] - 1] for i in range(lth)]
        return Permute(sol, self.key, flag=True)

    def __rmul__(self, other):
        return other * self

    def __div__(self, other):
        return self * (other.inverse())

    __truediv__ = __div__

    def __rdiv__(self, other):
        return other * (self.inverse())

    def __pow__(self, other):
        sol = self.__class__(self.data, self.key, flag=True)
        if not rational.isIntegerObject(other):
            raise TypeError("cannot pow operation with %s" % other)
        if other > 0:
            for i in range(other - 1):
                sol = self * sol
        else:
            inv = self.inverse()
            for i in range(abs(other) + 1):
                sol = inv * sol
        return sol

    def __call__(self, other):
        return self.permute(other)

    def setKey(self, key=None):
        """
        Set other key.
        The function may be used if you want to permute a different sequence with
        the same permutation.
        """
        if key == 0:
            self.key = range(len(data))
        elif key:
            if len(key) != len(self.key):
                raise TypeError, "The length of `key` should be equal to that of self.key."
            else:
                if key[0] in self.key: # key transformation
                    data = list(self.data)
                    keys = self.key
                    sol = [0] * len(data)
                    try:
                        for i in range(len(data)):
                            sol[key.index(keys[i])] = key.index(keys[data[i]-1]) + 1
                        self.data = sol
                    except ValueError:
                        pass
                self.key = key
        else:
            self.key = range(1, len(data) + 1)

    def getValue(self):
        """
        Return value of self.
        """
        return [self.key[self.data[i] - 1] for i in range(len(self.data))]

    def inverse(self):
        s_data = self.data
        sol = [0] * len(s_data)
        for i in range(len(s_data)):
            sol[s_data[i] - 1] = i+1
        return Permute(sol, self.key, flag=True)

    def getGroup(self):
        return PermGroup(self.key)

    def numbering(self):
        """
        Return the number of self by the following numbering.

        This is the inductive definition on the dimension.
        It is symmetrical arranging.

        Example,
        2-dimension [1,2], [2,1]
        3-dimension [1,2,3], [2,1,3], [1,3,2], [2,3,1], [3,1,2], [3,2,1]
        4-dimension [1,2,3,4], [2,1,3,4], [1,3,2,4], [2,3,1,4], [3,1,2,4],
                    [3,2,1,4], ..., [4,3,2,1]
        """
        s_data = self.data
        val = [-1] * (len(s_data))
        for i in range(len(s_data)):
            val[s_data[i] - 1] = 0
            for j in range(s_data[i], len(val)):
                if val[j] != -1:
                    val[j] += 1
        sol = 0
        val[0] = 1
        for j in range(len(val) - 1, -1, -1):
            sol = (j+1) * sol + val[j]
        return sol

    def order(self):
        """
        This method returns order for permutation element.
        """
        return self.ToCyclic().order()

    def ToTranspose(self):
        """
        This method returns
         2-dimensional cyclic type element of permutation group.

        The method uses recursion.
        """
        s_data = list(self.data)
        lth = len(s_data)
        if lth == 1:
            return ExPermute(1, [])
        else:
            sol = []
            if s_data[lth - 1] != lth:
                sol.append((s_data[lth - 1], lth))
                s_data[s_data.index(lth)] = s_data[lth - 1]
            sol.extend((Permute(s_data[:lth - 1]).ToTranspose()).data)
            return ExPermute(lth, sol, self.key, flag=True)

    def ToCyclic(self):
        """
        This method returns cyclic type element of permutation group.
        """
        s_data = self.data
        box = list(self.data)
        sol = []
        for i in range(len(s_data)):
            if box[i] != '*':
                p_sol = [(i+1)]
                box[i] = '*'
                j = i
                while s_data[j] != i+1:
                    p_sol.append(s_data[j])
                    j = s_data[j] - 1
                    box[j] = '*'
                if len(p_sol) != 1:
                    sol.append(tuple(p_sol))
        return ExPermute(len(s_data), sol, self.key, flag=True)

    def sgn(self):
        """
        This method returns sign.
        
        If self is even permutation, that is, self can be written as a composition
        of an even number of transpositions, it returns 1. Otherwise,that is, for odd
        permutation, it returns -1.
        """
        return self.ToCyclic().sgn()

    def types(self):
        """
        This method returns 'type' defined by each cyclic element length.
        """
        c_data = self.ToCyclic().data
        sol = [len(c_data[i]) for i in range(len(c_data))]
        sol.sort()
        return sol

    def ToMatrix(self):
        """
        This method returns the permutation matrix.
        """
        lth = len(self.data)
        A = matrix.SquareMatrix(lth)
        for j in range(lth):
            A[j+1, self.data[j]] = 1
        return A

    def permute(self, lists):
        """
        permute list following with self permutation.

        Warning: The method do not check the compatibility of `lists` and self.key (except dict type).
        """
        if len(lists) != len(self.data):
            raise TypeError, "The length of `lists` should be equal to that of self.key."
        if isinstance(lists, dict):
            sol = {}
            key = self.key
            for x in lists.keys():
                sol[key[self.data[key.index(x)] - 1]] = lists[x]
        elif isinstance(lists, (list, tuple)):
            sol = [0] * len(lists)
            for i in range(len(lists)):
                sol[self.data[i] - 1] = lists[i]
        return sol

    def __eq__(self, other):
        s_data = self.data
        o_data = other.data
        lth = len(s_data)
        if self.key != other.key or lth != len(o_data):
            return False
        for i in range(lth):
            if s_data[i] != o_data[i]:
                return False
        return True

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return repr(self.key)+" -> "+repr(self.getValue())

    def __str__(self):
        return str(self.key)+" -> "+str(self.getValue())


class ExPermute:

    """
    This is a class for cyclic type element of permutation group.
    Example, (5, [(1, 2), (3, 4)])
    This means (1, 2)(3, 4)=[2, 1, 4, 3, 5]
    """

    def __init__(self, dim, value, key=None, flag=False):
        if not (rational.isIntegerObject(dim) and isinstance(value, list)):
            raise TypeError("cannot convert ExPermute. `dim` should be an integer and `value` should be a list.")
        self.dim = dim
        data = value
        self.data = []
        if key == 0:
            self.key = range(dim)
            for x in data:
                ele = [ y + 1 for y in x ]
                self.data.append(tuple(ele))
        elif key:
            if isinstance(key, (list, tuple)):
                self.key = list(key)
                if dim != len(key):
                    raise TypeError("cannot convert ExPermute. The length of `key` should be equal to dim.")
            else:
                raise TypeError("cannot convert ExPermute. `key` should be a list or a tuple.")
            key = self.key
            if flag:
                self.data = data
            else:
                for x in data:
                    ele = [key.index(x[i]) + 1 for i in range(len(x))]
                    self.data.append(tuple(ele))
        else:
            self.data = data
            self.key = range(1, dim + 1)
        data = self.data
        for x in data:
            if not isinstance(x, tuple):
                raise TypeError("cannot convert ExPermute. `flag` should be False.")
            box = range(dim)
            for y in x:
                if (y > dim) or (y <= 0):
                    raise TypeError("cannot convert ExPermute. The map should be onto.")
                elif box[y-1] == -1:
                    raise ValueError("cannot convert ExPermute. The map should be one-to-one.")
                else:
                    box[y-1] = -1

    def __getitem__(self, other):
        try:
            idx = self.key.index(other)
        except ValueError:
            raise ValueError("The indices must be elements of self.key.")
        val = idx + 1
        for i in range(len(self.data) - 1, -1, -1):
            data_i = list(self.data[i])
            try:
                pla = data_i.index(val)
                val = data_i[(pla+1) % len(data_i)]
            except ValueError:
                pass
        return self.key[val - 1]

    def __mul__(self, other):
        if self.key != other.key or self.dim != other.dim:
            raise TypeError("cannot multiply an ExPermute by an ExPermute which has a different type.")
        sol = [x for x in self.data] + [x for x in other.data]
        return ExPermute(self.dim, sol, self.key, flag=True)

    def __rmul__(self, other):
        return other * self

    def __div__(self, other):
        return self * other.inverse()

    __truediv__ = __div__

    def __rdiv__(self, other):
        return other * self.inverse()

    def __pow__(self, other):
        sol = ExPermute(self.dim, self.data, self.key, flag=True)  # other instance
        if not rational.isIntegerObject(other):
            raise TypeError("cannot pow operation with %s" % other)
        if other > 0:
            for i in range(other - 1):
                sol = self * sol
        else:
            inv = self.inverse()
            for i in range(abs(other) + 1):
                sol = inv * sol
        return sol

    def __call__(self, other):
        return self.permute(other)

    def setKey(self, key=None):
        """
        Set other key.
        The function may be used if you want to permute a different sequence with 
        the same permutation.
        """
        if key == 0:
            self.key = range(self.dim)
        elif key:
            if len(key) != self.dim:
                raise TypeError, "The lenght of `key` should be equal to that of self.key."
            else:
                if key[0] in self.key: # key transformation
                    data = list(self.data)
                    keys = self.key
                    sol = []
                    try:
                        for x in data:
                            p_ele = []
                            for i in range(len(x)):
                                p_ele.append(key.index(keys[x[i]-1])+1)
                            sol.append(tuple(p_ele))
                        self.data = sol
                    except ValueError:
                        pass
                self.key = key
        else:
            self.key = range(1, self.dim + 1)

    def getValue(self):
        """
        Return value of self.
        """
        out = []
        for x in self.data:
            out.append(tuple([self.key[x[i] - 1] for i in range(len(x))]))
        return out

    def inverse(self):
        s_data = list(self.data)
        s_data.reverse()
        for i in range(len(s_data)):
            ele_data = list(s_data[i])
            if len(s_data[i]) > 2:
                ele_data.reverse()
            s_data[i] = tuple(ele_data)
        return ExPermute(self.dim, s_data, self.key, flag=True)

    def getGroup(self):
        return PermGroup(self.key)

    def order(self):
        """
        This method returns order for permutation element.
        """
        data = self.simplify().data
        sol = 1
        for x in data:
            sol = gcd.lcm(sol, len(x))
        return sol

    def ToNormal(self):
        """
        This method returns normal type element of permutation group.
        """
        dim = self.dim
        s_data = list(self.data)
        s_data.reverse()
        sol = ['*'] * dim
        for x in s_data:
            ele_data = list(x)
            ele_data.append(ele_data[0])
            trans_data = []
            for y in x:
                if sol[y-1] != '*':
                    trans_data.append(sol.index(y))
                else:
                    trans_data.append(y-1)
            for j in range(len(trans_data)):
                sol[trans_data[j]] = ele_data[j+1]
                if sol[trans_data[j]] == trans_data[j] + 1:
                    sol[trans_data[j]] = '*'
        for i in range(dim):
            if sol[i] == '*':
                sol[i] = i+1
        return Permute(sol, self.key, flag=True)

    def simplify(self):
        """
        This method returns more simple element.
        """
        return self.ToNormal().ToCyclic()

    def sgn(self):
        """
        This method returns sign for permutation element.

        If self is even permutation, that is, self can be written as a composition
        of an even number of transpositions, it returns 1. Otherwise,that is, for odd
        permutation, it returns -1.
        """
        sol = 1
        for x in self.data:
            if len(x) & 1 == 0:
                sol = -sol
        return sol

    def permute(self, lists):
        """
        permute list following with self permutation
        
        Warning: The method do not check the compatibility of `lists` and self.key (except dict type).
        """
        if len(lists) != self.dim:
            raise TypeError, "The length of `lists` should be equal to self.dim."
        if isinstance(lists, dict):
            sol = dict(lists)
            key = self.key
            data = self.data
            for i in range(len(data) - 1, -1, -1):
                data_i = data[i]
                first = key[data_i[0] - 1]
                for j in range(len(data_i) - 1):
                    idx = key[data_i[j+1] - 1]
                    sol[first], sol[idx] = sol[idx], sol[first]
        elif isinstance(lists, (list, tuple)):
            sol = list(lists)
            data = self.data
            for i in range(len(data) - 1, -1, -1):
                data_i = data[i]
                first = data_i[0] - 1
                for j in range(len(data[i]) - 1):
                    idx = data_i[j+1] - 1
                    sol[first], sol[idx] = sol[idx], sol[first]
        return sol

    def __eq__(self, other):
        if self.key != other.key or self.dim != other.dim:
            return False
        s_data = (self.simplify()).data
        o_data = (other.simplify()).data
        if len(s_data) != len(o_data):
            return False
        for i in range(len(s_data)):
            for j in range(len(s_data[i])):
                if s_data[i][j] != o_data[i][j]:
                    return False
        return True

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return repr(self.getValue()) + " <" + repr(self.key) + ">"

    def __str__(self):
        self.data = self.simplify().data
        return str(self.getValue()) + " <" + str(self.key) + ">"

class PermGroup:
    """
    This is a class for permutation group.
    """
    def __init__(self, key):
        if isinstance(key, (int, long)):
            self.key = range(1, key + 1)
        elif isinstance(key, (list, tuple)):
            self.key = list(key)
        elif isinstance(key, dict):
            self.key = dict.keys()
        else:
            raise TypeError, "cannot convert PermGroup. `key` should be an integer or a list/tuple/dict."

    def __repr__(self):
        return repr(self.key)

    def __str__(self):
        return str(self.key)

    def __eq__(self, other):
        if self.key == other.key:
            return True
        else:
            return False

    def __ne__(self, other):
        return not(self == other)

    def card(self):
        return self.grouporder()

    def createElement(self, seed):
        """
        Create Permute or ExPermute with seed.
        createElement(dict) -> Permute(dict)
        createElement(tuple) -> Permute(list(tuple), self.key)
        createElement(key_element_list) -> Permute(key_element_list, self.key)
        createElement([cyclic_tuple,...]) -> ExPermute(len(self.key), [cyclic_tuple,...], self.key)
        """
        if isinstance(seed, dict):
            if set(self.key) == set(dict.keys()):
                return Permute(seed)
            else:
                raise TypeError, "`seed`.key should be equal to self.key."
        elif isinstance(seed, tuple):
            return Permute(list(seed). self.key)
        elif isinstance(seed, list):
            if seed[0] in self.key:
                return Permute(seed, self.key)
            elif isinstance(seed[0], tuple):
                return ExPermute(len(self.key), seed, self.key)
        raise TypeError, "`seed` should be a dict/tuple/list."

    def identity(self):
        return Permute(self.key, self.key)

    def identity_c(self):
        """
        Return identity for cyclic type.
        """
        return ExPermute(len(self.key), [], self.key)

    def grouporder(self):
        return combinatorial.factorial(len(self.key))

    def randElement(self):
        """
        Create random Permute type element.
        """
        copy = list(self.key)
        sol = []
        while copy:
            sol.append(copy.pop(random.randrange(len(copy))))
        return Permute(sol, self.key)
