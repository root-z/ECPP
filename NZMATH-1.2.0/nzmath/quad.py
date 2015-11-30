from __future__ import division
import math
import copy
import warnings
import nzmath.gcd as gcd
import nzmath.arith1 as arith1
import nzmath.prime as prime
import nzmath.factor.misc as misc
import nzmath.bigrandom as bigrandom


class ReducedQuadraticForm(object):
    """
    The class is for reduced quadratic form.
    """
    def __init__(self, element, unit):
        self.element = element # form = (a_1, a_2, a_3)
        self.unit = unit

    def __repr__(self):
        return_str = 'ReducedQuadraticForm(%d, %d, %d)' % self.element
        return return_str

    def __str__(self):
        return_str = '[%d, %d, %d]' % self.element
        return return_str

    def __mul__(self, other):
        if not isinstance(other, ReducedQuadraticForm):
            return NotImplemented
        composition = compositePDF(self.element, other.element)
        return self.__class__(composition, self.unit)

    def __pow__(self, exp):
        if not isinstance(exp, (int, long)):
            raise TypeError("powering index must be an integer.")
        powered_form = powPDF(self.element, exp, self.unit)
        return self.__class__(powered_form, self.unit)

    def __truediv__(self, other):
        invel = other.inverse()
        composition = compositePDF(self.element, invel.element)
        return self.__class__(composition, self.unit)

    def inverse(self):
        """
        Return inverse of the form.
        """
        if self.element == self.unit:
            return copy.deepcopy(self)
        else:
            f_1, f_2, f_3 = self.element
            return self.__class__(reducePDF((f_1, -f_2, f_3)), self.unit)

    def __eq__(self, other):
        if self.__class__ != other.__class__:
            return False
        if self.element == other.element and self.unit == other.unit:
            return True
        else:
            return False

    def __hash__(self):
        val = sum([hash(ele) for ele in self.element])
        val += sum([hash(ele) for ele in self.unit])
        return val

    def __ne__(self, other):
        return not self.__eq__(other)

    def disc(self):
        """
        Return discriminant of the form.
        """
        return self.element[1]**2 - 4*self.element[0]*self.element[2]

    def repOfModule(self):
        ld = self.disc()
        a_m2 = 2*self.element[0]
        rb = -self.element[1]

        return_str = '%s + root(%s) / %s' % (rb, ld, a_m2)
        return return_str

    def __iter__(self):
        """
        Return iterator for form (a, b, c).
        """
        return iter(self.element)

    def __getitem__(self, index):
        """
        Return each component of form (a, b, c),
        """
        return self.element[index]


class ClassGroup(object):
    """
    The class is for class group.
    """
    def __init__(self, disc, classnum, elements=None):
        # element is an element of some class (for example ReducedQuadraticForm
        self.disc = disc
        self.rootoftree = []
        self.rootornot = 0
        if not elements:
            self.elements = []
        else:
            self.elements = copy.deepcopy(elements)
        self.classnum = classnum
        self.expunit = unit_form(disc)

    def __str__(self):
        return_str = "class of ClassGroup:\n"
        return_str = return_str + 'disc is %s\n' % self.disc
        return_str = return_str + 'rootoftree is %s' % self.rootoftree
        return return_str

    def __repr__(self):
        return "ClassGroup(%d, %d, %s)" % (self.disc, self.classnum, self.elements)

    def inserels(self, newlist):
        for newel in newlist:
            self.inserel(newel)

    def inserel(self, newel):
        newestl = copy.deepcopy(newel)
        self.elements.append(newestl)

    def inststree(self, newlist):
        for newel in newlist:
            self.insttree(newel)

    def insttree(self, newel0):
        newel = copy.deepcopy(newel0)
        disc = newel.element[1]**2 - 4*newel.element[0]*newel.element[2]
        if disc != self.disc:
            raise ValueError("this value is not an element of the discriminant")
        if self.rootornot == 0:
            self.rootoftree = [newel, [], []]
            self.rootornot = 1
            return True
        else:
            curntpnt = self.rootoftree
        while curntpnt != []:
            if newel.element == curntpnt[0].element:
                return True
            elif newel.element < curntpnt[0].element:
                curntpnt = curntpnt[1]
            else:
                curntpnt = curntpnt[2]

        curntpnt.append(newel)
        self.elements.append(newel)
        curntpnt.append([])
        curntpnt.append([])

    def search(self, tarel):
        curntpnt = self.rootoftree
        while (curntpnt != []):
            if tarel.element == curntpnt[0].element:
                return curntpnt[0]
            elif tarel.element < curntpnt[0].element:
                curntpnt = curntpnt[1]
            else:
                curntpnt = curntpnt[2]
        return False

    def retel(self):
        copyofroot = copy.deepcopy(self.rootoftree)
        tpa = []
        while copyofroot:
            curntpnt = copyofroot
            while True:
                if curntpnt[1]:
                    curntpnt = curntpnt[1]
                elif curntpnt[2]:
                    curntpnt = curntpnt[2]
                else:
                    tpa.append(curntpnt[0])
                    del curntpnt[:3]
                    break
        return tpa

def class_formula(disc, uprbd):
    """
    Return the approximation of class number 'h' with the given discriminant.
    h = sqrt(|D|)/pi (1 - (D/p)(1/p))^{-1} where p is less than ubound.
    """
    ht = math.sqrt(abs(disc)) / math.pi
    ml = number_unit(disc) / 2

    for factor in prime.generator_eratosthenes(uprbd):
        ml = ml * (1 - (kronecker(disc, factor) / factor))**(-1)
    return int(ht * ml + 0.5)

def class_number(disc, limit_of_disc=100000000):
    """
    Return class number with the given discriminant by counting reduced forms.
    Not only fundamental discriminant.
    """
    if disc & 3 not in (0, 1):
        raise ValueError("a discriminant must be 0 or 1 mod 4")
    if disc >= 0:
        raise ValueError("a discriminant must be negative")
    if -disc >= limit_of_disc:
        warnings.warn("the discriminant seems to have too big absolute value")

    h = 1
    b = disc & 1
    c_b = int(math.sqrt(-disc / 3))

    while b <= c_b:
        q = (b**2 - disc) >> 2
        a = b
        if a <= 1:
            a = 2
        while a**2 <= q:
            if q % a == 0 and gcd.gcd_of_list([a, b, q//a])[0] == 1:
                if a == b or a**2 == q or b == 0:
                    h += 1
                else:
                    h += 2
            a += 1
        b += 2
    return h

def class_group(disc, limit_of_disc=100000000):
    """
    Return the class number and the class group with the given discriminant
    by counting reduced forms. Not only fundamental discriminant.
    """
    if disc & 3 not in (0, 1):
        raise ValueError("a discriminant must be 0 or 1 mod 4")
    if disc >= 0:
        raise ValueError("a discriminant must be negative")
    if -disc >= limit_of_disc:
        warnings.warn("the discriminant seems to have too big absolute value")

    h = 1
    b = disc & 1
    c_b = int(math.sqrt(-disc / 3))

    ret_list = []
    f_a = 1
    if disc & 3 == 0:
        f_b = 0
        f_c = -(disc >> 2)
    else:
        f_b = 1
        f_c = -((disc - 1) >> 2)
    unit = (f_a, f_b, f_c)
    ret_list.append(ReducedQuadraticForm(unit, unit))

    while b <= c_b:
        q = (b**2 - disc) >> 2
        a = b
        if a <= 1:
            a = 2
        while a**2 <= q:
            if q % a == 0 and gcd.gcd_of_list([a, b, q//a])[0] == 1:
                f_c = (b**2 - disc) // (4 * a)
                if a == b or a**2 == q or b == 0:
                    h += 1
                    ret_list.append(ReducedQuadraticForm((a, b, f_c), unit))
                else:
                    h += 2
                    ret_list.append(ReducedQuadraticForm((a, b, f_c), unit))
                    ret_list.append(ReducedQuadraticForm((a, -b, f_c), unit))
            a += 1
        b += 2

    return h, ret_list


class ReducedQuadraticFormForBSGS(ReducedQuadraticForm):
    """
    The class is for reduced quadratic form for *bsgs functions.
    """
    def __init__(self, element, unit):
        ReducedQuadraticForm.__init__(self, element, unit)
        self.ind = -1
        self.alpha = []
        self.beta = []
        self.s_parent = 0
        self.g_parent = 0


def class_number_bsgs(disc):
    """
    Return the class number with the given discriminant.
    """
    if disc & 3 not in (0, 1):
        raise ValueError("a discriminant must be 0 or 1 mod 4")

    if disc >= 0:
        raise ValueError("a discriminant must be negative")

    lx = max(arith1.floorpowerroot(abs(disc), 5), 500 * (math.log(abs(disc)))**2)
    uprbd = int(class_formula(disc, int(lx)) * 3 / 2)
    lwrbd = (uprbd >> 1) + 1
    bounds = [lwrbd, uprbd]

    # get the unit
    element = RetNext(disc)
    ut = element.unit()

    # append the unit to subset of G
    sossp = ClassGroup(disc, 0, [])
    sogsp = ClassGroup(disc, 0, [])
    sossp.insttree(ut)
    sogsp.insttree(ut)

    h = 1 # order
    finished = False
    while not finished:
        mstp1 = bounds[1] - bounds[0]
        if mstp1 <= 1:
            q = 1
        else:
            q = arith1.floorsqrt(mstp1)
            if misc.primePowerTest(mstp1)[1] != 2:
                q += 1
        # get next element
        nt = element.retnext()
        # x is the set of elements of G
        x = [ut, nt ** h]
        if q > 2:
            x.extend([0] * (q - 2))
        # compute small steps
        if x[1] == ut:
            # compute the order of nt
            n = trorder(h, sossp, sogsp, nt, disc)
        else:
            n = trbabysp(q, x, bounds, sossp, sogsp, ut, h, nt, disc)

        # finished?
        finished, h, sossp, sogsp = isfinished_trbsgs(lwrbd, bounds, h, n, sossp, sogsp, nt, disc)

    return h

def class_group_bsgs(disc, classnum, qin):
    """
    Return the construction of the class group with the given discriminant.
    """
    if disc & 3 not in (0, 1):
        raise ValueError("a discriminant must be 0 or 1 mod 4")
    if disc >= 0:
        raise ValueError("a discriminant must be negative")

    matla = []
    lstofelg = []
    lpt = []

    # compute bounds
    qpt = qin[0] ** qin[1]
    uprbd = qpt + 1
    lwrbd = (uprbd >> 1) + 1
    if lwrbd > uprbd:
        raise TypeError("lower bound needs to be less than upper bound")
    if lwrbd <= (uprbd / 2):
        raise TypeError("upper bound / 2 needs to be more than lower bound")
    bounds = [lwrbd, uprbd]

    # get the unit
    uto = unit_form(disc)
    ut = ReducedQuadraticFormForBSGS(uto, uto)

    # append the unit to subset of G
    sossp = ClassGroup(disc, classnum, []) # a subset of G
    sogsp = ClassGroup(disc, classnum, []) # a subset of G
    utwi = copy.deepcopy(ut)
    utwi.alpha.append([0, ut, 1])
    utwi.beta.append([0, ut, 1])
    sossp.insttree(utwi)
    sogsp.insttree(utwi)
    finished = False

    # take a new element of the group.
    indofg = 1
    while not finished:
        # get next element
        while True:
            nt = rand_generator(disc, classnum, qin)
            if nt not in lstofelg:
                lstofelg.append(nt)
                break

        # compute small steps
        assert nt != ut, "random element may not be equal to unit"
        n, tmp_ss, tmp_gs = babyspcv(bounds, sossp, sogsp, utwi, nt, disc, classnum)
        matla.append(setind(n, indofg, tmp_ss, tmp_gs))
        finished, sossp, sogsp = isfinished_bsgscv(n, sossp, sogsp, nt, lpt, qpt, disc, classnum, indofg)
        indofg = indofg + 1
    return lstofelg, matla

##############################################################
# following functions are sub functions for above functions. #
##############################################################

class RetNext:
    """
    Next element producer.
    """
    def __init__(self, disc):
        self.disc = disc
        self.utroot = unit_form(disc)
        self._unit = ReducedQuadraticForm(self.utroot, self.utroot)
        self.cnt = 1
        self.previous = []
        self.elhash = range(int(math.sqrt(abs(disc) // 3)) + 2)

    def unit(self):
        """
        Return reduced quadratic unit form.
        """
        return self._unit

    def retnext(self):
        """
        Return the next random element.
        """
        while True:
            next_form = ReducedQuadraticForm(self._randomele1(), self.utroot)
            if not self.has_in_hash(next_form):
                self.add_to_hash(next_form)
                return next_form

    def _randomele1(self):
        """
        Return a reduced random form with the given discriminant.
        """
        while True:
            nextp = prime.nextPrime(self.cnt)
            self.cnt += 1
            if kronecker(self.disc, nextp) == 1:
                nxtfm = sqroot(self.disc, nextp)
                if not self.previous or nxtfm != self.previous:
                    self.previous = nxtfm
                    return nxtfm

    def add_to_hash(self, form):
        """
        Add a form to hash
        """
        key = form.element[0]
        if isinstance(self.elhash[key], (int, long)):
            self.elhash[key] = [form]
        else:
            self.elhash[key].append(form)

    def has_in_hash(self, form):
        """
        check hash
        """
        key = form.element[0]
        if isinstance(self.elhash[key], (int, long)):
            return False
        return form in self.elhash[key]


def disc(f):
    """
    Return the discriminant of the given quadratic form 'f'.
    f = [a, b, c]
    """
    if len(f) != 3:
        raise TypeError("form must be composed of 3 integers")
    for i in f:
        if not isinstance(i, (int, long)):
            raise TypeError("all components must be integers")
    return (f[1]*f[1] - 4*f[0]*f[2])

def reducePDF(f):
    """
    Return the reduced form of the given positive definite form 'f'.
    f = (f_1, f_2, f_3)
    """
    f_1, f_2, f_3 = f
    if f_1 < 0:
        raise ValueError("f_1 must be positive in quadratic form f=(f_1,f_2,f_3).")
    if (f_2**2 - 4*f_1*f_3) >= 0:
        raise ValueError("discriminant (D= f_2^2 - 4*f_1*f_3) must be negative.")
    if -f_1 < f_2 <= f_1:
        if f_1 > f_3:
            f_2 = -f_2
            f_1, f_3 = f_3, f_1
        else:
            if (f_1 == f_3) and (f_2 < 0):
                f_2 = -f_2
            return (f_1, f_2, f_3)
    while 1:
        q = f_2 // (2*f_1)
        r = f_2 - q*(2*f_1)
        if r > f_1:
            r -= 2*f_1
            q += 1
        f_3 -= ((f_2 + r)>>1)*q
        f_2 = r
        if f_1 > f_3:
            f_2 = -f_2
            f_1, f_3 = f_3, f_1
            continue
        else:
            if (f_1 == f_3) and (f_2 < 0):
                f_2 = -f_2
            return (f_1, f_2, f_3)

def sqrPDF(f):
    """
    Return the square of the given quadratic form 'f'.
    """
    # compute disc and etc
    D = disc(f)
    sogsp = arith1.floorpowerroot(int(abs(D / 4)), 4)
    (u, v, d_1) = euclid_exd(f[1], f[0])

    la = f[0] // d_1
    lb = f[1] // d_1
    lc = (-f[2] * u) % la
    c_1 = la - lc
    if c_1 < lc:
        lc = -c_1

    # partial reduction
    v_2, v_3, z, d, v = parteucl(la, lc, sogsp)

    if z == 0:
        g = (lb * v_3 + f[2]) // d
        a_2 = d ** 2
        c_2 = v_3 ** 2
        b_2 = f[1] + (d + v_3)**2 - a_2 - c_2
        c_2 = c_2 + g * d_1
        return reducePDF((a_2, b_2, c_2))

    e = (f[2] * v + lb * d) // la
    g = (e * v_2 - lb) // v
    b_2 = e * v_2 + v * g
    if d_1 > 1:
        b_2 = d_1 * b_2
        v = d_1 * v
        v_2 = d_1 * v_2

    a_2 = d ** 2
    c_2 = v_3 ** 2
    b_2 = b_2 + (d + v_3) ** 2 - a_2 - c_2
    a_2 = a_2 + e * v
    c_2 = c_2 + g * v_2
    return reducePDF((a_2, b_2, c_2))

def powPDF(f, exp, ut=None):
    """
    Return the powering 'exp' of the given quadratic form 'f'.
    """
    if ut is None:
        D = disc(f)
        ut = unit_form(D)

    if exp == 0:
        return ut
    elif exp == 1:
        return f
    elif f == ut:
        return ut
    if exp < 0:
        lexp = -exp
        sz = (f[0], -f[1], f[2])
    else:
        lexp = exp
        sz = f
    # Right-Left Binary algorithm
    sy = ut
    while True:
        if (lexp & 1) == 1:
            sy = compositePDF(sz, sy)
        lexp >>= 1
        if lexp == 0:
            return sy
        else:
            sz = sqrPDF(sz)

def compositePDF(f_1, f_2):
    """
    Return the reduced form of composition of the given forms 'f_1' and 'f_2'.
    'f_1' and 'f_2' are quadratic forms with same disc.
    """
    if gcd.gcd_of_list(f_1)[0] != 1:
        raise ValueError("coefficients of a quadratic form must be relatively prime")
    if gcd.gcd_of_list(f_2)[0] != 1:
        raise ValueError("coefficients of a quadratic form must be relatively prime")
    if disc(f_1) != disc(f_2):
        raise ValueError("two quadratic forms must have same discriminant")

    if f_1[0] > f_2[0]:
        f_1, f_2 = f_2, f_1

    s = (f_1[1] + f_2[1]) >> 1
    n = f_2[1] - s

    if f_2[0] % f_1[0] == 0:
        y_1 = 0
        d = f_1[0]
    else:
        u, v, d = euclid_exd(f_2[0], f_1[0])
        y_1 = u

    if s % d == 0:
        y_2 = -1
        x_2 = 0
        d_1 = d
    else:
        u, v, d_1 = euclid_exd(s, d)
        x_2 = u
        y_2 = -v

    v_1 = f_1[0] // d_1
    v_2 = f_2[0] // d_1
    r = (y_1*y_2*n - x_2*f_2[2]) % v_1

    b_3 = f_2[1] + 2*v_2*r
    a_3 = v_1*v_2
    c_3 = (f_2[2]*d_1 + r*(f_2[1] + v_2*r)) // v_1

    return reducePDF((a_3, b_3, c_3))

def unit_form(disc):
    """
    Return generated quadratic form with the given discriminant.
    """
    if disc & 3 == 0:
        a = 1
        b = 0
        c = disc // -4
    elif disc & 3 == 1:
        a = 1
        b = 1
        c = (disc - 1) // -4
    else:
        raise ValueError("discriminant is not 0 or 1 mod 4.")
    return (a, b, c)

def kronecker(a, b):
    """
    Compute the Kronecker symbol (a/b) using algo 1.4.10 in Cohen's book.
    """
    tab2 = (0, 1, 0, -1, 0, -1, 0, 1)
    if b == 0:
        if abs(a) != 1:
            return 0
        if abs(a) == 1:
            return 1
    if a & 1 == 0 and b & 1 == 0:
        return 0

    v, b = arith1.vp(b, 2)
    if v & 1 == 0:
        k = 1
    else:
        k = tab2[a & 7]
    if b < 0:
        b = -b
        if a < 0:
            k = -k
    while a:
        v, a = arith1.vp(a, 2)
        if v & 1:
            k *= tab2[b & 7]
        if a & b & 2:
            # both a and b are 3 mod 4
            k = -k
        r = abs(a)
        a = b % r
        b = r
    if b > 1:
        # a and be are not coprime
        return 0
    return k

def number_unit(disc):
    """
    Return the number of units with the given discriminant.
    """
    if disc < -4:
        return 2
    elif disc == -4:
        return 4
    elif disc == -3:
        return 6
    else:
        raise ValueError

def crt(inlist):
    """
    Chinese Remainder Theorem, Algo. 1.3.11 of Cohen's Book.
    """
    k = len(inlist)
    if k < 2:
        raise ValueError("nothing to be done for one element")
    ccj = [None] # gabage to simplify loop index
    ellist = list(inlist)
    ellist.sort()
    modulus = 1
    for j in range(1, k):
        modulus *= ellist[j - 1][1]
        d, tpl = gcd.gcd_of_list([modulus % ellist[j][1], ellist[j][1]])
        if d > 1:
            raise ValueError("moduli are not pairwise coprime")
        ccj.append(tpl[0])

    yj = [ellist[0][0] % ellist[0][1]]
    for j in range(1, k):
        ctp = yj[-1]
        for i in range(j - 2, -1, -1):
            ctp = yj[i] + ellist[i][1] * ctp
        yj.append(((ellist[j][0] - ctp) * ccj[j]) % ellist[j][1])
    outp = yj.pop()
    for j in range(k - 2, -1, -1):
        outp = yj.pop() + (ellist[j][1]) * outp
    return outp

def rand_generator(disc, classnum, qin):
    """
    Return the reduced random quadratic form with given discriminant and order t,
    where t = classnum / a ** b and qin = [a, b].
    """
    q = qin[0]**qin[1]
    unit = unit_form(disc)
    while True:
        elgt1 = randomele(disc, unit)
        elg1 = elgt1 ** (classnum // q)
        if elg1.element == elg1.unit:
            continue
        return elg1

def sqroot(disc, p):
    """
    Return a reduced quadratic form with the given discriminant.
    'disc' is a quadratic residue mod 'p'.
    """
    if p == 2: # if 8 | disc => (disc / 8) = 0, 8 not | disc but 4 | disc => 2
        if (disc & 7) == 0:
            bp = disc
        elif (disc & 3) == 0: # 4 - 4 * odd % 8 => 0
            bp = 2
        elif (disc & 7) == 1: # disc is odd and disc % 8 is 1
            bp = disc
        else: # disc is odd and disc & 3 is 1 => impossible (-5 / 2) = -1
            raise ValueError("disc is odd and disc & 3 is 1 => impossible (-5 / 2) = -1")
    else:
        bpf1 = arith1.modsqrt(disc, p)
        bpf2 = disc
        bp = crt([(bpf1, p), (bpf2, 4)])
    if bp > p:
        bp = 2 * p - bp

    fpt = reducePDF((p, bp, ((bp ** 2) - disc) // (4 * p)))
    return fpt

def randomele(disc, unit):
    """
    Return a reduced random form with the given discriminant and the given unit.
    Also random element is not unit.
    """
    limit = int(math.sqrt(-disc / 3))
    while True:
        a = bigrandom.randrange(1, limit + 1)
        ind = 0
        while ind < 2*a:
            b = bigrandom.randrange(a)
            if bigrandom.randrange(2):
                b = -b
            tp = disc - b**2
            if tp % (-4 * a) == 0:
                c = tp // (-4 * a)
                if gcd.gcd_of_list([a, b, c])[0] != 1:
                    continue
                red = reducePDF((a, b, c))
                if red != unit:
                    return ReducedQuadraticForm(red, unit)
            ind += 1


def isfundamental(disc):
    """
    Determine whether the given discriminant is fundamental or not.
    """
    if disc == 1:
        return False
    spt = misc.squarePart(abs(disc))
    if (disc & 3) == 1:
        return spt == 1
    elif (disc & 3) == 0:
        spt //= 2
        if spt != 1:
            return False
        discof = (disc >> 2) & 3
        return discof == 2 or discof == 3
    return False

def euclid_exd(a, b):
    """
    Return a tuple (u, v, d); they are the greatest common divisor d
    of two integers a and b and u, v such that d = a * u + b * v.
    """
    if not isinstance(a, (int, long)) or not isinstance(b, (int, long)):
        raise TypeError
    u = 1
    d = a
    if b == 0:
        v = 0
        return (u, v, d)
    else:
        v_1 = 0
        v_3 = b

        while 1:
            if v_3 == 0:
                v = (d - a*u) // b
                return (u, v, d)
            q = d // v_3
            t_3 = d % v_3
            t_1 = u - q*v_1
            u = v_1
            d = v_3
            v_1 = t_1
            v_3 = t_3

def parteucl(a, b, sogsp):
    """
    Do extended partial Euclidean algorithm on 'a' and 'b'.
    """
    v = 0
    d = a
    v_2 = 1
    v_3 = b
    z = 0

    while 1:
        if abs(v_3) > sogsp:
            # euclidean step
            q = d // v_3
            t_3 = d % v_3
            t_2 = v - (q * v_2)
            v = v_2
            d = v_3
            v_2 = t_2
            v_3 = t_3
            z = z + 1
        else:
            if z & 1:
                v_2 = -v_2
                v_3 = -v_3
            return (v_2, v_3, z, d, v)

def isfinished_bsgscv(n, sossp, sogsp, nt, lpt, qpt, disc, classnum, indofg):
    """
    Determine whether the bsgs algorithm is finished or not yet.
    This is a submodule called by the bsgs module.
    """
    lpt.append(n)
    sumn = 1
    for nid in lpt:
        sumn = sumn * nid
    if sumn == qpt:
        return True, sossp, sogsp
    elif sumn > qpt:
        raise ValueError

    if n == 1:
        return False, sossp, sogsp
    else:
        tpsq = misc.primePowerTest(n)
        if (tpsq[1] != 0) and ((tpsq[1] & 1) == 0):
            q = arith1.floorsqrt(n)
        else:
            q = arith1.floorsqrt(n) + 1

    ss = sossp.retel()
    new_sossp = ClassGroup(disc, classnum, [])
    tnt = copy.deepcopy(nt)
    for i in range(q):
        base = tnt ** i
        for ssi in ss:
            newel = base * ssi
            if new_sossp.search(newel) is False:
                newel.alpha = ssi.alpha[:]
                lenal = len(newel.alpha)
                sfind = indofg - lenal
                for sit in range(sfind):
                    newel.alpha.append([lenal + sit, 0, 0])
                newel.alpha.append([indofg, tnt, i])
                new_sossp.insttree(newel) # multiple of two elements of G

    y = nt ** q
    ltl = sogsp.retel()
    new_sogsp = ClassGroup(disc, classnum, [])
    for i in range(q + 1):
        base = y ** (-i)
        for eol in ltl:
            newel2 = base * eol
            if new_sogsp.search(newel2) is False:
                newel2.beta = eol.beta[:]
                lenbt = len(newel2.beta)
                gfind = indofg - lenbt
                for git in range(gfind):
                    newel2.beta.append([lenbt + git, 0, 0])
                newel2.beta.append([indofg, tnt, q * (-i)])
                new_sogsp.insttree(newel2) # multiple of two elements of G

    return False, new_sossp, new_sogsp

def ordercv(n, sossp, sogsp, nt, disc, classnum, tmp_ss, tmp_gs):
    """
    n: int
    """
    factor_n = misc.FactoredInteger(n)
    while True:
        c_s1 = ClassGroup(disc, classnum, []) # a subset of G
        lstp_ls = sossp.retel()
        sogsptp = sogsp.retel()
        for p, e in factor_n:
            base = nt ** (int(factor_n) // p)
            for ttp_ls in lstp_ls:
                tmp_c_s1 = base * ttp_ls
                tmp_c_s1.s_parent = ttp_ls
                c_s1.insttree(tmp_c_s1)
            for tmp_ell in sogsptp:
                rete = c_s1.search(tmp_ell)
                if rete != False:
                    factor_n //= p
                    tmp_ss = rete.s_parent
                    tmp_gs = tmp_ell
                    break
            else:
                continue
            break
        else:
            break
    return int(factor_n), tmp_ss, tmp_gs

def giantspcv(q, sz, y, c_s1, bounds, sogsp):
    """
    giant step called from babyspcv.

    q: int
    sz, y: element
    """
    n = bounds[0]
    # sz == x[1] ** n
    # y  == x[1] ** q
    while 1:
        for tpw in sogsp.retel():
            sz1 = sz * tpw
            sz1.g_parent = tpw
            rete = c_s1.search(sz1)
            if rete is not False:
                return n - rete.ind, rete.s_parent, sz1.g_parent
        # continue (sp. 5)
        sz = sz * y
        n = n + q
        if n - q + 1 > bounds[1]:
            raise ValueError("the order is larger than upper bound")

def babyspcv(bounds, sossp, sogsp, utwi, nt, disc, classnum):
    """
    Compute small steps
    """
    mstp1 = bounds[1] - bounds[0]
    if (mstp1 == 0) or (mstp1 == 1):
        q = 1
    else:
        tppm = misc.primePowerTest(mstp1)
        q = arith1.floorsqrt(mstp1)
        if (tppm[1] == 0) or (tppm[1] & 1):
            q += 1

    n_should_be_set = True
    # initialize
    c_s1 = ClassGroup(disc, classnum, []) # a subset of G

    # extracting i = 0 case of main loop
    for ttr in sossp.retel():
        tmpx = ttr
        tmpx.s_parent = ttr # tmpx belongs ttr in the set of smallstep
        # index of the element
        tmpx.ind = 0
        c_s1.insttree(tmpx)

    # main loop
    x_i = nt
    for i in range(1, q):
        for ttr in sossp.retel():
            tmpx = x_i * ttr
            tmpx.s_parent = ttr # tmpx belongs ttr in the set of smallstep
            if n_should_be_set and tmpx == utwi:
                n = i
                tmp_ss = tmpx.s_parent
                tmp_gs = utwi
                n_should_be_set = False
            # index of the element
            tmpx.ind = i
            c_s1.insttree(tmpx)
        x_i = nt * x_i
    assert x_i == nt ** q

    if n_should_be_set:
        sz = nt ** bounds[0]
        n, tmp_ss, tmp_gs = giantspcv(q, sz, x_i, c_s1, bounds, sogsp)
    return ordercv(n, sossp, sogsp, nt, disc, classnum, tmp_ss, tmp_gs)

def trbabysp(q, x, bounds, sossp, sogsp, ut, h, nt, disc):
    """
    Compute small steps.

    q, h: int
    ut: unit element
    nt: element
    """
    c_s1 = ClassGroup(disc, 0, []) # a subset of G
    n_should_be_set = True

    # extracting i = 0 case simplifies the main loop
    for tmpx in sossp.retel():
        tmpx.ind = 0
        c_s1.insttree(tmpx)

    # main loop
    x_i = x[1]
    for i in range(1, q):
        for ttr in sossp.retel():
            tmpx = x_i * ttr
            if n_should_be_set and tmpx == ut:
                n = i
                n_should_be_set = False
            tmpx.ind = i
            c_s1.insttree(tmpx)
            # sort ( if you want to sort it with your estimate,
            # you have to implement '__ge__' method of the class with your way.)
        x_i = x[1] * x_i

    if n_should_be_set:
        n = trgiantsp(q, x[1] ** bounds[0], x_i, c_s1, bounds, sogsp)
    return trorder(n * h, sossp, sogsp, nt, disc)

def trgiantsp(stride_index, pivot, stride, c_s1, bounds, sogsp):
    """
    Compute giant steps.

    stride_index: int
    pivot, stride: element
    """
    pivot_index = bounds[0]
    # pivot == x[1] ** pivot_index
    # stride = x[1] ** stride_index
    while True:
        for tpw in sogsp.retel():
            rete = c_s1.search(pivot * tpw)
            if rete is not False:
                return pivot_index - rete.ind
        pivot, pivot_index = pivot * stride, pivot_index + stride_index
        if pivot_index - stride_index + 1 > bounds[1]:
            raise ValueError("the order is larger than upper bound")

def trorder(n, sossp, sogsp, nt, disc):
    """
    Compute the order.

    n: int
    nt: element
    """
    factor_n = misc.FactoredInteger(n)
    while True:
        c_s1 = ClassGroup(disc, 0, [])
        lstp_ls = sossp.retel()
        sogsptp = sogsp.retel()
        for p, e in factor_n:
            # initialize c_s1
            base = nt ** (int(factor_n) // p)
            for ttp_ls in lstp_ls:
                c_s1.insttree(base * ttp_ls)
            # search in c_s1
            for tmp_ell in sogsptp:
                rete = c_s1.search(tmp_ell)
                if rete != False:
                    factor_n //= p
                    break
            else:
                continue
            break
        else:
            return int(factor_n)

def isfinished_trbsgs(lwrbd, bounds, h, n, sossp, sogsp, nt, disc):
    """
    Determine whether bsgs is finished or not yet.
    This is a submodule called by the bsgs module.

    lwrbd, h, n: int
    nt: element
    """
    h *= n
    if h >= lwrbd:
        result = True
    elif n == 1:
        result = False
    else:
        bounds[0] = (bounds[0] + n - 1) // n # ceil of lower bound // n
        bounds[1] = bounds[1] // n # floor of upper bound // n
        q = arith1.floorsqrt(n)
        if misc.primePowerTest(n)[1] != 2:
            q = arith1.floorsqrt(n) + 1 # ceil of sqrt
        sossp, sogsp = _update_subgrps(q, nt, sossp, sogsp, disc)
        result = False

    return result, h, sossp, sogsp

def _update_subgrps(q, element, sossp, sogsp, discriminant):
    """
    update sossp and sogsp
    """
    new_sossp = ClassGroup(discriminant, 0, [])
    new_sogsp = ClassGroup(discriminant, 0, [])

    unit = element ** 0
    ithpow = unit
    for i in range(q):
        for ssi in sossp.retel():
            newel = ithpow * ssi
            if new_sossp.search(newel) is False:
                new_sossp.insttree(newel)
        ithpow = ithpow * element
    assert ithpow == element ** q

    qithpow = unit
    for i in range(q + 1):
        for eol in sogsp.retel():
            newel = qithpow * eol
            if new_sogsp.search(newel) is False:
                new_sogsp.insttree(newel)
        qithpow *= ithpow

    return new_sossp, new_sogsp

def setind(n, indofg, tmp_ss, tmp_gs):
    """
    """
    lgtinlst = indofg
    if lgtinlst == 1:
        return [n]
    tmp_mt = [n]
    for idofel in range(lgtinlst):
        if idofel == 0:
            continue
        try:
            if type(tmp_ss.alpha[idofel][1]) != int:
                ioind = tmp_ss.alpha[idofel][2]
            else:
                ioind = 0
        except IndexError:
            ioind = 0
        except:
            raise ValueError
        try:
            if type(tmp_gs.beta[idofel][1]) != int:
                joind = tmp_gs.beta[idofel][2]
            else:
                joind = 0
        except IndexError:
            joind = 0
        except:
            raise ValueError
        tmp_mt.append(ioind - joind)
    return tmp_mt
