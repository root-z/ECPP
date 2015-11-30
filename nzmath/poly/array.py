import nzmath.arith1 as arith1

def check_zero_poly(coefficients):
    """
    This function checks whether all elements of coefficients equal
    zero or not.  If all elements of coefficients equal zero, this
    function returns True.  Else this function returns False.
    """
    for i in coefficients:
        if i != 0:
            return False
    return True

def arrange_coefficients(coefficients):
    """
    This function arranges coefficient.
    For example, [1,2,0,3,0] => [1,2,0,3].
    """
    if check_zero_poly(coefficients):
        return [0]
    while coefficients[-1] == 0:
        coefficients = coefficients[:-1]
    return coefficients

class ArrayPoly():
    """
    Polynomial with integer number coefficients.
    Coefficients has to be a initializer for list.
    """

    def __init__(self, coefficients = [0]):
        """
        Initialize the polynomial.
        coefficients:initializer for polynomial coefficients
        """
        self.coefficients = arrange_coefficients(coefficients)
        self.degree = len(self.coefficients) - 1

    def coefficients_to_dict(self):
        """
        Return coefficients as dict.
        """
        if self.coefficients == [0]:
            return {0:0}
        dict_coefficients = {}
        for i in range(self.degree + 1):
            if self.coefficients[i] != 0:
                dict_coefficients.update({i:self.coefficients[i]})
        return dict_coefficients

    def __repr__(self):
        poly_repr = self.coefficients_to_dict()
        return "ArrayPoly : %s" % poly_repr

    def __str__(self):
        poly_str = self.coefficients_to_dict()
        return "polynomial : %s" % poly_str

    def __add__(self, other):
        """
        self + other
        """
        add = []
        deg1 = self.degree
        deg2 = other.degree
        if deg1 >= deg2:
            long_coefficients = self.coefficients[:]
            short_coefficients = other.coefficients[:]
            deg = deg2 + 1
        else:
            long_coefficients = other.coefficients[:]
            short_coefficients = self.coefficients[:]
            deg = deg1 + 1
        for i in range(deg):
            add.append(long_coefficients[i] + short_coefficients[i])
        add += long_coefficients[deg:]
        return ArrayPoly(add)

    def __sub__(self, other):
        """
        self - other
        """
        sub = []
        deg1 = self.degree
        deg2 = other.degree
        if deg1 >= deg2:
            long_coefficients = self.coefficients[:]
            short_coefficients = other.coefficients[:]
            deg = deg2 + 1
        else:
            long_coefficients = other.coefficients[:]
            short_coefficients = self.coefficients[:]
            deg = deg1 + 1
        for i in range(deg):
            sub.append(long_coefficients[i] - short_coefficients[i])
        sub += long_coefficients[deg:]
        return ArrayPoly(sub)

    def scalar_mul(self, scalar):
        """
        Return the result of scalar multiplicaton.
        """
        scalar_mul = [i * scalar for i in self.coefficients]
        return ArrayPoly(scalar_mul)

    def upshift_degree(self, slide):
        """
        Return the polynomial obtained by shifting upward all terms
        with degrees of 'slide'.
        """
        if slide == 0:
            return ArrayPoly(self.coefficients[:])
        up_degree = [0] * slide
        return ArrayPoly(up_degree + self.coefficients[:])

    def downshift_degree(self, slide):
        """
        Return the polynomial obtained by shifting downward all terms
        with degrees of 'slide'.
        """
        if slide == 0:
            return ArrayPoly(self.coefficients[:])
        elif slide > self.degree:
            return ArrayPoly()
        down_degree= self.coefficients[slide:]
        return ArrayPoly(down_degree)

    def __eq__(self, other):
        """
        self == other
        """
        self_dict_coefficients = self.coefficients_to_dict()
        other_dict_coefficients = other.coefficients_to_dict()
        return self_dict_coefficients == other_dict_coefficients

    def __ne__(self, other):
        """
        self != other
        """
        return not self.__eq__(other)

    def __mul__(self, other):
        """
        self * other
        """
        total = [0] * (self.degree + other.degree + 1)
        for i in range(self.degree + 1):
            for j in range(other.degree + 1):
                total[i + j] += self.coefficients[i] * other.coefficients[j]
        return ArrayPoly(total)

    def power(self):
        """
        self * self
        """
        total = [0] * (self.degree + self.degree + 1)
        for i in range(self.degree + 1):
            total[i + i] += self.coefficients[i] * self.coefficients[i]
            for j in range(i + 1, self.degree + 1):
                total[i + j] += ((self.coefficients[i] * self.coefficients[j]) << 1)
        return ArrayPoly(total)

    def split_at(self, split_point):
        """
        Return tuple of two polynomials, which are splitted at the
        given degree.  The term of the given degree, if exists,
        belongs to the lower degree polynomial.
        """
        split = self.coefficients[:]
        if split_point == 0:
            return (ArrayPoly(), ArrayPoly(split))
        elif split_point >= self.degree:
            return (ArrayPoly(split), ArrayPoly())
        split1 = split[:split_point + 1]
        split2 = split[:]
        for i in range(split_point + 1):
            split2[i] = 0
        return (ArrayPoly(split1), ArrayPoly(split2))

    def FFT_mul(self, other):
        """
        self * other by Fast Fourier Transform.
        """
        coefficients1 = [abs(i) for i in self.coefficients]
        coefficients2 = [abs(i) for i in other.coefficients]
        bound1 = max(coefficients1)
        bound2 = max(coefficients2)
        bound = bound1 * bound2 * (max(self.degree, other.degree) + 1)
        bound = ceillog(bound, 2)
        bound = 1 << bound
        if bound < (self.degree + other.degree + 1):
            bound = self.degree + other.degree + 1
            bound = ceillog(bound, 2)
            bound = 1 << bound
        fft1 = ArrayPoly(self.coefficients[:])
        fft2 = ArrayPoly(other.coefficients[:])
        fft_mul1 = FFT(fft1, bound)
        fft_mul2 = FFT(fft2, bound)
        fft_mul = []
        for i in range(bound):
            fft_mul.append(fft_mul1[i] * fft_mul2[i])
        #print fft_mul
        total = reverse_FFT(fft_mul, bound)
        #print total
        return ArrayPoly(total)


class ArrayPolyMod(ArrayPoly):
    """
    Polynomial with modulo n coefficients, n is a nutural number.
    Coefficients has to be a initializer for list.
    """

    def __init__(self, coefficients, mod):
        """
        Initialize the polynomial.
        coefficients:initializer for polynomial coefficients
        mod:initializer for coefficients modulo mod
        """
        if mod <= 0:
            raise ValueError("Please input a positive integer in 'mod'.")
        mod_coefficients = [i % mod for i in coefficients]
        ArrayPoly.__init__(self, mod_coefficients)
        self.mod = mod

    def __repr__(self):
        poly_repr = self.coefficients_to_dict()
        return "polynomial_mod(%d):%s" % (self.mod, poly_repr)

    def __str__(self):
        poly_str = self.coefficients_to_dict()
        return "polynomial_mod(%d):%s" % (self.mod, poly_str)

    def __add__(self, other):
        """
        self + other
        """
        if self.mod != other.mod:
            raise ValueError("mod mismatch")
        add = ArrayPoly.__add__(self, other)
        return ArrayPolyMod(add.coefficients, self.mod)

    def __sub__(self, other):
        """
        self - other
        """
        if self.mod != other.mod:
            raise ValueError("mod mismatch")
        sub = ArrayPoly.__sub__(self, other)
        return ArrayPolyMod(sub.coefficients, self.mod)

    def scalar_mul(self, scalar):
        """
        Return the result of scalar multiplicaton.
        """
        scalar_mul = ArrayPoly.scalar_mul(self, scalar)
        return ArrayPolyMod(scalar_mul.coefficients, self.mod)

    def upshift_degree(self, slide):
        """
        Return the polynomial obtained by shifting upward all terms
        with degrees of 'slide'.
        """
        up_degree = ArrayPoly.upshift_degree(self, slide)
        return ArrayPolyMod(up_degree.coefficients, self.mod)

    def downshift_degree(self, slide):
        """
        Return the polynomial obtained by shifting downward all terms
        with degrees of 'slide'.
        """
        down_degree = ArrayPoly.downshift_degree(self, slide)
        return ArrayPolyMod(down_degree.coefficients, self.mod)

    def __eq__(self, other):
        """
        self == other
        """
        if self.mod != other.mod:
            raise ValueError("mod mismatch")
        eq = ArrayPoly.__eq__(self, other)
        return eq

    def __ne__(self, other):
        """
        self != other
        """
        if self.mod != other.mod:
            raise ValueError("mod mismatch")
        ne = ArrayPoly.__ne__(self, other)
        return ne

    def __mul__(self, other):
        """
        self * other
        """
        if self.mod != other.mod:
            raise ValueError("mod mismatch")
        total = [0] * (self.degree + other.degree + 1)
        for i in range(self.degree + 1):
            for j in range(other.degree + 1):
                total[i + j] = (total[i + j] + self.coefficients[i] * other.coefficients[j]) % self.mod
        return ArrayPolyMod(total, self.mod)

    def power(self):
        """
        self * self
        """
        total = [0] * (self.degree + self.degree + 1)
        for i in range(self.degree + 1):
            total[i + i] = (total[i + i] + self.coefficients[i] * self.coefficients[i]) % self.mod
            for j in range(i + 1, self.degree + 1):
                total[i + j] = (total[i + j] + ((self.coefficients[i] * self.coefficients[j]) << 1)) % self.mod
        return ArrayPolyMod(total, self.mod)

    def split_at(self, split_point):
        """
        Return tuple of two polynomials, which are splitted at the
        given degree.  The term of the given degree, if exists,
        belongs to the lower degree polynomial.
        """
        if split_point == 0:
            return (ArrayPolyMod([0], self.mod), ArrayPolyMod(self.coefficients, self.mod))
        elif split_point >= self.degree:
            return (ArrayPolyMod(self.coefficients, self.mod), ArrayPolyMod([0], self.mod))
        split = ArrayPoly.split_at(self, split_point)
        split1 = split[0].coefficients
        split2 = split[1].coefficients
        return (ArrayPolyMod(split1, self.mod), ArrayPolyMod(split2, self.mod))

    def FFT_mul(self, other):
        """
        self * other by Fast Fourier Transform.
        """
        if self.mod != other.mod:
            raise ValueError("mod mismatch")
        bound1 = max(self.coefficients)
        bound2 = max(other.coefficients)
        bound = bound1 * bound2 * (max(self.degree, other.degree) + 1)
        bound = ceillog(bound, 2)
        bound = 1 << bound
        if bound < (self.degree + other.degree + 1):
            bound = self.degree + other.degree + 1
            bound = ceillog(bound, 2)
            bound = 1 << bound
        fft1 = ArrayPolyMod(self.coefficients[:], self.mod)
        fft2 = ArrayPolyMod(other.coefficients[:], other.mod)
        fft_mul1 = FFT(fft1, bound)
        fft_mul2 = FFT(fft2, bound)
        fft_mul = []
        for i in range(bound):
            fft_mul.append(fft_mul1[i] * fft_mul2[i])
        total = reverse_FFT(fft_mul, bound)
        return ArrayPolyMod(total, self.mod)

#
#Some functions for FFT(Fast Fourier Transform).
#
def min_abs_mod(a, b):
    """
    This function returns absolute minimum modulo of a over b.
    """
    if a >= 0:
        return a - b * ((a + a + b) // (b + b))
    a = -a
    return -(a - b * ((a + a + b) // (b + b)))

def bit_reverse(n, bound):
    """
    This function returns the result reversed bit of n.
    bound:number of significant figures of bit.
    """
    bit = []
    while n > 0:
        bit.append(n & 1)
        n = n >> 1
    bound_bitlen = ceillog(bound, 2)
    if bound_bitlen > len(bit):
        bit.extend([0] * (bound_bitlen - len(bit)))
    bit.reverse()
    total = 0
    count = 0
    for i in bit:
        if i != 0:
            total += 1 << count
        count += 1
    return total

def ceillog(n, base=2):
    """
    Return ceiling of log(n, 2)
    """
    floor = arith1.log(n, base)
    if base ** floor == n:
        return floor
    else:
        return floor + 1

def perfect_shuffle(List):
    """
    This function returns list arranged by divide-and-conquer method.
    """
    length = len(List)
    shuffle = [0] * length
    for i in range(length):
        rebit = bit_reverse(i, length)
        shuffle[rebit] = List[i]
    return shuffle

def FFT(f, bound):
    """
    Fast Fourier Transform.
    This function returns the result of valuations of f by fast fourier transform
    against number of bound different values.
    """
    count = 1 << (bound >> 1)
    mod = count + 1
    if isinstance(f, ArrayPoly):
        coefficients = f.coefficients[:]
    else:
        coefficients = f[:]
    coefficients.extend([0] * (bound - len(coefficients)))
    shuf_coefficients = perfect_shuffle(coefficients)
    shuf_coefficients = [min_abs_mod(i, mod) for i in shuf_coefficients]
    bound_log = arith1.log(bound, 2)
    for i in range(1, bound_log + 1):
        m = 1 << i
        wm = count
        wm = wm % mod
        w = 1
        m1 = m >> 1
        for j in range(m1):
            for k in range(j, bound, m):
                m2 = k + m1
                plus = shuf_coefficients[k]
                minus = w * shuf_coefficients[m2]
                shuf_coefficients[k] = (plus + minus) % mod
                shuf_coefficients[m2] = (plus - minus) % mod
            w = w * wm
        shuf_coefficietns = [min_abs_mod(i, mod) for i in shuf_coefficients]
        count = arith1.floorsqrt(count)
    return shuf_coefficients

def reverse_FFT(values, bound):
    """
    Reverse Fast Fourier Tfransform.
    """
    mod = (1 << (bound >> 1)) + 1
    shuf_values = values[:]
    reverse_FFT_coeffficients = FFT(shuf_values, bound)
    inverse = arith1.inverse(bound, mod)
    reverse_FFT_coeffficients = [min_abs_mod(inverse * i, mod) for i in reverse_FFT_coeffficients]
    reverse_FFT_coeffficients1 = reverse_FFT_coeffficients[:1]
    reverse_FFT_coeffficients2 = reverse_FFT_coeffficients[1:]
    reverse_FFT_coeffficients2.reverse()
    reverse_FFT_coeffficients_total = reverse_FFT_coeffficients1 + reverse_FFT_coeffficients2
    return reverse_FFT_coeffficients_total
