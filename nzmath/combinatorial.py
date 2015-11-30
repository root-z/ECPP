"""
Combinatorial functions
"""

import bisect
import itertools
from nzmath.rational import Integer, Rational


def binomial(n, m):
    """
    The binomial coefficient.
    binomial(n, m) returns n ! / ((n - m) ! * m !).

    n must be a positive integer and m must be a non-negative integer.
    For convinience, binomial(n, n+i) = 0 for positive i, and
    binomial(0,0) = 1.

    In other cases, it raises an exception.
    """
    if not isinstance(n, (int, long)):
        raise TypeError("integer is expected, %s detected." % n.__class__)
    if not isinstance(m, (int, long)):
        raise TypeError("integer is expected, %s detected." % m.__class__)
    if n == m >= 0 or m == 0 and n > 0:
        return 1
    if n <= 0:
        raise ValueError("non-positive number: %d" % n)
    if m < 0:
        raise ValueError("negative number: %d" % m)
    if n < m:
        return 0
    if m*2 > n:
        m = n - m
    retval = n
    for i in range(1, m):
        retval *= (n - i)
        retval /= (i + 1)
    return Integer(retval)

def factorial(n):
    """
    Return n! for non negative integer n.
    """
    if not isinstance(n, (int, long)):
        raise TypeError("integer is expected, %s detected." % n.__class__)
    elif n < 0:
        raise ValueError("argument must not be a negative integer.")
    elif n < 1500:
        # naive loop is best for small n.
        result = 1
        for i in xrange(2, n + 1):
            result *= i
        return Integer(result)
    # The following algorithm keeps temporary results rather small as
    # possible in order to make the function run faster than the naive
    # loop.
    l = range(1, n + 1)
    while len(l) > 1:
        for i in xrange(len(l) >> 1):
            l[i] *= l.pop()
    return Integer(l.pop())

def bernoulli(n):
    """
    Return n-th Bernoulli number.
    """
    if n != 1 and n & 1:
        return Integer(0)
    B = {0:Integer(1),
         1:Rational(-1, 2)}
    for i in range(2, n + 1, 2):
        a = B[0] + (i + 1) * B[1]
        for j in range(2, i, 2):
            a += binomial(i + 1, j) * B[j]
        B[i] = -a / (i + 1)
    return B[n]

def catalan(n):
    """
    Return n-th Catalan number.
    """
    return binomial(2*n, n) // (n+1)

def euler(n):
    """
    Return n-th Euler number.
    """
    if n & 1:
        return Integer(0)
    E = {0:Integer(1)}
    for i in range(2, n + 1, 2):
        E[i] = sum([-binomial(i, j) * E[j] for j in range(0, i, 2)])
    return E[n]

def fallingfactorial(n, m):
    """
    Return the falling factorial; n to the m falling, i.e. n(n-1)..(n-m+1).

    For Example:
    >>> fallingfactorial(7, 3)
    210
    """
    r = 1
    for i in range(n, n-m, -1):
        r *= i
    return r

def risingfactorial(n, m):
    """
    Return the rising factorial; n to the m rising, i.e. n(n+1)..(n+m-1).

    For example:
    >>> risingfactorial(7, 3)
    504
    """
    r = 1
    for i in range(n, n+m):
        r *= i
    return r

def multinomial(n, parts):
    """
    Return multinomial coefficient.

    parts MUST be a sequence of natural numbers and n==sum(parts) holds.
    """
    if n != sum(parts):
        raise ValueError("sum of parts must be equal to n.")
    for part in parts:
        if not isinstance(part, (int, long)) or part < 0:
            raise ValueError("parts must be a sequence of natural numbers.")
    # binomial case
    if len(parts) == 2:
        return binomial(n, parts[0])
    # other cases
    result = factorial(n)
    for part in parts:
        if part >= 2: # 0! = 1! = 1 are negligible
            result //= factorial(part)
    return result

def stirling1(n, m):
    """
    Stirling number of the first kind.

    Let s denote the Stirling number, (x)_n falling factorial, then
      (x)_n = \sum_{i=0}^{n} s(n, i) * x**i
    and s satisfies the recurrence relation:
      s(n, m) = s(n-1, m-1) - (n-1)*s(n-1, m)
    """
    if n == m:
        return 1
    elif n < m or n*m < 0:
        return 0
    elif m == 0:
        return 0
    elif m == 1:
        if n & 1:
            sign = 1
        else:
            sign = -1
        return sign * factorial(n - 1)
    elif m == n - 1:
        return -binomial(n, 2)
    else:
        # recursively calculate by s(n, m) = s(n-1, m-1) - (n-1)*S(n-1, m)
        slist = (1,) # S(1, 1) = 1
        for i in range(1, n):
            l, u = max(1, m - n + i + 1), min(i + 2, m + 1)
            if l == 1 and len(slist) < n - l:
                nlist = [-i * slist[0]]
            else:
                nlist = []
            for j in range(len(slist) - 1):
                nlist.append(slist[j] - i * slist[j + 1])
            if len(slist) <= u - l:
                nlist.append(slist[-1])
            slist = tuple(nlist)
        return slist[0]

def stirling2(n, m):
    """
    Stirling number of the second kind.

    Let S denote the Stirling number, (x)_i falling factorial, then:
      x**n = \sum_{i=0}^{n} S(n, i) * (x)_i
    S satisfies:
      S(n, m) = S(n-1, m-1) + m*S(n-1, m)
    """
    if n == m:
        return 1
    elif n < m or n * m <= 0:
        return 0
    elif n < 0:
        return stirling2_negative(-m, -n)
    elif m == 1:
        return 1
    elif m == 2:
        return 2**(n - 1) - 1
    elif m == n - 1:
        return n * m >> 1
    else:
        r = 0
        b = m
        l = 1
        while 1:
            if l & 1:
                r -= b * l**n
            else:
                r += b * l**n
            if l == m:
                break
            b = (b * (m - l)) // (l + 1)
            l += 1
        if m & 1:
            r = -r
        return r // factorial(m)

def stirling2_negative(n, m):
    """
    Stiring number of the second kind extended to negative numbers.

    Let S# denote the extended Stirling number, S the original, then
    S#(n, m) = S(-m, -n) and extended by the recurrence relation:
      S#(n, m) = S#(n-1, m-1) + (n-1)*S#(n-1, m)
    """
    if n == m:
        return 1
    elif n < m or n * m <= 0:
        return 0
    elif n < 0:
        return stirling2(-m, -n)
    elif m == 1:
        return factorial(n - 1)
    else: # return S(n-1,m-1)+(n-1)*S(n-1,m)
        slist = [1]
        i = 1
        while i < n:
            l, u = max(1, m - n + i + 1), min(i + 2, m + 1)
            if len(slist) < u - l:
                nlist = [i * slist[0]]
            else:
                nlist = []
            for j in range(len(slist) - 1):
                nlist.append(slist[j] + i * slist[j + 1])
            if len(slist) <= u - l:
                nlist.append(slist[-1])
            slist = nlist
            i += 1
        return slist[-1]

def bell(n):
    """
    Bell number.

    The Bell number b is defined by:
      b(n) = \sum_{i=0}^{n} S(n, i)
    where S denotes Stirling number of the second kind.
    """
    return sum([stirling2(n, i) for i in range(n + 1)])


# generators
def combination_index_generator(n, m):
    """
    Generate indices of m elment subsets of n element set.

    The number of generated indices is binomial(n, m).

    For example,
    combinationIndexGenerator(5,3) generates the following lists:
        [0, 1, 2]
        [0, 1, 3]
        [0, 1, 4]
        [0, 2, 3]
        [0, 2, 4]
        [0, 3, 4]
        [1, 2, 3]
        [1, 2, 4]
        [1, 3, 4]
        [2, 3, 4]
    """
    assert n >= m > 0
    idx = range(m)
    while True:
        yield list(idx)
        for i in range(1, m+1):
            if idx[-i] != n-i:
                idx[-i] += 1
                for j in range(i-1, 0, -1):
                    idx[-j] = idx[-j-1] + 1
                break
        else:
            raise StopIteration


def permutation_generator(n):
    """
    Generate all permutations of n elements as lists.

    The number of generated lists is n!, so be careful to use big n.

    For example,
    permutationGenerator(3) generates the following lists:
        [0, 1, 2]
        [0, 2, 1]
        [1, 0, 2]
        [1, 2, 0]
        [2, 0, 1]
        [2, 1, 0]
    """
    # traverse tree by depth first
    perm = range(n)
    unused = []
    while True:
        # leaf is reached, thus yield the value.
        yield list(perm)

        # track back until node with subtree yet to be traversed
        last = n - 1
        unused.append(perm[-1])
        while last and perm[last - 1] > unused[-1]:
            last -= 1
            unused.append(perm[last])

        # exhausted
        if not last:
            break

        # assert unused == sorted(unused)
        # replace with just bigger than perm[last - 1]
        index = bisect.bisect(unused, perm[last - 1])
        unused[index], perm[last - 1] = perm[last - 1], unused[index]
        # replace remaining part
        perm[last:] = unused
        del unused[:]


def dyck_word_generator(n, alphabet=(0, 1)):
    """
    Generate all Dyck words of length 2*n as tuples.

    The Dyck words are words on a two character alphabet.
    The number of each character in a word is equal, 
    and the number of the second character never exceeds the first
    in any initial parts of the word.

    The number of generated words is the n-th Catalan number.

    The alphabet is {0, 1} by default, but you can pass it into the
    optional argument 'alphabet'.

    For example,
    >>> for word in dyck_word_generator(3, alphabet=("(", ")")):
    ...     print "".join(word)
    ... 
    ()()()
    ()(())
    (())()
    (()())
    ((()))
    >>> 
    """
    if n == 0:
        yield ()
        raise StopIteration()
    else:
        # memo[i] is a list of Dyck words of length 2*i
        memo = [[()]]
        alpha, beta = (alphabet[0],), (alphabet[1],)

    # prepare up to n-1
    for i in range(1, n):
        memo.append([])
        for j in range(i):
            for p in memo[j]:
                prefix = alpha + p + beta
                for q in memo[i - j - 1]:
                    memo[i].append(prefix + q)
        
    # output
    for j in range(n):
        for p in memo[j]:
            prefix = alpha + p + beta
            for q in memo[n - j - 1]:
                yield prefix + q


# partition related functions and classes
def partition_generator(n, maxi=None):
    """
    Generate partitions of n.
    If maxi is given, then parts are limited to at most maxi.
    """
    if not n:
        yield ()
        raise StopIteration
    if maxi is None or maxi > n:
        maxi = n
    partition = [maxi]
    rest = n - maxi
    while True:
        key = partition[-1]
        while rest >= key:
            partition.append(key)
            rest -= key
        if rest:
            partition.append(rest)

        yield tuple(partition)

        try:
            # wind up all 1's.
            first_one = partition.index(1)
            if not first_one:
                # partition was [1]*n means all partitions have been generated.
                raise StopIteration
            rest = len(partition) - first_one
            del partition[first_one:]
        except ValueError:
            # 1 is not found
            rest = 0

        partition[-1] -= 1
        rest += 1


class PartitionDriver(object):
    """
    PartitionDriver knows how to construct partitions.

    This class can generate all partitions without any additional
    condition, yet it should be usable as a base class for subclasses
    to generate conditioned partitions, such as limited maximum, with
    odd parts, with distinct parts, etc.
    """
    def __init__(self):
        """
        Initialize the class.

        Public attributes are partition and rest, though they are not
        expected to be modified directly.
        """
        self.partition = []
        self.rest = 0

    def set_parameters(self, integer):
        """
        Set parameters for generate:
          integer: the number to be partitioned

        The method signature may differ class to class.
        """
        self.partition = []
        self.rest = integer

    def ispreconditioned(self):
        """
        Check whether preconditions are satisfied or not.

        For example, partition into evens needs rest be even.
        """
        return True

    def pull(self):
        """
        Pull the maximum part obtainable from rest to partition.

        When you override the method, please try to make it precise.
        If it is unavoidable to do some failing attempt, raising a
        LookupError is the signal to tell that no parts can be pulled.
        """
        if self.partition:
            if self.rest >= self.partition[-1]:
                part = self.partition[-1]
            else:
                part = self.rest
        else:
            part = self.rest

        self.partition.append(part)
        self.rest -= part

    def backtrack(self):
        """
        Back track the history from either broken or complete partition.
        """
        try:
            if self.partition[-1] == 1:
                # wind up all 1's.
                first_one = self.partition.index(1)
                self.rest = len(self.partition) - first_one
                del self.partition[first_one:]
            else:
                self.rest = 0
        except Exception:
            print self.partition
            print self.rest

    def push(self):
        """
        On each end of backtrack, the tail of partition has to be
        decreased and the decrement is cancelled by an increment of
        rest.
        """
        self.partition[-1] -= 1
        self.rest += 1

    def shouldterminate(self):
        """
        Knowing that the generation should terminate, return True to
        stop the loop.
        """
        return not self.partition

    def generate(self, *args):
        """
        Generate all partitions of integer

        Variable arguments *args is used for inheritance.
        """
        self.set_parameters(*args)
        if not self.ispreconditioned():
            raise StopIteration

        while True:
            while self.rest:
                try:
                    self.pull()
                except LookupError:
                    break

            if not self.rest:
                yield tuple(self.partition)

            self.backtrack()
            if self.shouldterminate():
                break
            self.push()


class LimitedMaximumPartitionDriver(PartitionDriver):
    """
    Only limit the muximum of parts.
    """
    def __init__(self):
        PartitionDriver.__init__(self)
        self.limit = 0

    def set_parameters(self, integer, limit):
        """
        Set parameters for generate:
          integer: the number to be partitioned
          limit: the maximum of parts
        """
        self.partition = []
        self.rest = integer
        self.limit = limit

    def pull(self):
        """
        Pull the maximum part obtainable from rest to partition.
        """
        if self.partition:
            if self.rest >= self.partition[-1]:
                part = self.partition[-1]
            else:
                part = self.rest
        elif self.rest >= self.limit:
            part = self.limit
        else:
            part = self.rest

        self.partition.append(part)
        self.rest -= part


class OddPartitionDriver(PartitionDriver):
    """
    All parts are odd.
    """
    def __init__(self):
        PartitionDriver.__init__(self)

    def pull(self):
        """
        Pull the maximum odd part obtainable from rest to partition.
        """
        if self.partition and self.rest >= self.partition[-1]:
            part = self.partition[-1]
        elif self.rest & 1:
            part = self.rest
        else:
            part = self.rest - 1

        self.partition.append(part)
        self.rest -= part

    def push(self):
        """
        On each end of backtrack, the tail of partition has to be
        decreased and the decrement is cancelled by an increment of
        rest keeping oddity of parts.
        """
        self.partition[-1] -= 2
        self.rest += 2


class OddMaximumPartitionDriver(LimitedMaximumPartitionDriver, OddPartitionDriver):
    def __init__(self):
        # Calling LimitedMaximum suffices, since extra work isn't needed for Odd
        LimitedMaximumPartitionDriver.__init__(self)

    def set_parameters(self, integer, limit):
        """
        Set parameters for generate:
          integer: the number to be partitioned
          limit: the maximum of parts
        """
        if not (limit & 1):
            limit -= 1
        LimitedMaximumPartitionDriver.set_parameters(self, integer, limit)

    def pull(self):
        """
        Pull the maximum odd part obtainable from rest to partition.
        """
        if self.partition and self.rest >= self.partition[-1]:
            part = self.partition[-1]
        elif self.rest >= self.limit:
            part = self.limit
        elif self.rest & 1:
            part = self.rest
        else:
            part = self.rest - 1

        self.partition.append(part)
        self.rest -= part


def partition_into_odd_generator(n, maxi=None):
    if maxi is None:
        driver = OddPartitionDriver()
        return driver.generate(n)
    else:
        driver = OddMaximumPartitionDriver()
        return driver.generate(n, maxi)


def partition_numbers_upto(n):
    """
    Return the partition numbers for 0 to '''n''' (inclusive).
    """
    penta = list(itertools.izip(
        itertools.takewhile(lambda k: k <= n, _pentagonal()),
        itertools.cycle((1, 1, -1, -1))))
    p = [1]
    for i in range(1, n + 1):
        s = 0
        for k, sign in penta:
            if k > i:
                break
            s += sign * p[i - k]
        p.append(s)
    return p

def _pentagonal():
    """
    Generates pentagonal and skew pentagonal numbers.
    (1, 2, 5, 7, 12, 15, ...)
    """
    j = 1
    while True:
        yield j*(3*j - 1) >> 1
        yield j*(3*j + 1) >> 1
        j += 1

def partition_number(n):
    """
    Return the partition number for '''n'''.
    """
    return partition_numbers_upto(n)[-1]

def partition_conjugate(partition):
    """
    Return the conjugate partition of 'partition'.

    For example:
    >>> partition_conjugate((5, 3, 1))
    (3, 2, 2, 1, 1)
    """
    conj = []
    last_index = len(partition) - 1
    for i, part in enumerate(partition):
        if i == last_index:
            conj.extend([i + 1] * part)
        elif part != partition[i + 1]:
            conj.extend([i + 1] * (part - partition[i + 1]))
    conj.reverse()
    return tuple(conj)


# aliases
combinationIndexGenerator = combination_index_generator
partitionGenerator = partition_generator
permutationGenerator = permutation_generator
DyckWordGenerator = dyck_word_generator
