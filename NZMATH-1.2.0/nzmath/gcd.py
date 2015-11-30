"""
funtions related to the greatest common divisor of integers.
"""

import nzmath.arygcd as arygcd

def gcd(a, b):
    """
    Return the greatest common divisor of 2 integers a and b.
    """
    a, b = abs(a), abs(b)
    while b:
        a, b = b, a % b
    return a

def binarygcd(a, b):
    """
    Return the greatest common divisor of 2 integers a and b
    by binary gcd algorithm.
    """
    # use arygcd version
    a, b = abs(a), abs(b)
    return arygcd.binarygcd(a, b)

def extgcd(x, y):
    """
    Return a tuple (u, v, d); they are the greatest common divisor d
    of two integers x and y and u, v such that d = x * u + y * v.
    """
    # Crandall & Pomerance "PRIME NUMBERS", Algorithm 2.1.4
    a, b, g, u, v, w = 1, 0, x, 0, 1, y
    while w:
        q, t = divmod(g, w)
        a, b, g, u, v, w = u, v, w, a-q*u, b-q*v, t
    if g >= 0:
        return (a, b, g)
    else:
        return (-a, -b, -g)

def gcd_of_list(integers):
    """
    Return a list [d, [c1, ..., cn]] for a list of integers [x1, ..., xn]
    such that d = c1 * x1 + ... + cn * xn.
    """
    the_gcd = 0
    total_length = len(integers)
    coeffs = []
    coeffs_length = 0
    for integer in integers:
        multiplier, new_coeff, the_gcd = extgcd(the_gcd, integer)
        if multiplier != 1:
            for i in range(coeffs_length):
                coeffs[i] *= multiplier
        coeffs.append(new_coeff)
        coeffs_length += 1
        if the_gcd == 1:
            coeffs.extend([0] * (total_length - coeffs_length))
            break
    return [the_gcd, coeffs]

def lcm(a, b):
    """
    Return the least common multiple of given 2 integers.
    If both are zero, it raises an exception.
    """
    return a // gcd(a, b) * b

def coprime(a, b):
    """
    Return True if a and b are coprime, False otherwise.

    For Example:
    >>> coprime(8, 5)
    True
    >>> coprime(-15, -27)
    False
    >>>
    """
    return gcd(a, b) == 1

def pairwise_coprime(int_list):
    """
    Return True if all integers in int_list are pairwise coprime,
    False otherwise.

    For example:
    >>> pairwise_coprime([1, 2, 3])
    True
    >>> pairwise_coprime([1, 2, 3, 4])
    False
    >>>
    """
    int_iter = iter(int_list)
    product = int_iter.next()
    for n in int_iter:
        if not coprime(product, n):
            return False
        product *= n
    return True
