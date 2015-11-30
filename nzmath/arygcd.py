"""
binary-like gcd algorithms for rational, gauss, and eisenstein integer
"""

import math


def bit_num(a):
    """
    return the number of bits
    """
    if a == 0:
        return 1
    else:
        return int(math.log(a, 2)) + 1


def binarygcd(a, b):
    """
    Return the greatest common divisor of 2 integers a and b
    by binary gcd algorithm.
    """
    if a < b:
        a, b = b, a
    if b == 0:
        return a
    a, b = b, a % b
    if b == 0:
        return a
    k = 0
    while not a&1 and not b&1:
        k += 1
        a >>= 1
        b >>= 1
    while not a&1:
        a >>= 1
    while not b&1:
        b >>= 1
    while a:
        while not a & 1:
            a >>= 1
        if abs(a) < abs(b):
            a, b = b, a
        a = (a - b) >> 1
    return b << k


def arygcd_i(a1, a2, b1, b2):
    """
    Return the greatest common divisor of 2 gauss-integers a1+a2*i and b1+b2*i
    by (1+i)-ary gcd algorithm.
    """
    if a1 == 0 and a2 == 0:
       return b1, b2
    elif b1 == 0 and b2 == 0:
       return a1, a2
    ap, bp = 0, 0
    while (a1-a2)&1 == 0:
        a1, a2 = (a1+a2) >> 1, (a2-a1) >> 1
        ap += 1
    while (b1-b2)&1 == 0:
        b1, b2 = (b1+b2) >> 1, (b2-b1) >> 1
        bp += 1
    k = min(ap, bp)
    while a1 != 0 or a2 != 0:
        while (a1-a2) & 1 == 0:
            a1, a2 = (a1+a2) >> 1, (a2-a1) >> 1
        norm_a, norm_b = _ap_norm_i(a1, a2, b1, b2)
        if norm_a < norm_b:
            a1, a2, b1, b2 = b1, b2, a1, a2
        a1, a2 = _FormAdj_i(a1, a2)
        b1, b2 = _FormAdj_i(b1, b2)
        a1, a2 = a1-b1, a2-b2
    if 0 in (ap, bp):
        return b1, b2
    else:
        s, t = _n_pow_i(1, 1, k)
        return (b1*s)-(b2*t), (b1*t)+(b2*s)

def _ap_norm_i(a, b, c, d):
    """
    Return approximately norm of 2 gaussInteger a+b*i and c+d*i  
    """
    a, b, c, d = abs(a), abs(b), abs(c), abs(d)
    max_dig = max(bit_num(a), bit_num(b), bit_num(c), bit_num(d))
    if max_dig > 6:
        sft_num = max_dig - 6
        a, b, c, d = a >> sft_num, b >> sft_num, c >> sft_num, d >> sft_num
        return a*a + b*b, c*c + d*d
    else:
        return a*a + b*b, c*c + d*d

def _n_pow_i(a, b, n):
    """
    return (1+i)**k 
    """
    x = a
    y = b

    for i in range(1, n):
        x1 = (x*a) - (y*b)
        y1 = (y*a) + (x*b)
        x = x1
        y = y1
    return x, y

def _FormAdj_i(a, b):
    """
    transform gaussInteger a+b*i ->  form 1+2(1+i)*(x+y*i)
    """
    if a & 1 == 0:
        a, b = -b, a
    if (b - a + 1) & 3 == 0:
        return a, b
    else:
        return -a, -b


def arygcd_w(a1, a2, b1, b2):
    """
    Return the greatest common divisor of 2 eisensteinIntegers a1+a2*w and b1+b2*w
    by (1-w)-ary gcd algorithm.
    """
    if a1 == 0 and a2 == 0:
       return b1, b2
    elif b1 == 0 and b2 == 0:
       return a1, a2
    ap, bp = 0, 0
    while (a1 + a2) % 3 == 0:
        a1, a2 = (a1 + a1 - a2) / 3, (a1 + a2) / 3
        ap += 1
    while (b1 + b2) % 3 == 0:
        b1, b2 = (b1 + b1 - b2) / 3, (b1 + b2) / 3
        bp += 1
    k = min(ap, bp)
    while a1 != 0 or a2 != 0:
        while (a1 + a2) % 3 == 0:
            a1, a2 = (a1 + a1 - a2) / 3, (a1 + a2) / 3
        nrm_a, nrm_b = _ap_norm_w(a1, a2, b1, b2)
        if nrm_a < nrm_b:
            a1, a2, b1, b2 = b1, b2, a1, a2
        a1, a2 = _FormAdj_w(a1, a2)
        b1, b2 = _FormAdj_w(b1, b2)
        a1, a2 = a1 - b1, a2 - b2
    k1, k2 = _n_pow_w(k)
    return k1*b1 - k2*b2, k1*b2 + k2*b1 - k2*b2

def _ap_norm_w(a, b, c, d):
    """
    Return approximately norm of 2 eisensteinInteger a+b*w and c+d*w  
    """
    a, b, c, d = abs(a), abs(b), abs(c), abs(d)
    max_dig = max(bit_num(a), bit_num(b), bit_num(c), bit_num(d))
    if max_dig > 6:
        sft_num = max_dig - 6
        a, b, c, d = a >> sft_num, b >> sft_num, c >> sft_num, d >> sft_num
        subst1, subst2 = (a - b) >> sft_num, (c - d) >> sft_num
        return a*a + b*b + subst1*subst1, c*c + d*d + subst2*subst2
    else:
        return a*a + b*b + (a - b)*(a - b), c*c + d*d + (c -d)*(c -d)

def _n_pow_w(n):
    """
    return (1-w)**k
    """
    x, y = divmod(n, 2)
    k1 = 3**x
    if y == 1:
        k2 = -k1
    else:
        k2 = 0
    return k1, k2

def _FormAdj_w(a1, a2):
    """
    transform eisensteinInteger a1+a2*w ->  form 1+3*(x+y*w)
    """
    if a1 % 3 == 0:
        if a2 % 3 == -1 or a2 % 3 == 2:
            return a1 - a2, a1
        else:
            return a2 - a1, -a1
    elif a1 % 3 == 1:
        if a2 % 3 == 1:
            return a2, a2 -a1
        else:
            return a1, a2
    else:
        if a2 % 3 == -1 or a2 % 3 == 2:
            return -a2, a1 -a2
        else:
            return -a1, -a2
