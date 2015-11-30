import random
import nzmath.arith1 as arith1
import nzmath.arygcd as arygcd

def c_root_p(a, p):
    """
    Return the cubic root modulo p.
    (i.e. a^3 = x (mod p))
    """
    if (a % p) == 0:
        return [0]
    if p == 2 or p == 3 :
        return [a % p]
    if (p % 3) == 2:
        return [pow(a, (((2 * p) - 1) // 3), p)]
    p_div_3, p_mod_3 = divmod((p - 1), 3)
    # Compute e,q
    e = 0
    temp = p_div_3
    tempmod = p_mod_3

    while tempmod == 0:
        e += 1
        temp, tempmod = divmod(temp, 3)
    q = (p - 1) // (3 ** e)
    search_range = (p - 1) >> 1
    h = 2
    while pow(h, p_div_3, p) == 1:
        h = random.randrange(2, search_range)
    sym = pow(h, p_div_3, p)
    g = pow(h, q, p)
    # Initialize
    y = g
    r = e
    if q % 3 == 2:
        x = pow(a, (q - 2) // 3, p)
    else:
        x = pow(a, (((2 * q) - 2) // 3), p)

    b = (pow(a, 2, p) * pow(x, 3, p)) % p
    x = (a * x) % p
    while (b % p) != 1:
        # Find exponent
        b_pow = pow(b, 3, p)
        m = 1
        while b_pow != 1:
            b_pow = (b_pow ** 3) % p
            m += 1
        if m == r:
            raise ValueError("there is no cubic root mod p")
        # Reduce exponent
        if sym == pow(b, pow(3, m - 1, p), p):
            t = pow(y, 2, p)
            sym = pow(sym, 2, p)
        else:
            t = y
        t = pow(t, pow(3, r - m - 1, p), p)
        y = pow(t, 3, p)
        r = m
        x = (x * t) % p
        b = (b * y) % p
    return [ x, (x * sym) % p, (x * pow(sym, 2, p)) % p ]

def c_residue(a, p):
    """
    Check whether rational integer a is cubic residue mod prime p.
    If a % p == 0, then return 0,
    elif cubic residue, then return 1,
    elif cubic non-residue, then return -1.
    """
    if p == 3:
        if a % p == 0:
            return 0
        else:
            return 1
    elif p % 3 == 1:
        b1, b2 = decomposite_p(p)
        return c_symbol(a, 0, b1, b2)
    else:
        return c_symbol(a, 0, p, 0)

def c_symbol(a1, a2, b1, b2):
    """
    Return the (Jacobi) cubic residue symbol of 2 Eisenstein Integers ((a1+a2*w)/(b1+b2*w)).
    assume that b1+b2*w is not divisable 1-w.
    """
    r1, r2, j1 = _divides(a1, a2)
    r1, r2, i1 = _FormAdj_w(r1, r2)
    d1, d2, i2 = _FormAdj_w(b1, b2)
    m = (d1 - 1) / 3
    n = d2 / 3
    t = ((m * j1) + (-(m + n) * i1)) % 3
    a1, a2 = r1, r2
    b1, b2 = d1, d2
    nrm_a, nrm_b = arygcd._ap_norm_w(a1, a2, b1, b2)
    if nrm_a < nrm_b:
        tmpre, tmpim = a1, a2
        a1, a2 = b1, b2
        b1, b2 = tmpre, tmpim
        m = (b1 - 1) / 3
        n = b2 / 3
    while not((a1 == b1) and (a2 == b2)):
        sub_re_ab, sub_im_ab = a1 - b1, a2 - b2
        tmp_sub_re_ab, tmp_sub_im_ab, j = _divides(sub_re_ab, sub_im_ab)
        r1, r2, i = _FormAdj_w(tmp_sub_re_ab, tmp_sub_im_ab)
        t = (t + (m * j) + (-(m + n) * i)) % 3
        a1, a2 = r1, r2
        nrm_a, nrm_b = arygcd._ap_norm_w(a1, a2, b1, b2)
        if nrm_a < nrm_b:
            tmpre, tmpim = a1, a2
            a1, a2 = b1, b2
            b1, b2 = tmpre, tmpim
            m = (b1 - 1) / 3
            n = b2 / 3
    if not ((a1 == 1) and (a2 == 0)) :
        return 0
    elif t  == 0:
        return 1
    else:
        return -1

def _FormAdj_w(a1, a2):
    """
    Transform Eisenstein Integer a1+a2*w  ->  (-w)^i * (x+y*w)
    x+y*w is primary element
    assume that a1+a2*w is not divisible 1-w
    """
    if a1 % 3 == 0:
        if ((a2 % 3) == 2) or ((a2 % 3) == -1):
            return a1 - a2, a1, 1
        elif ((a2 % 3) == 1) or ((a2 % 3) == -2):
            return a2 - a1, -a1, 4
    elif ((a1 % 3) == 1) or ((a1 % 3) == -2):
        if ((a2 % 3) == 1) or ((a2 % 3) == -2):
            return a2, a2 - a1, 5
        elif a2 % 3 == 0:
            return a1, a2, 0
    else:
        if ((a2 % 3) == 2) or ((a2 % 3) == -1):
            return -a2, a1 - a2, 2
        elif a2 % 3 == 0:
            return -a1, -a2, 3

def _divides(a1, a2):
    """
    divide 1-w
    (i.e. a1+a2*w -> (1-w)^k * (x+y*w))
    """
    if (a1 == 0) and (a2 == 0):
        return a1, a2, 0
    j = 0
    while (a1 + a2) % 3 == 0:
        tmpa1 = a1
        a1 = ((a1 + a1) - a2) / 3
        a2 = (tmpa1 + a2) / 3
        j += 1
    return a1, a2, j

def decomposite_p(p):
    """
    Decomposite p = 1(mod 3) into its primary prime factors in z[w]
    """
    s, t = cornacchia(3, p)
    if t % 3 == 0:
        a = s + t
        b = 2 * t
    if (t % 3) != 0:
        a = 2 * t
        b = s + t
    return a, b

def cornacchia(d, p):
    """
    Return the solution of x^2 + d * y^2 = p .
    p be a prime and d be an integer such that 0 < d < p.
    """
    if (d <= 0) or (d >= p):
        raise ValueError("invalid input")
    k = arith1.legendre(-d, p)
    if k == -1:
        raise ValueError("no solution")
    x0 = arith1.modsqrt(-d, p)
    if x0 < (p / 2):
        x0 = p - x0
    a = p
    b = x0
    l = arith1.floorsqrt(p)
    while b > l:
        a, b = b, a % b
    c, r = divmod(p - b * b, d)
    if r:
        raise ValueError("no solution")
    t = arith1.issquare(c)
    if t == 0:
        raise ValueError("no solution")
    else:
        return (b, t)
