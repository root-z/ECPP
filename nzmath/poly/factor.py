from __future__ import division
import nzmath.arith1 as arith1
import nzmath.prime as prime
import nzmath.rational as rational
import nzmath.combinatorial as combinatorial
import nzmath.intresidue as intresidue
import nzmath.finitefield as finitefield
import nzmath.poly.uniutil as uniutil
import nzmath.poly.ring as poly_ring
import nzmath.poly.hensel as hensel

def zassenhaus(f):
    """
    zassenhaus(f) -> list of factors of f.

    Factor a squarefree monic integer coefficient polynomial f with
    Berlekamp-Zassenhaus method.
    """
    # keep leading coefficient
    lf = f.leading_coefficient()

    # p-adic factorization
    p, fp_factors = padic_factorization(f)
    if len(fp_factors) == 1:
        return [f]

    # purge leading coefficient from factors
    for i,g in enumerate(fp_factors):
        if g.degree() == 0:
           del fp_factors[i]
           break

    # lift to Mignotte bound
    blm = upper_bound_of_coefficient(f)
    bound = p**(arith1.log(2*blm,p)+1)

    # Hensel lifting
    lf_inv_modq = intresidue.IntegerResidueClass(lf, bound).inverse()
    fq = f.coefficients_map(lambda c: (lf_inv_modq*c).minimumAbsolute()) # fq is monic
    fq_factors, q = hensel.lift_upto(fq, fp_factors, p, bound)

    return brute_force_search(f, fq_factors, bound)

def padic_factorization(f):
    """
    padic_factorization(f) -> p, factors

    Return a prime p and a p-adic factorization of given integer
    coefficient squarefree polynomial f. The result factors have
    integer coefficients, injected from F_p to its minimum absolute
    representation. The prime is chosen to be 1) f is still squarefree
    mod p, 2) the number of factors is not greater than with the
    successive prime.
    """
    num_factors = f.degree()
    stock = None
    for p in prime.generator():
        fmodp = uniutil.polynomial(
            f.terms(),
            finitefield.FinitePrimeField.getInstance(p))
        if f.degree() > fmodp.degree():
            continue
        g = fmodp.getRing().gcd(fmodp,
                                fmodp.differentiate())
        if g.degree() == 0:
            fp_factors = fmodp.factor()
            if not stock or num_factors > len(fp_factors):
                stock = (p, fp_factors)
                if len(fp_factors) == 1:
                    return stock
                num_factors = len(fp_factors)
            else:
                break
    p = stock[0]
    fp_factors = []
    for (fp_factor, m) in stock[1]:
        assert m == 1 # since squarefree
        fp_factors.append(minimum_absolute_injection(fp_factor))
    return (p, fp_factors)

def brute_force_search(f, fp_factors, q):
    """
    brute_force_search(f, fp_factors, q) -> [factors]

    Find the factorization of f by searching a factor which is a
    product of some combination in fp_factors.  The combination is
    searched by brute force.
    """
    factors = []
    d, r = 1, len(fp_factors)
    while 2*d <= r:
        found, combination = find_combination(f, d, fp_factors, q)
        if found:
            factors.append(found)
            for picked in combination:
                fp_factors.remove(picked)
            f = f.pseudo_floordiv(found).primitive_part()
            r -= d
        else:
            d += 1
    factors.append(f.primitive_part())
    return factors

def padic_lift_list(f, factors, p, q):
    """
    padicLift(f, factors, p, q) -> lifted_factors

    Find a lifted integer coefficient polynomials such that:
      f = G1*G2*...*Gm (mod q*p),
      Gi = gi (mod q),
    from f and gi's of integer coefficient polynomials such that:
      f = g1*g2*...*gm (mod q),
      gi's are pairwise coprime
    with positive integers p dividing q.
    """
    ZpZx = poly_ring.PolynomialRing(
        intresidue.IntegerResidueClassRing.getInstance(p))
    gg = arith1.product(factors)
    h = ZpZx.createElement([(d, c // q) for (d, c) in (f - gg).iterterms()])
    lifted = []
    for g in factors:
        gg = gg.pseudo_floordiv(g)
        g_mod = ZpZx.createElement(g)
        if gg.degree() == 0:
            break
        u, v, w = extgcdp(g, gg, p)
        if w.degree() > 0:
            raise ValueError("factors must be pairwise coprime.")
        v_mod = ZpZx.createElement(v)
        t = v_mod * h // g_mod
        lifted.append(g + minimum_absolute_injection(v_mod * h - g_mod * t)*q)
        u_mod = ZpZx.createElement(u)
        gg_mod = ZpZx.createElement(gg)
        h = u_mod * h + gg_mod * t
        lifted.append(g + minimum_absolute_injection(h)*q)
    return lifted

def extgcdp(f, g, p):
    """
    extgcdp(f,g,p) -> u,v,w

    Find u,v,w such that f*u + g*v = w = gcd(f,g) mod p.
    """
    zpz = intresidue.IntegerResidueClassRing.getInstance(p)
    f_zpz = uniutil.polynomial(f, zpz)
    g_zpz = uniutil.polynomial(g, zpz)
    zero, one = f_zpz.getRing().zero, f_zpz.getRing().one
    u, v, w, x, y, z = one, zero, f_zpz, zero, one, g_zpz
    while z:
        q = w // z
        u, v, w, x, y, z = x, y, z, u - q*x, v - q*y, w - q*z
    if w.degree() == 0 and w[0] != zpz.one:
        u = u.scalar_exact_division(w[0]) # u / w
        v = v.scalar_exact_division(w[0]) # v / w
        w = w.getRing().one # w / w
    u = minimum_absolute_injection(u)
    v = minimum_absolute_injection(v)
    w = minimum_absolute_injection(w)
    return u, v, w

def minimum_absolute_injection(f):
    """
    minimum_absolute_injection(f) -> F

    Return an integer coefficient polynomial F by injection of a Z/pZ
    coefficient polynomial f with sending each coefficient to minimum
    absolute representatives.
    """
    coefficientRing = f.getCoefficientRing()
    if isinstance(coefficientRing, intresidue.IntegerResidueClassRing):
        p = coefficientRing.m
    elif isinstance(coefficientRing, finitefield.FinitePrimeField):
        p = coefficientRing.getCharacteristic()
    else:
        raise TypeError("unknown ring (%s)" % repr(coefficientRing))
    half = p >> 1
    g = {}
    for i, c in f.iterterms():
        if c.n > half:
            g[i] = c.n - p
        else:
            g[i] = c.n
    return uniutil.polynomial(g, rational.theIntegerRing)


def upper_bound_of_coefficient(f):
    """
    upper_bound_of_coefficient(polynomial) -> int

    Compute Landau-Mignotte bound of coefficients of factors, whose
    degree is no greater than half of the given polynomial.  The given
    polynomial must have integer coefficients.
    """
    weight = 0
    for c in f.itercoefficients():
        weight += abs(c)**2
    weight = arith1.floorsqrt(weight) + 1
    degree = f.degree()
    lc = f[degree]
    m = (degree >> 1) + 1
    bound = 1
    for i in range(1, m):
        b = combinatorial.binomial(m - 1, i) * weight + \
            combinatorial.binomial(m - 1, i - 1) * lc
        if bound < b:
            bound = b
    return bound

def find_combination(f, d, factors, q):
    """
    find_combination(f, d, factors, q) -> g, list

    Find a combination of d factors which divides f (or its
    complement).  The returned values are: the product g of the
    combination and a list consisting of the combination itself.
    If there is no combination, return (0,[]).
    """
    lf = f.leading_coefficient()
    ZqZX = poly_ring.PolynomialRing(
        intresidue.IntegerResidueClassRing.getInstance(q))

    if d == 1:
        for g in factors:
            product = minimum_absolute_injection(ZqZX.createElement(lf*g))
            if divisibility_test(lf*f, product):
                return (product.primitive_part(), [g])
    else:
        for idx in combinatorial.combinationIndexGenerator(len(factors), d):
            picked = [factors[i] for i in idx]
            product = lf * arith1.product(picked)
            product = minimum_absolute_injection(ZqZX.createElement(product))
            if divisibility_test(lf*f, product):
                return (product.primitive_part(), picked)
    return 0, [] # nothing found

def divisibility_test(f, g):
    """
    Return boolean value whether f is divisible by g or not, for polynomials.
    """
    if g[0] and f[0] % g[0]:
        return False
    if isinstance(f, uniutil.FieldPolynomial) and f % g:
        return False
    elif isinstance(f, uniutil.UniqueFactorizationDomainPolynomial) and f.pseudo_mod(g):
        return False
    return True

def integerpolynomialfactorization(f):
    """
    integerpolynomialfactorization -> list of (factors,index) of f.

    Factor a integer coefficient polynomial f with
    Berlekamp-Zassenhaus method.
    """
    cont = f.content()
    prim = f.primitive_part()
    F = [prim]
    G = prim
    c = 0
    one = G.getRing().one
    while (G.differentiate() and F[c] != one):
        deriv = G.differentiate()
        F.append(F[c].subresultant_gcd(deriv))
        c = c + 1
        G = F[c]
    sqfree_part = F[0].pseudo_floordiv(F[0].subresultant_gcd(F[1])).primitive_part()
    N = zassenhaus(sqfree_part)

    if cont != 1:
        result = [(one.scalar_mul(cont) ,1)]
    else:
        result = []

    F.reverse()
    e = len(F)
    for factor in N:
        for deg, deriv in enumerate(F):
            if not (deriv.pseudo_mod(factor)):
                result.append((factor, (e-deg)))
                break
    return result
