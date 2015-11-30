"""
prime ideal decomposition over an (absolute) number field
"""

import nzmath.arith1 as arith1
import nzmath.rational as rational
import nzmath.matrix as matrix
import nzmath.vector as vector
import nzmath.algfield as algfield
import nzmath.finitefield as finitefield
import nzmath.poly.uniutil as uniutil
import nzmath.module as module
import nzmath.algorithm as algorithm


def prime_decomp(p, polynomial):
    """
    Return prime decomposition of (p) over Q[x]/(polynomial).
    This method returns a list of (P_k, e_k, f_k), 
    where P_k is an instance of Ideal_with_generator, 
          e_k is the ramification index of P_k, 
          f_k is the residue degree of P_k.
    """
    field = algfield.NumberField(polynomial)
    field_disc = field.field_discriminant()
    if (field.disc() // field_disc) % p != 0:
        return _easy_prime_decomp(p, polynomial)
    else:
        return _main_prime_decomp(p, polynomial)

def _easy_prime_decomp(p, polynomial):
    """
    prime decomposition by factoring polynomial
    """
    f = algfield.fppoly(polynomial, p)
    d = f.degree()
    factor = f.factor()
    p_alg = algfield.BasicAlgNumber([[p] + [0] * (d - 1), 1], polynomial)
    if len(factor) == 1 and factor[0][1] == 1:
        return [(module.Ideal_with_generator([p_alg]),
            1, factor[0][0].degree())]
    else:
        dlist = []
        for i in range(len(factor)):
            basis_list = []
            for j in range(d):
                if factor[i][0][j] == 0:
                    basis_list.append(0)
                else:
                    basis_list.append(factor[i][0][j].toInteger())
            prime_ideal = module.Ideal_with_generator([p_alg,
                          algfield.BasicAlgNumber([basis_list, 1], polynomial)])
            dlist.append((prime_ideal, factor[i][1], factor[i][0].degree()))
        return dlist

def _main_prime_decomp(p, polynomial):
    """
    main step of prime decomposition
    """
    field = algfield.NumberField(polynomial)
    n = field.degree
    int_basis = field.integer_ring()
    base_multiply = module._base_multiplication(int_basis, field)
    H_lst = _squarefree_decomp(p, field, base_multiply) # step 2-5
    done_list = []
    for j in range(len(H_lst)):
        H_j = H_lst[j]
        # step 7
        if H_j != None and H_j.column == n: #rank(H_j)
            continue
        L = [H_j]
        while L:
            H = L.pop()
            A, gamma_multiply = _separable_algebra(
                p, H, base_multiply) # step 8, 9
            sol = _check_H(p, H, field, gamma_multiply) # step 10-11
            if isinstance(sol, tuple):
                done_list.append((sol[0], j + 1, sol[1]))
            else:
                L += _splitting_squarefree(p, H, base_multiply, 
                    A, gamma_multiply, sol) # step 12-15
    return done_list

def _squarefree_decomp(p, field, base_multiply):
    """
    decompose (p) as (p)= A_1 (A_2)^2 (A_3)^2 ..., where each A_i is squarefree
    """
    #CCANT Algo 6.2.9 step 2-5
    n = field.degree
    int_basis = field.integer_ring()
    # step 2
    log_n_p = arith1.log(n, p)
    q = p ** log_n_p
    if  q < n:
        q *= p
    base_p = _kernel_of_qpow(int_basis, q, p, field)
    I_p_over_pO = base_p
    # step 3
    K_lst = [I_p_over_pO]
    while K_lst[-1]: # K[-1] is zero matrix in F_p
        K_lst.append(_mul_mod_pO(K_lst[0], K_lst[-1], base_multiply))
    i = len(K_lst)
    # step 4
    J_lst = [K_lst[0]]
    for j in range(1, i):
        J_lst.append(_div_mod_pO(K_lst[j], K_lst[j - 1], base_multiply, p))
    # step 5
    H_lst = []
    for j in range(i - 1):
        H_lst.append(_div_mod_pO(J_lst[j], J_lst[j + 1], base_multiply, p))
    H_lst.append(J_lst[-1])
    return H_lst

def _kernel_of_qpow(basis_matrix, q, p, field):
    """
    compute kernel of w^q.
    (this method is same as round2._kernel_of_qpow)
    """
    omega = _matrix_to_algnumber_list(basis_matrix, field)
    omega_pow = []
    for omega_i in omega:
        omega_pow.append(omega_i ** q)
    omega_pow_mat = _algnumber_list_to_matrix(omega_pow, field)
    A_over_Z = basis_matrix.inverse(omega_pow_mat) #theta repr -> omega repr
    A = A_over_Z.map(
         lambda x: finitefield.FinitePrimeFieldElement(x.numerator, p))
    return A.kernel()

def _mul_mod_pO(self, other, base_multiply):
    """
    ideal multiplication modulo pO using base_multiply
    """
    #CCANT Algorithm 6.2.5
    n = base_multiply.column
    new_self = _matrix_mul_by_base_mul(self, other, base_multiply)
    return new_self.image()

def _div_mod_pO(self, other, base_multiply, p):
    """
    ideal division modulo pO by using base_multiply
    """
    n = other.row
    m = other.column 
    F_p = finitefield.FinitePrimeField(p)
    #Algo 2.3.7
    if self == None:
        self_new = matrix.zeroMatrix(n, F_p)
        r = 0
    else:
        self_new = self
        r = self.column
    if r == m:
        return matrix.unitMatrix(n, F_p)
    # if r > m, raise NoInverseImage error
    X = matrix.Subspace.fromMatrix(other.inverseImage(self_new))
    B = X.supplementBasis()
    other_self = other * B # first r columns are self, others are supplement
    other_self.toFieldMatrix()

    gamma_part = other_self.subMatrix(
        range(1, n + 1), range(r + 1, m + 1))
    omega_gamma = other_self.inverseImage(_matrix_mul_by_base_mul(
        matrix.unitMatrix(n, F_p), gamma_part, base_multiply))
    vect_list = []
    for k in range(1, n + 1):
        vect = vector.Vector([F_p.zero] * ((m - r) ** 2))
        for i in range(1, m - r + 1):
            for j in range(1, m - r + 1):
                vect[(m - r) * (i - 1) + j] = _mul_place(
                    k, i, omega_gamma, m - r)[j + r]
        vect_list.append(vect)
    return matrix.FieldMatrix( (m - r) ** 2, n, vect_list).kernel()

def _separable_algebra(p, H, base_multiply):
    """
    return A=O/H (as matrix representation)
    """
    # CCANT Algo 6.2.9 step 8-9
    F_p = finitefield.FinitePrimeField(p)
    n = base_multiply.row
    if H == None:
        H_new = matrix.zeroMatrix(n, 1, F_p)
        r = 0
    else:
        H_new = H
        r = H.column
    # step 8
    H_add_one = H_new.copy()
    if H_add_one.column == n:
        raise ValueError
    H_add_one.extendColumn(vector.Vector([F_p.one] + [F_p.zero] * (n - 1)))
    H_add_one.toSubspace(True)
    full_base = H_add_one.supplementBasis()
    full_base.toFieldMatrix()
    # step 9
    A = full_base.subMatrix(
        range(1, n + 1), range(r + 1, n + 1))
    gamma_mul = full_base.inverseImage(
        _matrix_mul_by_base_mul(A, A, base_multiply))
    gamma_mul = gamma_mul.subMatrix(
        range(r + 1, n + 1), range(1, gamma_mul.column + 1))
    return A, gamma_mul

def _check_H(p, H, field, gamma_mul):
    """
    If H express the prime ideal, return (H, f),
    where f: residual degree.
    Else return some column of M_1 (not (1,0,..,0)^t).
    """
    F_p = finitefield.FinitePrimeField(p)
    if H == None:
        f = field.degree
    else:
        f = H.row - H.column # rank(A), A.column
    # CCANT Algo 6.2.9 step 10-11
    # step 10
    B = matrix.unitMatrix(f, F_p)
    M_basis = []
    for j in range(1, f + 1):
        alpha_pow = _pow_by_base_mul(B[j], p, gamma_mul, f)
        alpha_pow -= B[j]
        M_basis.append(alpha_pow)
    M_1 = matrix.FieldMatrix(f, f, M_basis).kernel()
    # step 11
    if M_1.column > 1: # dim(M_1) > 1
        return M_1[M_1.column]
    else:
        H_simple = _two_element_repr_prime_ideal(
                   p, H.map(lambda x: x.getResidue()), field, f)
        p_alg = algfield.BasicAlgNumber([[p] + [0] * (field.degree - 1), 1], 
                field.polynomial)
        sol = module.Ideal_with_generator([p_alg, H_simple])
        return (sol, f)

def _two_element_repr_prime_ideal(p, hnf, field, f):
    """
    return alpha such that (p, alpha) is same ideal to 
    the ideal represented by hnf (matrix) w.r.t. the integer ring of field.
    we assume the ideal is a prime ideal whose norm is p^f, where p is prime.
    """
    # assume that first basis of the integral basis equals 1
    int_basis = field.integer_ring()
    
    # step 1
    n = field.degree
    Z = rational.theIntegerRing
    beta = (p * matrix.unitMatrix(n, Z).subMatrix(
        range(1, n + 1), range(2, n + 1))).copy()
    beta.extendColumn(hnf)
    m = beta.column
    # step 2
    R = 0
    P_norm = p ** f
    while True:
        R += 1
        lmd = [R] * m
        while True:
            # step 3
            alpha = lmd[0] * beta[1]
            for j in range(1, m):
                alpha += lmd[j] * beta[j + 1]
            alpha = _vector_to_algnumber(int_basis * alpha, field)
            n = alpha.norm() / P_norm
            r = divmod(int(n), p)[1]
            if r:
                return alpha
            # step 4
            for j in range(m - 1, 0, -1):
                if lmd[j] != -R:
                    break
            lmd[j] -= 1
            for i in range(j + 1, m):
                lmd[i] = R
            # step 5
            for j in range(m):
                if lmd[j] != 0:
                    break
            else:
                break

def _splitting_squarefree(p, H, base_multiply, A, gamma_mul, alpha):
    """
    split squarefree part
    """
    F_p = finitefield.FinitePrimeField(p)
    if H == None:
        n = base_multiply.row
        f = n
        H_new = matrix.zeroMatrix(n, 1, F_p)
    else:
        n = H.row
        f = n - H.column
        H_new = H.copy()
    # step 12
    one_gamma = vector.Vector([F_p.one] + [F_p.zero] * (f - 1)) # w.r.t. gamma_1
    minpoly_matrix_alpha = matrix.FieldMatrix(f, 2, [one_gamma, alpha])
    ker = minpoly_matrix_alpha.kernel()
    alpha_pow = alpha
    while ker == None:
        alpha_pow = _vect_mul_by_base_mul(alpha_pow, alpha, gamma_mul, f)
        minpoly_matrix_alpha.extendColumn(alpha_pow)
        ker = minpoly_matrix_alpha.kernel()
    minpoly_alpha = uniutil.FinitePrimeFieldPolynomial(
        enumerate(ker[1].compo), F_p)
    # step 13
    m_list = minpoly_alpha.factor()
    # step 14
    L = []
    for (m_s, i) in m_list:
        M_s = H_new.copy()
        beta_s_gamma = _substitution_by_base_mul(
            m_s, alpha, gamma_mul, one_gamma, f)
        # beta_s_gamma (w.r.t. gamma) -> beta_s (w.r.t. omega)
        beta_s = A * beta_s_gamma
        omega_beta_s = _matrix_mul_by_base_mul(
            matrix.unitMatrix(n, F_p), beta_s.toMatrix(True), base_multiply)
        M_s.extendColumn(omega_beta_s)
        L.append(M_s.image())
    return L

#################################################
# Utility function
#################################################
def _vector_to_algnumber(vector_repr, field):
    """
    vector (w.r.t. theta) -> omega
    """
    int_vector, denom = module._toIntegerMatrix(vector_repr)
    int_vector = int_vector[1].compo # Submodule -> list
    omega = algfield.BasicAlgNumber(
            [int_vector, denom], field.polynomial)
    return omega

def _algnumber_to_vector(omega, field):
    """
    omega -> vector (w.r.t. theta)
    """
    n = field.degree
    vect = [rational.Rational(
        omega.coeff[j], omega.denom) for j in range(n)]
    return vector.Vector(vect)

def _matrix_to_algnumber_list(matrix_repr, field):
    """
    matrix (w.r.t. theta) -> [omega_i]
    """
    int_matrix, denom = module._toIntegerMatrix(matrix_repr)
    omega = []
    for i in range(1, int_matrix.column + 1):
        omega_i = algfield.BasicAlgNumber(
            [int_matrix[i], denom], field.polynomial)
        omega.append(omega_i)
    return omega

def _algnumber_list_to_matrix(omega, field):
    """
    [omega_i] -> matrix (w.r.t. theta)
    """
    vect = []
    for omega_i in omega:
        vect_i = _algnumber_to_vector(omega_i, field)
        vect.append(vector.Vector(vect_i))
    return matrix.Matrix(omega[0].degree, len(omega), vect)

def _vect_mul_by_base_mul(self, other, base_mul, max_j=False):
    """
    return self * other with respect to basis (gamma_i),
    where self, other are vector.
    (base_mul)_ij express gamma_i gamma_j
    """
    n = base_mul.row
    sol_val = vector.Vector([0] * n)
    for i1 in range(1, n + 1):
        for i2 in range(1, n + 1):
            if not max_j:
                sol_val += self[i1] * other[
                    i2] * module._symmetric_element(i1, i2, base_mul)
            else:
                sol_val += self[i1] * other[
                    i2] * _mul_place(i1, i2, base_mul, max_j)
    return sol_val

def _matrix_mul_by_base_mul(self, other, base_mul, max_j=False):
    """
    return self * other with respect to basis (gamma_i),
    where self, other are matrix.
    (base_mul)_ij express gamma_i gamma_j
    """
    n = base_mul.row
    sol = []
    for k1 in range(1, self.column + 1):
        for k2 in range(1, other.column + 1):
            sol_ele = _vect_mul_by_base_mul(
                self[k1], other[k2], base_mul, max_j)
            sol.append(sol_ele)
    return matrix.FieldMatrix(n, len(sol), sol)

def _mul_place(i, j, mul_matrix, max_j):
    """
    find the vector in mul_matrix corresponding (i, j),
    where 1 <= j <= max_j.
    """
    return mul_matrix[(i - 1) * max_j + j]

def _pow_by_base_mul(self, other, base_mul, max_j=False):
    """
    return self * other with respect to basis (gamma_i),
    where self, other are matrix.
    (base_mul)_ij express gamma_i * gamma_j
    This function uses right-left binary method.
    """
    func = algorithm.powering_func(
        lambda x,y:_vect_mul_by_base_mul(x, y, base_mul, max_j))
    return func(self, other)

def _substitution_by_base_mul(P, ele, base_mul, one_gamma, max_j=False):
    """
    return P(ele) with respect to basis (gamma_i),
    where P is a polynomial, ele is a vector w.r.t. gamma_i.
    (base_mul)_ij express gamma_i * gamma_j
    This function uses digital method (Horner method).
    """
    coeff = reversed(P.sorted)
    mul = lambda x, y : _vect_mul_by_base_mul(x, y, base_mul, max_j)
    return algorithm.digital_method(
        coeff, ele, lambda x, y : x + y, mul,
        lambda a, x : a * x,
        algorithm.powering_func(one_gamma, mul),
        vector.Vector([base_mul.coeff_ring.zero] * base_mul.row),
        one_gamma
        )
