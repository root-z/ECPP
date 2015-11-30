import math
import time
import logging
import nzmath.arith1 as arith1
import nzmath.gcd as gcd
import nzmath.prime as prime
import nzmath.factor.find as find

# TODO: add psyco to config entry
##         try:
##             import psyco
##             psyco.full()
##         except ImportError:
##             pass

_log = logging.getLogger('nzmath.factor.mpqs')

class QS(object):
    def __init__(self, n, sieverange, factorbase):
        self.number = n
        self.sqrt_n = int(math.sqrt(n))
        for i in [2,3,5,7,11,17,19]:
            if n % i == 0:
                raise ValueError("This number is divided by %d" % i)

        self.digit = arith1.log(self.number, 10) + 1
        self.Srange = sieverange
        self.FBN = factorbase

        self.move_range = range(self.sqrt_n-self.Srange, self.sqrt_n+self.Srange+1)

        i = 0
        k = 0
        factor_base = [-1]
        FB_log = [0]
        while True:
            ii = primes_table[i]
            if arith1.legendre(self.number, ii) == 1:
                factor_base.append(ii)
                FB_log.append(primes_log_table[i])
                k += 1
                i += 1
                if k == self.FBN:
                    break
            else:
                i += 1

        self.FB = factor_base
        self.FB_log = FB_log
        self.maxFB = factor_base[-1]
        N_sqrt_list = []
        for i in self.FB:
            if i != 2 and i != -1:
                e = int(math.log(2*self.Srange, i))
                N_sqrt_modp = sqroot_power(self.number, i, e)
                N_sqrt_list.append(N_sqrt_modp)
        self.solution = N_sqrt_list  #This is square roots of N on Z/pZ, p in factor base.

        poly_table = []
        log_poly = []
        minus_val = []
        for j in self.move_range:
            jj = (j**2)-self.number
            if jj < 0:
                jj = -jj
                minus_val.append(j-self.sqrt_n+self.Srange)
            elif jj == 0:
                jj = 1
            lj = int((math.log(jj)*30)*0.97)  # 0.97 is an erroe
            poly_table.append(jj)
            log_poly.append(lj)
        self.poly_table = poly_table  # This is Q(x) value , x in [-M+sqrt_n,M+sqrt_n].
        self.log_poly = log_poly      # This is log(Q(x)) value.
        self.minus_check = minus_val # This is "x" that Q(x) is minus value.

    def run_sieve(self):
        T = time.time()
        M = self.Srange
        start_location = []
        logp = [0]*(2*M+1)
        j = 2
        for i in self.solution:
            k = 0
            start_p = []
            while k < len(i):
                q = int((self.sqrt_n)/(self.FB[j]**(k+1)))
                s_1 = q*(self.FB[j]**(k+1))+i[k][0]
                s_2 = q*(self.FB[j]**(k+1))+i[k][1]
                while True:
                    if s_1 < self.sqrt_n-M:
                        s_1 = s_1+(self.FB[j]**(k+1))
                        break
                    else:
                        s_1 = s_1-(self.FB[j]**(k+1))
                while True:
                    if s_2 < self.sqrt_n-M:
                        s_2 = s_2+(self.FB[j]**(k+1))
                        break
                    else:
                        s_2 = s_2-(self.FB[j]**(k+1))
                start_p.append([s_1-self.sqrt_n+M,s_2-self.sqrt_n+M])

                k += 1
            start_location.append(start_p)
            j += 1
        self.start_location = start_location

        if self.poly_table[0] & 1 == 0:
            i = 0
            while i <= 2*M:
                j = 1
                while True:
                    if self.poly_table[i] % (2**(j+1)) == 0:
                        j+=1
                    else:
                        break
                logp[i] += self.FB_log[1]*j
                i += 2
        else:
            i = 1
            while i <= 2*M:
                j = 1
                while True:
                    if self.poly_table[i] % (2**(j+1)) == 0:
                        j += 1
                    else:
                        break
                logp[i] += self.FB_log[1]*j
                i += 2
        L = 2
        for j in self.start_location:
            k = 0
            while k < len(j):
                s_1 = j[k][0]
                s_2 = j[k][1]
                h_1 = 0
                h_2 = 0
                while s_1+h_1 <= 2*M:
                    logp[s_1+h_1] += self.FB_log[L]
                    h_1 += self.FB[L]**(k+1)
                while s_2+h_2 <= 2*M:
                    logp[s_2+h_2] += self.FB_log[L]
                    h_2 += self.FB[L]**(k+1)
                k += 1
            L += 1

        self.logp = logp
        smooth = []
        for t in range(2*M+1):
            if logp[t] >= self.log_poly[t]:
                poly_val = self.poly_table[t]
                index_vector = []
                for p in self.FB:
                    if p == -1:
                        if t in self.minus_check:
                            index_vector.append(1)
                        else:
                            index_vector.append(0)
                    else:
                        r = 0
                        while poly_val % (p**(r+1)) == 0:
                            r += 1
                        v = r & 1
                        index_vector.append(v)
                smooth.append([index_vector, (poly_val, t+self.sqrt_n-M)])
        _log.info(" Sieve time is %f sec" % (time.time()-T))
        _log.info(" Found smooth numbers are %d / %d" % (len(smooth), len(self.FB)))
        self.smooth = smooth
        return smooth


class MPQS(object):
    def __init__(self, n, sieverange=0, factorbase=0, multiplier=0):
        self.number = n
        _log.info("The number is %d MPQS starting" % n)

        if prime.primeq(self.number):
            raise ValueError("This number is Prime Number")
        for i in [2,3,5,7,11,13]:
            if n % i == 0:
                raise ValueError("This number is divided by %d" % i)

        self.sievingtime = 0
        self.coefficienttime = 0
        self.d_list = []
        self.a_list = []
        self.b_list = []

        #Decide prameters for each digits
        self.digit = arith1.log(self.number, 10) + 1

        if sieverange != 0:
            self.Srange = sieverange
            if factorbase != 0:
                self.FBN = factorbase
            elif self.digit < 9:
                self.FBN = parameters_for_mpqs[0][1]
            else:
                self.FBN = parameters_for_mpqs[self.digit-9][1]
        elif factorbase != 0:
            self.FBN = factorbase
            if self.digit < 9:
                self.Srange = parameters_for_mpqs[0][0]
            else:
                self.Srange = parameters_for_mpqs[self.digit-9][0]
        elif self.digit < 9:
            self.Srange = parameters_for_mpqs[0][0]
            self.FBN = parameters_for_mpqs[0][1]
        elif self.digit > 53:
            self.Srange = parameters_for_mpqs[44][0]
            self.FBN = parameters_for_mpqs[44][1]
        else:
            self.Srange = parameters_for_mpqs[self.digit-9][0]
            self.FBN = parameters_for_mpqs[self.digit-9][1]

        self.move_range = range(-self.Srange, self.Srange+1)

        # Decide k such that k*n = 1 (mod4) and k*n has many factor base
        if multiplier == 0:
            self.sqrt_state = []
            for i in [3,5,7,11,13]:
                s = arith1.legendre(self.number, i)
                self.sqrt_state.append(s)

            if self.number % 8 == 1 and self.sqrt_state == [1,1,1,1,1]:
                k=1
            else:
                index8 = (self.number & 7) >> 1
                j = 0
                while self.sqrt_state != prime_8[index8][j][1]:
                    j += 1
                k = prime_8[index8][j][0]
        else:
            if n & 3 == 1:
                k = 1
            else:
                if multiplier == 1:
                    raise ValueError("This number is 3 mod 4 ")
                else:
                    k = multiplier
        self.number = k*self.number
        self.multiplier = k

        _log.info("%d - digits Number" % self.digit)
        _log.info("Multiplier is %d" % self.multiplier)

        # Table of (log p) , p in FB
        i = 0
        k = 0
        factor_base = [-1]
        FB_log = [0]
        while k < self.FBN:
            ii = primes_table[i]
            if arith1.legendre(self.number,ii) == 1:
                factor_base.append(ii)
                FB_log.append(primes_log_table[i])
                k += 1
            i += 1

        self.FB = factor_base
        self.FB_log = FB_log
        self.maxFB = factor_base[-1]

        # Solve x^2 = n (mod p^e)
        N_sqrt_list = []
        for i in self.FB:
            if i != 2 and i != -1:
                e = int(math.log(2*self.Srange, i))
                N_sqrt_modp = sqroot_power(self.number, i, e)
                N_sqrt_list.append(N_sqrt_modp)
        self.Nsqrt = N_sqrt_list

    def make_poly(self):
        """
        Make coefficients of f(x)= ax^2+b*x+c
        """
        T = time.time()
        if self.d_list == []:
            d = int(math.sqrt((math.sqrt(self.number)/(math.sqrt(2)*self.Srange))))
            if d& 1 == 0:
                if (d+1)& 3 == 1: #case d=0 mod4
                    d += 3
                else:
                    d += 1       #case d=2 mod4
            elif d& 3 == 1:       #case d=1 mod4
                d += 2
                                 #case d=3 mod4
        else:
            d = self.d_list[-1]

        while d in self.d_list or not prime.primeq(d) or \
              arith1.legendre(self.number,d) != 1 or d in self.FB:
            d += 4
        a = d**2
        h_0 = pow(self.number, (d-3) >> 2, d)
        h_1 = (h_0*self.number) % d
        h_2 = ((arith1.inverse(2,d)*h_0*(self.number - h_1**2))//d) % d
        b = (h_1 + h_2*d) % a
        if b& 1 == 0:
            b = b - a

        self.d_list.append(d)
        self.a_list.append(a)
        self.b_list.append(b)

        # Get solution of  F(x) = 0 (mod p^i)
        solution = []
        i = 0
        for s in self.Nsqrt:
            k = 0
            p_solution = []
            ppow = 1
            while k < len(s):
                ppow *= self.FB[i+2]
                a_inverse = arith1.inverse(2*self.a_list[-1], ppow)
                x_1 = ((-b + s[k][0])*a_inverse) % ppow
                x_2 = ((-b + s[k][1])*a_inverse) % ppow
                p_solution.append([x_1, x_2])
                k += 1
            i += 1
            solution.append(p_solution)
        self.solution = solution
        self.coefficienttime += time.time() - T

    def run_sieve(self):
        self.make_poly()
        T = time.time()
        M = self.Srange
        a = self.a_list[-1]            #
        b = self.b_list[-1]            # These are coefficients of F(x)
        c = (b**2-self.number)//(4*a)   #
        d = self.d_list[-1]            #

        self.poly_table = []  # This is F(x) value , x in [-M,M].
        self.log_poly = []    # This is log(F(x)) value.
        self.minus_check = [] # This is "x" that F(x) is minus value.
        for j in self.move_range:
            jj = (a * j + b) * j + c
            if jj < 0:
                jj = -jj
                self.minus_check.append(j + M)
            elif jj == 0:
                jj = 1
            lj = int((math.log(jj)*30)*0.95)  # 0.95 is an erroe
            self.poly_table.append(jj)
            self.log_poly.append(lj)

        y = arith1.inverse(2*d, self.number)
        start_location = []
        logp = [0]*(2*M+1)
        j = 2
        for i in self.solution:
            start_p = []
            ppow = 1
            for k in range(len(i)):
                ppow *= self.FB[j]
                q = -M // ppow
                s_1 = (q + 1) * ppow + i[k][0]
                s_2 = (q + 1) * ppow + i[k][1]
                while s_1 + M >= ppow:
                    s_1 = s_1 - ppow
                while s_2 + M >= ppow:
                    s_2 = s_2 - ppow
                start_p.append([s_1+M, s_2+M])
            start_location.append(start_p)
            j += 1
        self.start_location = start_location

        i = self.poly_table[0] & 1
        while i <= 2*M:
            j = 1
            while self.poly_table[i] % (2**(j+1)) == 0:
                j += 1
            logp[i] += self.FB_log[1]*j
            i += 2
        L = 2
        for plocation in self.start_location:
            for k in range(len(plocation)):
                s_1 = plocation[k][0]
                s_2 = plocation[k][1]
                ppow = self.FB[L]**(k+1)
                while s_1 <= 2*M:
                    logp[s_1] += self.FB_log[L]
                    s_1 += ppow
                while s_2 <= 2*M:
                    logp[s_2] += self.FB_log[L]
                    s_2 += ppow
            L += 1

        self.logp = logp
        smooth = []
        for t in range(2*M+1):
            if logp[t] >= self.log_poly[t]:
                poly_val = self.poly_table[t]
                index_vector = []
                H = (y*(2*a*(t-self.Srange)+b))%self.number
                for p in self.FB:
                    if p == -1:
                        if t in self.minus_check:
                            index_vector.append(1)
                        else:
                            index_vector.append(0)
                    else:
                        r = 0
                        while poly_val % (p**(r+1)) == 0:
                            r += 1
                        v = r & 1
                        index_vector.append(v)
                smooth.append([index_vector, (poly_val, H)])
        self.sievingtime += time.time() - T
        return smooth

    def get_vector(self):
        T = time.time()
        P = len(self.FB)
        if P < 100:
            V = -5
        else:
            V = 0
        smooth = []
        i = 0
        while P*1 > V:
            n = self.run_sieve()
            V += len(n)
            smooth += n
            i += 1
            if i % 50 == 0:
                if i == 50:
                    _log.info("*/ Per 50 times changing poly report /* ")
                _log.info("Time of deciding coefficient = %f sec" % self.coefficienttime)
                _log.info("Sieving Time = %f sec" % self.sievingtime)
                _log.info("Find smooth numbers are %d / %d" % (V, P))
        if P < 100:
            V += 5
        self.smooth = smooth
        _log.info("*/ Total %d times changing poly report /*" % i)
        _log.info("Time of deciding coefficient = %f sec" % self.coefficienttime)
        _log.info("Sieving Time = %f sec" % self.sievingtime)
        _log.info("Found smooth numbers are %d / %d" % (V, P))
        _log.info("Total time of getting enough smooth numbers = %f sec" % (time.time()-T))
        return smooth


class Elimination(object):
    def __init__(self, smooth):
        self.vector = []
        self.history = []
        i = 0
        for vec in smooth:
            self.vector.append(vec[0])
            self.history.append({i:1})
            i += 1
        self.FB_number = len(self.vector[0])
        self.row_size = len(self.vector)
        self.historytime = 0

    def vector_add(self, i, j):
        V_i = self.vector[i]
        V_j = self.vector[j]
        k = 0
        while k < len(V_i):
            if V_i[k] == 1:
                if V_j[k] == 1:
                    V_j[k] = 0
                else:
                    V_j[k] = 1
            k += 1

    def transpose(self):
        Transe_vector = []
        i = 0
        while i < self.FB_number:
            j = 0
            vector = []
            while j < self.row_size:
                vector.append(self.vector[j][i])
                j += 1
            Transe_vector.append(vector)
            i += 1
        self.Transe = Transe_vector

    def history_add(self, i, j):
        T = time.time()
        H_i = self.history[i].keys()
        H_j = self.history[j].keys()
        for k in H_i:
            if k in H_j:
                del self.history[j][k]
            else:
                self.history[j][k] = 1
        self.historytime += time.time() - T

    def gaussian(self):
        T = time.time()
        pivot = []
        FBnum = self.FB_number
        Smooth = len(self.vector)
        for j in range(self.FB_number):
            for k in range(Smooth):
                if k in pivot or not self.vector[k][j]:
                    continue
                pivot.append(k)
                V_k = self.vector[k]
                for h in range(Smooth):
                    if h in pivot or not self.vector[h][j]:
                        continue
                    self.history_add(k, h)
                    V_h = self.vector[h]
                    for q in range(j, FBnum):
                        if V_k[q]:
                            V_h[q] = not V_h[q]
                break
        self.pivot = pivot

        zero_vector = []
        for check in range(Smooth):
            if check not in pivot:
                g = 0
                while g < FBnum:
                    if self.vector[check][g] == 1:
                        break
                    g += 1
                if g == FBnum:
                    zero_vector.append(check)
        _log.info("Time of Gaussian Elimination = %f sec" % (time.time() - T))
        return zero_vector

def qs(n, s, f):
    """
    This is main function of QS
    Arguments are (composite_number, sieve_range, factorbase_size)
    You must input these 3 arguments.
    """
    a = time.time()
    Q = QS(n, s, f)
    _log.info("Sieve range is [ %d , %d ] , Factorbase size = %d , Max Factorbase %d" % (Q.move_range[0], Q.move_range[-1], len(Q.FB), Q.maxFB))
    Q.run_sieve()
    V = Elimination(Q.smooth)
    A = V.gaussian()
    _log.info("Found %d linearly dependent relations" % len(A))
    answerX_Y = []
    N_factors = []
    for i in A:
        B = V.history[i].keys()
        X = 1
        Y = 1
        for j in B:
            X *= Q.smooth[j][1][0]
            Y *= Q.smooth[j][1][1]
            Y = Y % Q.number
        X = sqrt_modn(X, Q.number)
        answerX_Y.append(X-Y)
    for k in answerX_Y:
        if k != 0:
            factor = gcd.gcd(k, Q.number)
            if factor not in N_factors and factor != 1 and \
               factor != Q.number and prime.primeq(factor) == 1:
                N_factors.append(factor)
    N_factors.sort()
    _log.info("Total time = %f sec" % (time.time()-a))
    _log.info(str(N_factors))

def mpqs(n, s=0, f=0, m=0):
    """
    This is main function of MPQS.
    Arguments are (composite_number, sieve_range, factorbase_size, multiplier)
    You must input composite_number at least.
    """
    T = time.time()
    M = MPQS(n, s, f, m)
    _log.info("Sieve range is [ %d , %d ] , Factorbase size = %d , Max Factorbase %d" % (M.move_range[0], M.move_range[-1], len(M.FB), M.maxFB))
    M.get_vector()
    N = M.number // M.multiplier
    V = Elimination(M.smooth)
    A = V.gaussian()
    _log.info("Found %d linerly dependent relations" % len(A))
    answerX_Y = []
    N_prime_factors = []
    N_factors = []
    output = []
    for i in A:
        B = V.history[i].keys()
        X = 1
        Y = 1
        for j in B:
            X *= M.smooth[j][1][0]
            Y *= M.smooth[j][1][1]
            Y = Y % M.number
        X = sqrt_modn(X, M.number)
        if X != Y:
            answerX_Y.append(X-Y)
    NN = 1
    for k in answerX_Y:
        factor = gcd.gcd(k, N)
        if factor not in N_factors and factor != 1 and factor != N \
               and factor not in N_prime_factors:
            if prime.primeq(factor):
                NN = NN*factor
                N_prime_factors.append(factor)
            else:
                N_factors.append(factor)

    _log.info("Total time = %f sec" % (time.time() - T))

    if NN == N:
        _log.debug("Factored completely!")
        N_prime_factors.sort()
        for p in N_prime_factors:
            N = N // p
            i = arith1.vp(N, p, 1)[0]
            output.append((p, i))
        return output
    elif NN != 1:
        f = N // NN
        if prime.primeq(f):
            N_prime_factors.append(f)
            _log.debug("Factored completely !")
            N_prime_factors.sort()
            for p in N_prime_factors:
                N = N // p
                i = arith1.vp(N, p, 1)[0]
                output.append((p, i))
            return output

    for F in N_factors:
        for FF in N_factors:
            if F != FF:
                Q = gcd.gcd(F, FF)
                if prime.primeq(Q) and Q not in N_prime_factors:
                    N_prime_factors.append(Q)
                    NN = NN*Q

    N_prime_factors.sort()
    for P in N_prime_factors:
        i, N = arith1.vp(N, P)
        output.append((P, i))

    if  N == 1:
        _log.debug("Factored completely!! ")
        return output

    for F in N_factors:
        g = gcd.gcd(N, F)
        if prime.primeq(g):
            N_prime_factors.append(g)
            N = N // g
            i = arith1.vp(N, g, 1)[0]
            output.append((g, i))
    if N == 1:
        _log.debug("Factored completely !! ")
        return output
    elif prime.primeq(N):
        output.append((N, 1))
        _log.debug("Factored completely!!! ")
        return output
    else:
        N_factors.sort()
        _log.error("Sorry, not factored completely")
        return output, N_factors



########################################################
#                                                      #
# Following functions are subfunction for main program #
#                                                      #
########################################################

def eratosthenes(n):
    if n < 2:
        raise ValueError("Input value must be bigger than 1")
    else:
        primes = [2]
        sieves = [1] * (((n+1) >> 1) - 1)
        k = 3
        i = 0
        sieves_len = len(sieves)
        while k <= math.sqrt(n):
            if sieves[i] == 1:
                j = 1
                while i+j*k <= sieves_len-1:
                    sieves[i+j*k] = 0
                    j = j+1
                k, i = k+2, i+1
            else:
                k, i = k+2, i+1
        y = 3
        for x in sieves:
            if x == 1:
                primes.append(y)
            y = y+2
        return primes

def prime_mod8(n):
    """
    This is table for choosing multiplier which makes N to have
    factorbase(2,3,5,7,11,13)
    """
    primes = eratosthenes(n)
    PrimeList = {1:[], 3:[], 5:[], 7:[]}
    LegendreList = {1:[], 3:[], 5:[], 7:[]}
    sp = [2, 3, 5, 7, 11, 13]
    for p in primes:
        if p not in sp:
            leg = [arith1.legendre(p, q) for q in sp[1:]]
            if leg not in PrimeList[p & 7]:
                LegendreList[p & 7].append(leg)
                PrimeList[p & 7].append([p, leg])
    return [PrimeList[1], PrimeList[3], PrimeList[5], PrimeList[7]]

def eratosthenes_log(n):
    primes = eratosthenes(n)
    primes_log = []
    for i in primes:
        l = int(math.log(i)*30) #30 is scale
        primes_log.append(l)
    return primes_log

def sqrt_modn(n, modulo):
    import nzmath.factor.methods as methods
    factorOfN = methods.trialDivision(n)
    prod = 1
    for p, e in factorOfN:
        prod = (prod * pow(p, e >> 1, modulo)) % modulo
    return prod

def sqroot_power(a, p, n):
    """
    return squareroot of a mod p^k for k = 2,3,...,n
    """
    x = arith1.modsqrt(a, p)
    answer = [[x, p - x]]
    ppower = p
    inverse = arith1.inverse(x << 1, p)
    for i in range(n - 1):
        x += (a - x ** 2) // ppower * inverse % p * ppower
        ppower *= p
        answer.append([x, ppower - x])
    return answer

#################
# Initial items #
#################
primes_table = eratosthenes(10**5)
primes_log_table = eratosthenes_log(10**5)
prime_8 = prime_mod8(8090)
parameters_for_mpqs = [[100,20], #9 -digits
                      [100,21], #10
                      [100,22], #11
                      [100,24], #12
                      [100,26], #13
                      [100,29], #14
                      [100,32], #15
                      [200,35], #16
                      [300,40], #17
                      [300,60], #18
                      [300,80], #19
                      [300,100], #20
                      [300,100], #21
                      [300,120], #22
                      [300,140], #23
                      [600,160], #24
                      [900,180], #25
                      [1200,200], #26
                      [1000,220], #27
                      [2000,240], #28
                      [2000,260], #29
                      [2000,325], #30
                      [2000,355], #31
                      [2000,375], #32
                      [3000,400], #33
                      [2000,425], #34
                      [2000,550], #35
                      [3000,650], #36
                      [5000,750], #37
                      [4000,850], #38
                      [4000,950], #39
                      [5000,1000], #40
                      [14000,1150], #41
                      [15000,1300], #42
                      [15000,1600], #43
                      [15000,1900], #44
                      [15000,2200], #45
                      [20000,2500], #46
                      [25000,2500], #47
                      [27500,2700], #48
                      [30000,2800], #49
                      [35000,2900], #50
                      [40000,3000], #51
                      [50000,3200], #52
                      [50000,3500]] #53


###
### only find a factor
###

def mpqsfind(n, s=0, f=0, m=0, verbose=False):
    """
    This is main function of MPQS.
    Arguments are (composite_number, sieve_range, factorbase_size, multiplier)
    You must input composite_number at least.
    """
    # verbosity
    if verbose:
        _log.setLevel(logging.DEBUG)
        _log.debug("verbose")
    else:
        _log.setLevel(logging.NOTSET)

    starttime = time.time()
    M = MPQS(n, s, f, m)
    _log.info("Sieve range is [%d, %d]" % (M.move_range[0], M.move_range[-1]))
    _log.info("Factorbase size = %d, Max Factorbase %d" % (len(M.FB), M.maxFB))
    M.get_vector()
    N = M.number // M.multiplier
    V = Elimination(M.smooth)
    A = V.gaussian()
    _log.info("Found %d linearly dependent relations" % len(A))
    differences = []
    for i in A:
        B = V.history[i].keys()
        X = 1
        Y = 1
        for j in B:
            X *= M.smooth[j][1][0]
            Y *= M.smooth[j][1][1]
            Y = Y % M.number
        X = arith1.floorsqrt(X) % M.number
        if X != Y:
            differences.append(X-Y)

    for diff in differences:
        divisor = gcd.gcd(diff, N)
        if 1 < divisor < N:
            _log.info("Total time = %f sec" % (time.time() - starttime))
            return divisor
