'''
from libnum
'''

def jacobi(a, n):
    """
    Return Jacobi symbol (or Legendre symbol if n is prime)
    """
    s = 1
    while True:
        if n < 1: raise ValueError("Too small module for Jacobi symbol: " + str(n))
        if n & 1 == 0: raise ValueError("Jacobi is defined only for odd modules")
        if n == 1: return s
        a = a % n
        if a == 0: return 0
        if a == 1: return s

        if a & 1 == 0:
            if n % 8 in (3, 5):
                s = -s
            a >>= 1
            continue

        if a % 4 == 3 and n % 4 == 3:
            s = -s

        a, n = n, a
    return
