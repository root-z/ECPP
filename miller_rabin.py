"""
Implementation of Miller-Rabin test.
"""
import random


def miller_rabin(n, k):
    """
    Main method
    Args:
        n:
        k: number of loops (randomly choosing witnesses)

    Returns:

    """
    for i in range(0, k):
        a = random.randrange(2, n-1)
        if not miller_rabin_once(n, a):
            return False
    return True


def miller_rabin_once(n, a):
    """
    Single Witness loop.
    Args:
        n:
        a:

    Returns:

    """
    if n % 2 == 0 or 1 < gcd(a, n)[0] < n:
        return False
    q, k = twokq(n-1)

    a = pow(a, q, n)
    if a == 1:
        return True
    for i in range(0, k):
        if a == -1 % n:
            return True
        a = pow(a, 2, n)
    return False


def twokq(n):
    k = 0
    while n%2 == 0:
        n //= 2
        k += 1
    return (n, k)


def explain(x):
    if not x:
        print('composite')
    else:
        print('prime')


def gcd(a, b):
    xminus2 = 1
    yminus2 = 0
    xminus1 = 0
    yminus1 = 1
    q = 0
    rminus2 = a
    rminus1 = b

    while rminus1 != 0:
        q = rminus2 // rminus1
        newx = xminus2 - q * xminus1
        newy = yminus2 - q * yminus1
        newr = rminus2 - q * rminus1
        #print(newx, newy, newr, q)

        xminus2 = xminus1
        xminus1 = newx
        yminus2 = yminus1
        yminus1 = newy
        rminus2 = rminus1
        rminus1 = newr

    return (rminus2, xminus2, yminus2)

if __name__=='__main__':
    print(miller_rabin(838041647, 20))
