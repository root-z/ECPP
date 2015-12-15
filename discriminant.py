from nzmath.arith1 import issquare
from hilbert import hilbert

def gen_discriminant(start=0):
    # really should consider generating list of discriminants
    # find a fundamental discriminant. Use start as a starting point.
    d = start
    # compute the odd part
    d_found = False
    while not d_found:
        '''find the next fundamental discriminant'''
        d -= 1
        odd = odd_part(-d)
        # check if the odd part is square free.
        if issquare(odd) not in {0, 1} :
            continue
        if not ((-d) % 16 in {3, 4, 7, 8, 11, 15}):
            continue
        return d
        # Have the next discriminant. See if it is good


def odd_part(n):
    # compute the odd part of a number
    odd = n
    while odd % 2 == 0:
        # odd = odd // 2
        odd //= 2
    return odd


if __name__=='__main__':
    print odd_part(8)
    print issquare(1)
    d = 0
    while -d < 100:
        d = gen_discriminant(d)
        print d