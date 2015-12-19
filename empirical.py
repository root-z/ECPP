from ecpp import atkin_morain
from nzmath.ecpp import ecpp
from miller_rabin import miller_rabin

if __name__=='__main__':
    n = 2**100 + 1
    print miller_rabin(n, 10)
    print ecpp(n)
    print atkin_morain(n)
