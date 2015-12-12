'''
author: Z.Lin
'''

import random
from jacobi import jacobi

def atkin_morain():
    pass


def generate_curve(p):
    #calculate quadratic nonresidue
    g = gen_QNR(p)
    #find discriminant
    
    
def find_discriminant(p):
    '''
    really should consider generating list of discriminants
    '''
    # find a fundamental discriminant
    d = 0
    # compute the odd part
    dFound = False
    while (not dFound):
        '''find the next fundamental discriminant'''
        d -= 1
        odd = odd_part(-d)
        #check if the odd part is square free.
        if (not square_free_odd(odd)):
            continue
        if not ((-d) % 16 in {3, 4, 7, 8, 11, 15}):
            continue
        #Have the next discriminant. See if it is good
        if (jacobi(d % p , p) != 1):
            continue
        #connarchia

def odd_part(n):
    #compute the odd part of a number
    oddPart = n
    while (oddPart % 2 == 0):
        oddPart = oddPart//2
    return oddPart
 
def square_free_odd(n):
    #check if the odd part is square free. No efficient algorithm is known.
    x = 1
    xSquare = x * x
    while (xSquare <= n):
        if xSquare == n:
            return False
        else:
            x += 2
            xSquare = x * x
    return True
    
 
'''
generating quardratic residue
p -- prime
'''
def gen_QNR(p):
    #generate random number g
    g = random.randrange(2, p)
    #1. g has to be quadratic nonresidue
    #2. repeat if p = 1 (mod 3) and g^((p-1)/3) = 1 (mod p)
    while(jacobi(g, p) != -1 or (p % 3 == 1 and pow(g, (p-1)/3, p) == 1)):
        g = random.randrange(2, p)
    return g
    
if __name__=='__main__':
    print(gen_QNR(17))
