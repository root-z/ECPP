'''
author: Z.Lin
'''

import random
from jacobi import jacobi

def atkin_morain():
    pass

'''
p -- a probable prime
'''
def generateCurve(p):
    g= genQNR(p)

'''
generating quardratic residue
p -- prime
'''
def genQNR(p):
    #generate random number g
    g = random.randrange(2, p)
    #1. g has to be quadratic nonresidue
    #2. repeat if p = 1 (mod 3) and g^((p-1)/3) = 1 (mod p)
    while(jacobi(g, p) != -1 or (p % 3 == 1 and pow(g, (p-1)/3, p) == 1)):
        g = random.randrange(2, p)
    return g
    
if __name__=='__main__':
    print(genQNR(17))
