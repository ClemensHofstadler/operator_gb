from __future__ import absolute_import

from .cofactor import Cofactor

from copy import copy

############################################################################
# S-polynomials
############################################################################
class CritPair:
    def __init__(self, amb , gi, gj):
    
        A = gi.parent()
        
        self.__degree = amb.degree()
        Ai,Ci,Aj,Cj = amb.AC()
        g1 = gi.lrmul(Ai,Ci)
        g2 = gj.lrmul(Aj,Cj)

        if g1 == g2:
            self.__degree = -1
        else:
            i,j = amb.ij()
            g1.append_cofactor(Cofactor(1,Ai,i,Ci,A))
            g2.append_cofactor(Cofactor(1,Aj,j,Cj,A))
            
            self.__f = g1
            self.__g = g2
            
############################################################################
    def degree(self): return self.__degree
    def f(self): return self.__f
    def g(self): return self.__g
    def fg(self): return [self.__f, self.__g]
############################################################################
    def __repr__(self):
        return "(" + str(self.__f) + ", " + str(self.__g) + ")"
