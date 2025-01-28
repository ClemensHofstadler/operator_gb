
from __future__ import absolute_import
from sage.all import *

from .auxiliary import * 
from .nc_monomial import NCMonomial

        
############################################################################
############################################################################
# Cofactor
############################################################################
############################################################################
class Cofactor:
    def __init__(self,c,*args):

        # in this case, we get c a i b
        if args:
            self.__c = c
            self.__a = args[0]
            self.__i = args[1]
            self.__b = args[2]
            self.__parent = args[3]
        # in this case, we get another cofactor and can just copy
        else:        
            self.__c = c.__c
            self.__a = c.__a
            self.__i = c.__i
            self.__b = c.__b
            self.__parent = c.__parent
############################################################################    
    def c(self): return self.__c
    def a(self): return self.__a
    def i(self): return self.__i
    def b(self): return self.__b
    def aib(self): return self.__a, self.__i, self.__b
    def caib(self): return self.__c, self.__a, self.__i, self.__b
    def parent(self): return self.__parent
############################################################################     
    def __repr__(self):
        c = ''
        cc = self.__c
        if cc < 0: 
            c += '-'
            cc = -cc
        if cc != 1: c += str(cc) + "*"
        a = c + my_pretty_print(self.__a, self.__parent)
        i = str(self.__i)
        b = my_pretty_print(self.__b, self.__parent)
        
        return '(' + ','.join([a,i,b]) + ')'
############################################################################
    def pretty_print(self,f):
        c = ''
        cc = self.__c
        if cc > 0:
            c += ' + '
        elif cc < 0: 
            c += ' - '
            cc = -cc
        if cc != 1: c += str(cc) + "*"
        if self.__a:
            a = c + my_pretty_print(self.__a, self.__parent) + "*"
        else:
            a = c
        if self.__b:
            b = "*" + my_pretty_print(self.__b, self.__parent)
        else:
            b = ''
        
        return a + "(" + str(f) + ")" + b  
############################################################################
    def __copy__(self): return Cofactor(self)
############################################################################    
    def __eq__(self,other):
         return self.__c == other.__c and self.__i == other.__i and \
                self.__a == other.__a and self.__b == other.__b
    def __ne__(self,other):
        return not (self == other)
############################################################################
    def __hash__(self):
         return hash(tuple(self.__c, self.__a, self.__i, self.__b))
############################################################################
    def __mul__(self, other):
        out = copy(self)
        out.__c *= other
        return out
############################################################################ 
    def __truediv__(self,other):
        out = copy(self)
        out.__c /= other
        return out
############################################################################
    def multiply_by(self,c):
        self.__c *= c
############################################################################
    def lmul(self, l):
        out = copy(self)
        out.__a = l + out.__a
        return out
############################################################################
    def rmul(self, r):
        out = copy(self)
        out.__b += r
        return out
############################################################################
    def lrmul(self, l, r):
        out = copy(self)
        out.__a = l + out.__a
        out.__b += r
        return out
############################################################################
    def to_free_algebra(self, F, f):
        c = self.__c
        a = F(my_pretty_print(self.__a, self.__parent))
        b = F(my_pretty_print(self.__b, self.__parent))
        return (c*a, f, b)