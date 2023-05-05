
from __future__ import absolute_import
from sage.all import *

from .auxiliary import flatten, simplify_str
from .orderings import deglex,multilex
from .nc_monomial import NCMonomial
import operator_gb.nc_polynomial
from .translator import Translator

    
############################################################################
############################################################################
# FreeAlgebra
############################################################################
############################################################################
class MyFreeAlgebra(Parent):
    def __init__(self,K,X):
                
        if K not in {QQ}:
            raise TypeError("coefficient field %s not supported" % str(K))
        if not len(X): raise ValueError("Need at least one variable")
        
        T = Translator(X)
        self.__translator = T       
        X = self.set_order(X) 
        self.__F = FreeAlgebra(K,X)
        self.__gens_dict = {str(x):i for i,x in enumerate(T.internal_gens())}
        
############################################################################    
    def base_ring(self): return self.__F.base_ring()
    def F(self): return self.__F
    def gens(self): return [str(x) for x in self.__F.gens()]
    def gens_dict(self): return self.__gens_dict
    def order(self): return self.__order
    def name_order(self): return self.__name_order
    def blocks(self): return self.__blocks
    def translator(self): return self.__translator
############################################################################      
    def __eq__(self,other):
        return self.__F == other.__F and self.__name_order == other.__name_order and self.__blocks == other.__blocks
    def __ne__(self,other): return not self.__eq__(other)
############################################################################      
    def __repr__(self):
        s = str(self.__F) + " with "
        if self.__name_order == 'deglex':
            for x in self.gens():
                s += x + " < "
        else: 
            idx = -1
            T = self.__translator
            block = {T(x,to_internal=False) for x in self.__blocks[idx]}
            for x in self.gens():
                if x in block:
                    s += x + " < "
                else:
                    s = s[:-3] + " << " + x + " < "
                    idx -= 1
                    block = {T(x,to_internal=False) for x in self.__blocks[idx]}
        return s[:-3]
############################################################################
    def __call__(self,f):
        F = self.__F
        # from NCPolynomial to FreeAlgebra
        if isinstance(f,operator_gb.nc_polynomial.NCPolynomial):
            s = ""
            for c,m in zip(f.coefficients(),f.monomials()):
                mon = m.pretty_print()
                s += str(c) + "*" + mon + "+"
            return F(s + "0")
        
        # from string to NCPolynomial
        elif isinstance(f,str):
            return self(F(f))
        # from FreeAlgebra to NCPolynomial
        else:
            #if f.parent() != self.__F:
            #    raise ValueError("Parent algebras do not match")
            T = self.__translator
            d = f.monomial_coefficients()
            d = [(d[key],NCMonomial(simplify_str(T(str(key))),self)) for key in d]
            d.sort(key=lambda p : p[1])
            mons = [m for c,m in d]
            coeffs = [c for c,m in d]
            if not mons:
                operator_gb.nc_polynomial.NCPolynomial.zero(self)
            return operator_gb.nc_polynomial.NCPolynomial(coeffs,mons)
############################################################################   
    def set_order(self,X):
        
        self.__order = list(X)
    
        if isinstance(X[0],list):
            T = self.__translator
            self.__name_order = 'multilex'
            self.__cmp = multilex
            self.__blocks = list(reversed([{T(str(x)) for x in block} for block in X if block]))
            return flatten(X)
        else:
            self.__name_order = 'deglex' 
            self.__cmp = deglex
            self.__blocks = None
            return X
        
############################################################################
    def cmp(self,m1,m2): return self.__cmp(m1,m2,self.__gens_dict,self.__blocks)
 ############################################################################   
    def zero(self):
        return self.base_ring().zero() 
