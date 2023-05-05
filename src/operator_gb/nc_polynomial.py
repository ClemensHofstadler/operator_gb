# coding: utf-8
"""
NCPoly
================

A class that provides basic funcionality for noncommutative polynomials

AUTHOR:

- Clemens Hofstadler (2023-02-09)

"""

from sage.all import prod
from copy import copy,deepcopy

import operator_gb.free_algebra
from .nc_monomial import NCMonomial
from .cofactor import Cofactor

############################################################################
############################################################################
# NCPolynomial
############################################################################
############################################################################
class NCPolynomial:

    def __init__(self,coeffs,mons=None):
        if mons is not None:
            self.__lm = mons[-1]
            self.__mons = mons
            self.__coeffs = coeffs  
            self.__cofactors = []
        else:
            P = coeffs.parent()
            A = operator_gb.free_algebra.MyFreeAlgebra(P.base_ring(), P.gens())
            f = A(coeffs)
            self.__lm = f.__lm
            self.__mons = f.__mons
            self.__coeffs = f.__coeffs
            self.__cofactors = []
###########################################################################
    def parent(self): return self.__lm.parent()
    def lm(self): return self.__lm
    def lc(self): return self.__coeffs[-1]
    def monomials(self): return self.__mons
    def coefficients(self): return self.__coeffs
    def cofactors(self): return self.__cofactors
############################################################################
    def variables(self,internal=False):
        T = self.parent().translator()
        out = set()
        for m in self.__mons:
            out.update(set(str(m)))
        return [T(x,to_internal=internal) for x in out]
############################################################################
    def coefficients_monomials(self): return zip(self.__coeffs, self.__mons)
############################################################################
    def __copy__(self): return NCPolynomial(copy(self.__coeffs),deepcopy(self.__mons))
############################################################################
    def to_native(self): return self.parent()(self)
############################################################################
    def __repr__(self):
        s = ""
        for c,m in zip(self.__coeffs,self.__mons):
            # sign
            if c == 0: continue
            if c > 0: s += " + "
            else: s += " - "
            
            if abs(c) != 1: 
                s += str(abs(c)) + "*"        
            s += m.pretty_print()
        
        # remove first + if existent
        if not s: s = "0"
        elif s[1] == "+": s = s[3:]
        else: s = s[1:] 
        
        return s
###########################################################################
    def __eq__(self, other):
        if len(self.__coeffs) != len(other.__coeffs): return False
        return self.__lm == other.__lm and self.__coeffs == other.__coeffs and self.__mons == other.__mons
    def __ne__(self,other):
        return not (self == other)
############################################################################
    def  __hash__(self):
        return hash((tuple(self.__mons),tuple(self.__coeffs)))
############################################################################
    def __neg__(self):
        out = self.__copy__()
        out.__coeffs = [-c for c in self.__coeffs]
        return out   
############################################################################
    def __add__(self,other):
        if self.parent() != other.parent():
            raise ValueError("Elements have to be defined over the same ring.")
        
        d = {m:c for c,m in self.coefficients_monomials()}
        for c,m in other.coefficients_monomials():
            if m in d: d[m] += c
            else: d[m] = c
        
        coeffs,mons = dict_to_poly(d)
        if not mons: return NCPolynomial.zero(self.parent())
        return NCPolynomial(coeffs,mons)
############################################################################
    def  __sub__(self,other): return self.__add__(other.__neg__())
############################################################################    
    def __mul__(self,other):
        if self.parent() != other.parent():
            raise ValueError("Elements have to be defined over the same ring.")
        
        d = dict()
        for c1,m1 in self.coefficients_monomials():
            for c2,m2 in other.coefficients_monomials():
                m1m2 = m1*m2
                if m1m2 in d: d[m1m2] += c1*c2
                else: d[m1m2] = c1*c2
        
        coeffs,mons = dict_to_poly(d)
        if not mons: return NCPolynomial.zero(self.parent())
        return NCPolynomial(coeffs,mons)
###########################################################################
    def coeff_mul(self,c):
        out = self.__copy__()
        out.__coeffs = [coeff * c for coeff in self.__coeffs]
        return out
###########################################################################
    def degree(self): return len(self.__lm)
###########################################################################
    @staticmethod
    def zero(P): return NCPolynomial([P.zero()],[NCMonomial('',P)])
############################################################################                
    def subs(self,x,f):
        T = self.parent().translator()
        to_replace = T(x)
        out = self.__copy__()
        i = 0
        while i < len(out.__mons):
            m = out.__mons[i]
            if to_replace not in str(m):
                i += 1
                continue
            else:
                out.__mons.remove(m)
                coeff = out.__coeffs.pop(i)
                a,b = str(m).split(to_replace,1)
                out.__mons += [m.lrmul(a,b) for m in f.__mons]
                out.__coeffs += [coeff * c for c in f.__coeffs]
        
        # clean up
        d = dict()
        for c,m in out.coefficients_monomials():
            if m in d: d[m] += c
            else: d[m] = c
        
        coeffs,mons = dict_to_poly(d)
        if not mons: return NCPolynomial.zero(self.parent())
        return NCPolynomial(coeffs,mons)        
############################################################################                
    def set_to_zero(self):
        P = self.parent()
        self.__coeffs = [P.zero()]
        self.__mons = [NCMonomial('',P)]
        self.__lm = self.__mons[-1]
############################################################################
    def is_zero(self): 
        z = self.__lm.parent().zero()
        return self.__coeffs == [z]*len(self.__coeffs)
############################################################################
    def change_parent(self,A):
        if self.is_zero():
            return NCPolynomial.zero(A)
        return A(self.__repr__())
############################################################################    
    def make_monic(self):
        """
        Really update self
        """
        lc = self.__coeffs[-1]
        if lc != 1: self.__coeffs = [c / lc for c in self.__coeffs]
        return lc
############################################################################
    def append_cofactor(self, c):
        self.__cofactors.append(c)
############################################################################
    def reset_cofactors(self):
        self.__cofactors = []
############################################################################
    def append_cofactors(self, c):
        self.__cofactors += c
############################################################################
    def lrmul(self, l, r):
        mons = [m.lrmul(l,r) for m in self.__mons]
        coeffs = copy(self.__coeffs)
        return NCPolynomial(coeffs,mons)
############################################################################   
    def is_left_multiple_of(self, g):
        return self.is_one_sided_multiple_of(g,'left')
############################################################################   
    def is_right_multiple_of(self, g):
        return self.is_one_sided_multiple_of(g,'right')
############################################################################   
    def is_one_sided_multiple_of(self, g, side):
        f = self
        # special case: at least one argument is zero
        if f.is_zero(): return True
        if g.is_zero(): return False
        
        assert side == 'left' or side == 'right'
                
        A = self.parent()   
        if not isinstance(g,NCPolynomial): g = A(g)
        lm_g = str(g.lm())
        coeffs,mons = [],[]
        
        while not f.is_zero():
            lm_f = str(f.lm())
            
            if side == 'left':
                if not lm_f.endswith(lm_g): return False
                l = lm_f[:-len(lm_g)]
                r = ''
                mons.append(NCMonomial(l,A))
            else:
                if not lm_f.startswith(lm_g): return False
                l = ''
                r = lm_f[len(lm_g):]
                mons.append(NCMonomial(r,A))
            
            c = f.lc() / g.lc()
            coeffs.append(c)
            f = f - g.lrmul(l,r).coeff_mul(c)
        
        quo = NCPolynomial(coeffs, mons)
        return quo    
############################################################################
    def normal_cofactors(self,gens):
        F = self.__lm.parent().F()
        cofactors = [c.to_free_algebra(F,gens[c.i()]) for c in self.__cofactors]   
        return cofactors 
############################################################################
    def expand_cofactors(self,gens):
        cofactors = self.normal_cofactors(gens)
        return sum(map(prod,cofactors))
############################################################################                
    def simplify_cofactors(self):
        d = dict()
        A = self.parent()
        for cofactor in self.__cofactors:
            c,a,i,b = cofactor.caib()
            if (a,i,b) in d: d[(a,i,b)] += c
            else: d[(a,i,b)] =  c
        
        self.__cofactors = [Cofactor(c,a,i,b,A) for (a,i,b),c in d.items() if c]

def dict_to_poly(d):
    d = list(d.items())
    d.sort(key=lambda p : p[0])
    coeffs = [c for m,c in d if c != 0]
    mons   = [m for m,c in d if c != 0]
    return coeffs, mons