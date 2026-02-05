# coding: utf-8
"""
Implementation of noncommutative ideals in free algebras

AUTHOR:

- Clemens Hofstadler (2023-03-01): initial version

"""

# ****************************************************************************
#                          Copyright (C) 2023
#      Clemens Hofstadler(clemens.hofstadler@jku.at)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import itertools
from copy import deepcopy

from .f4 import F4
from .free_algebra import MyFreeAlgebra
from .nc_polynomial import NCPolynomial
from .nc_ideal_right import NCIdeal_right
from .normal_form import reduced_form, interreduce
from .quiver import Quiver

############################################################################
#  NCIdeal
############################################################################
class NCIdeal:
    def __init__(self,gens,order=None):
        
        F = gens[0].parent()
        X = order if order else F.gens()
        A = MyFreeAlgebra(F.base_ring(),X)        
        
        self.__gens = gens
        self.__internal_gens = [f if isinstance(f,NCPolynomial) else A(f) for f in gens]
        self.__parent = A
        self.__G = []
        self.__algo = None
############################################################################
    def gens(self): return self.__gens
    def parent(self): return self.__parent
    def internal_gens(self): return self.__internal_gens
    def order(self): return self.__parent.order()
    def algo(self): return self.__algo
    def base_ring(self): return self.__parent().base_ring()
############################################################################    
    def __repr__(self):
        return "NCIdeal %s of %s" % (str(tuple(self.__gens)),str(self.__parent))
############################################################################
    def groebner_basis(self,maxiter=10,maxdeg=-1,trace_cofactors=True,criterion=True,reset=True,verbose=0):
        if not self.__algo:
            self.__algo = F4(self)
        self.__G = self.__algo.compute_basis(maxiter=maxiter,maxdeg=maxdeg,trace_cofactors=trace_cofactors,criterion=criterion,reset=reset,verbose=verbose)
        return self.__G
############################################################################
    def right_groebner_basis(self,degbound,relevant_variables=None,quiver=None,maxiter=10,maxdeg=-1,trace_cofactors=True,criterion=True,reset=True,verbose=0):
        G = self.groebner_basis(maxiter=maxiter,maxdeg=maxdeg,trace_cofactors=trace_cofactors,criterion=criterion,reset=reset,verbose=verbose)
        G = interreduce(G)
                
        T = self.__parent.translator()
        if relevant_variables:
            relevant_variables = {str(x) for x in relevant_variables}
        else:
            relevant_variables = set(T.gens())
        
        if quiver:
            X = quiver.to_dict()
            for k,v in X.items(): X[k] = [(t,T(x)) for t,x in v]
        else:
            X = {0 : [(0,T(x)) for x in relevant_variables]}
        
        right_GB = deepcopy(G)
        if quiver:
            # first part is target of g
            new_part = [(quiver.signature(g).pop()[1],g) for g in right_GB]
        else:
            new_part = [(0,g) for g in right_GB]
        min_deg = min(g[1].degree() for g in new_part)
        
        for _ in range(min_deg,degbound):
            new_part = [g for g in new_part if g[1].degree() < degbound]
            # tg = target g, t = target x
            new_part = [(t,g.lrmul(x,'')) for tg,g in new_part for t,x in X[tg] if x in relevant_variables]
            right_GB += [g for _,g in new_part]
                           
        return right_GB

############################################################################
    def __add__(self,other):
        if self.__parent != other.__parent:
            raise ValueError("Ideals have to be defined over the same ring.")
        
        return NCIdeal(self.__gens + other.__gens, order=self.__order)
############################################################################        
    def reduced_form(self,f,maxiter=10,maxdeg=-1,trace_cofactors=True,criterion=True,reset=True,verbose=0):
        G = self.groebner_basis(maxiter=maxiter,maxdeg=maxdeg,trace_cofactors=trace_cofactors,criterion=criterion,reset=reset,verbose=verbose)
        if not isinstance(f,NCPolynomial):
            g = self.parent()(f)
        else:
            g = f
        h = reduced_form(G, g, trace_cofactors=trace_cofactors)
        return h
############################################################################        
    def ideal_membership(self,f,maxiter=10,maxdeg=-1,criterion=True,reset=True,verbose=0):
        h = self.reduced_form(f, maxiter=maxiter, maxdeg=maxdeg, criterion=criterion, reset=reset, verbose=verbose)
        if h.is_zero():
            print("Ideal membership verified!")
            return h.cofactors()
        else:
            print("Ideal membership test inconclusive!")
            print("Reduced form is %s " % str(h))
            return None
############################################################################        
    def find_equivalent_expression(self,exp,order=None,maxiter=10,maxdeg=-1,trace_cofactors=True,criterion=True,reset=True,verbose=0,heuristic='groebner',prefix=None,suffix=None,degbound=5,relevant_variables=None,quiver=None):
            
        if order is None:
            I = self
        else:
            I = NCIdeal(self.gens(),order)   
        
        if not isinstance(exp,NCPolynomial):
            exp = I.parent()(exp)  
                
        if heuristic == 'naive':
            return I._brute_force_heuristic(exp,degbound,relevant_variables,quiver,prefix,suffix,
                                            maxiter=maxiter,maxdeg=maxdeg,trace_cofactors=False,
                                            criterion=criterion,reset=reset,verbose=verbose)
        elif heuristic == 'groebner':
            G = I.groebner_basis(maxiter=maxiter,maxdeg=maxdeg,trace_cofactors=trace_cofactors,
                                criterion=criterion,reset=reset,verbose=verbose)
            G += interreduce(G)
            G = set(G)
        elif heuristic == 'subalgebra':
            raise NotImplementedError
        elif heuristic == 'right-ideal':
            if not prefix: raise ValueError("Prefix required")
            gens = [exp, prefix]
            J = NCIdeal_right(gens, order=I.order())
            G = I.intersect_with_one_sided_ideal(J,degbound=degbound,relevant_variables=relevant_variables,quiver=quiver,
                                                    maxiter=maxiter,maxdeg=maxdeg,trace_cofactors=trace_cofactors,criterion=criterion,
                                                    reset=reset,verbose=verbose)
        elif heuristic == 'left-ideal':
            raise NotImplementedError
        else:
            raise NotImplementedError
        
        exp_mon = set(exp.monomials())
        return [g for g in G if exp_mon.issubset(g.monomials())]

############################################################################
    def intersect(self,other,maxiter=10,maxdeg=-1,trace_cofactors=True,criterion=True,reset=True,verbose=0):
        r"""
        Enumerate a Gröbner basis of the intersection of two two-sided
        ideals of noncommutative polynomials.    

        INPUT: 
    
        - ``other`` -- an NCIdeal
    
        OUTPUT: A (partial) Gröbner basis of the interection `(I) \cap (J)`.
            
        EXAMPLES::
     
            sage: from OperatorGB import *
            sage: F.<x,y> = FreeAlgebra(QQ,2)
            sage: I = NCIdeal([x])
            sage: J = NCIdeal([y])
            sage: I.intersect(J)
            [x*y, y*x]
        """
        if self.__parent != other.__parent:
            raise ValueError("Ideals have to be defined over the same ring.")
        
        A = self.__parent
        # change order
        # maybe change this ?
        new_var = 'tmp_var'
        order = self.order()
        if isinstance(order[0],list):
            order = order + [[new_var]]
        else:
            order = [order,[new_var]]
        
        # set up new ideal
        B = MyFreeAlgebra(A.base_ring(),order)
        t = B.F()(new_var)
        commutator = [B(t*x - x*t) for x in B.F().gens() if x != t]
        one_minus_t = B(1-t)
        t = B(t)
        gens = [t * f.change_parent(B) for f in self.internal_gens()] + [one_minus_t * f.change_parent(B) for f in other.internal_gens()]
        IJ = NCIdeal(gens+commutator,order)
        
        # enumerate Gröbner basis
        G = IJ.groebner_basis(maxiter=maxiter,maxdeg=maxdeg,trace_cofactors=trace_cofactors,criterion=criterion,reset=reset,verbose=verbose)
        
        # select elements free of new_var
        G = [g.change_parent(A) for g in G if new_var not in g.variables()]
    
        return G
############################################################################
    def intersect_with_subalgebra(self,subalgebra,order=None,split_at=0,maxiter=10,maxdeg=-1,trace_cofactors=True,criterion=True,reset=True,verbose=0):
        r"""
        Enumerate a Gröbner basis of the intersection of two two-sided
        ideals of noncommutative polynomials.    

        INPUT: 
    
        - ``subalgebra`` -- a list of noncommutative polynomials
    
        OUTPUT: A (partial) Gröbner basis of the interection `(I) \cap (J)`.
            
        EXAMPLES::
     
            sage: from OperatorGB import *
            sage: F.<x,y> = FreeAlgebra(QQ,2)

        """        
        # change order
        # maybe change this ?
        new_vars = ['tmp_var_' + str(i) for i in range(len(subalgebra))]
        if not order:
            order = self.order()
        if isinstance(order[0],list):
            order = [new_vars] + order
        else:
            order = [new_vars,order]
        
        # set up new ideal
        A = self.__parent
        B = MyFreeAlgebra(A.base_ring(),order)
        T = [B(t) for t in new_vars]
        subalgebra = [B(f) for f in subalgebra]
        gens = [t - f for t,f in zip(T[:split_at],subalgebra[:split_at])]
        gens += [f.change_parent(B) for f in self.__internal_gens]
        gens += [t - f for t,f in zip(T[split_at:],subalgebra[split_at:])]
        IJ = NCIdeal(gens,order)
        
        # enumerate Gröbner basis
        G = IJ.groebner_basis(maxiter=maxiter,maxdeg=maxdeg,trace_cofactors=trace_cofactors,criterion=criterion,reset=reset,verbose=verbose)
        
        # select elements free of old variables
        new_vars_set = set(new_vars)
        G = [g for g in G if set(g.variables()).issubset(new_vars_set)]
        for i in range(len(G)):
            g = G[i]
            for t,f in zip(new_vars,subalgebra):
                g = g.subs(t,f)
            G[i] = g.change_parent(A)
        return G
############################################################################
    def intersect_with_one_sided_ideal(self,other,degbound=5,relevant_variables=None,quiver=None,maxiter=10,maxdeg=-1,trace_cofactors=True,criterion=True,reset=True,verbose=0):
        r"""
        Enumerate a Gröbner basis of the intersection of ...
        """
    
        if self.__parent != other.parent():
            raise ValueError("Ideals have to be defined over the same ring.")
        
        right_GB = self.right_groebner_basis(degbound,relevant_variables=relevant_variables,quiver=quiver, 
                        maxiter=maxiter,maxdeg=maxdeg,trace_cofactors=False,criterion=criterion,reset=reset,verbose=verbose)
        
        # change order
        # maybe change this ?
        new_var = 'tmp_var'
        order = self.order()
        if isinstance(order[0],list):
            order = order + [[new_var]]
        else:
            order = [order,[new_var]]
        
        # set up new ideal
        A = self.__parent
        B = MyFreeAlgebra(A.base_ring(),order)
        t = B.F()(new_var)
        one_minus_t = B(1-t)
        t = B(t)
        gens = [t * f.change_parent(B) for f in right_GB] + [one_minus_t * f.change_parent(B) for f in other.internal_gens()]
        IJ = NCIdeal_right(gens,order)
        
        # enumerate Gröbner basis
        G = IJ.groebner_basis()
        
        # select elements free of new_var
        G = [g.change_parent(A) for g in G if new_var not in g.variables()]
    
        return G       
############################################################################        
    def apply_left_cancellability(self,a,b,maxiter=10,heuristic='subalgebra',order=None,degbound=5,split_at=0,maxdeg=-1,trace_cofactors=True,criterion=True,reset=True,verbose=0,quiver=None):
        
        if order is None:
            I = self
        else:
            I = NCIdeal(self.gens(),order)   
        
        ab = a*b
        if isinstance(ab,NCPolynomial):
            ab = I.parent()(ab)  
        
        if heuristic == 'subalgebra':
            subalgebra = [ab] + self.parent().gens()
            G = I.intersect_with_subalgebra(subalgebra,split_at=split_at,maxiter=maxiter,maxdeg=maxdeg,trace_cofactors=trace_cofactors,criterion=criterion,reset=reset,verbose=verbose)
        elif heuristic == 'two-sided':
            J = NCIdeal([ab], order=I.order())
            G = I.intersect(J,maxiter=maxiter,maxdeg=maxdeg,trace_cofactors=trace_cofactors,criterion=criterion,reset=reset,verbose=verbose)
        elif heuristic == 'one-sided':
            J = NCIdeal_right([ab], order=I.order())
            G = I.intersect_with_one_sided_ideal(J,degbound=degbound,quiver=quiver,maxiter=maxiter,maxdeg=maxdeg,
                                                trace_cofactors=trace_cofactors,criterion=criterion, reset=reset,verbose=verbose)  
        else:
            raise NotImplementedError("Heuristic %s not implemented" % heuristic)
        
        # pick those elements that are in (ab)_\rho
        # and compute x such that g = abx for g in G
        G = [g.is_right_multiple_of(ab) for g in G]
        G = [g for g in G if g not in {True,False}]
        
        # return bx
        A = self.parent()
        if not isinstance(b, NCPolynomial): b = A(b)
        return [b * x for x in G]
############################################################################        
    def apply_right_cancellability(self,a,b,maxiter=10,heuristic='subalgebra',order=None,degbound=5,split_at=0,maxdeg=-1,trace_cofactors=True,criterion=True,reset=True,verbose=0,quiver=None):
        
        if order is None:
            I = self
        else:
            I = NCIdeal(self.gens(),order)   
        
        ab = a*b
        if isinstance(ab,NCPolynomial):
            ab = I.parent()(ab) 
        
        if heuristic == 'subalgebra':
            subalgebra = [ab] + self.parent().gens()
            G = I.intersect_with_subalgebra(subalgebra,split_at=split_at,maxiter=maxiter,maxdeg=maxdeg,trace_cofactors=trace_cofactors,criterion=criterion,reset=reset,verbose=verbose)
        elif heuristic == 'two-sided':
            J = NCIdeal([ab], order=I.order())
            G = I.intersect(J,maxiter=maxiter,maxdeg=maxdeg,trace_cofactors=trace_cofactors,criterion=criterion,reset=reset,verbose=verbose)
        elif heuristic == 'one-sided':
            raise NotImplementedError("Heuristic %s not implemented yet" % heuristic)
        else:
            raise NotImplementedError("Heuristic %s not implemented" % heuristic)
        
        # pick those elements that are in (ab)_\lambda
        # and compute x such that g = xab for g in G
        G = [g.is_left_multiple_of(ab) for g in G]
        G = [g for g in G if g not in {True,False}]
        
        # return xa
        A = self.parent()
        if not isinstance(a, NCPolynomial): a = A(a)
        return [x * a for x in G]

############################################################################        
    def _brute_force_heuristic(self,f,degbound,relevant_variables,Q,prefix,suffix,maxiter=10,maxdeg=-1,trace_cofactors=True,criterion=True,reset=False,verbose=0):  
        
        G = self.groebner_basis(maxiter=maxiter,maxdeg=maxdeg,trace_cofactors=trace_cofactors,criterion=criterion,reset=reset,verbose=verbose)
        A = self.__parent
        T = A.translator()
        
        if relevant_variables:
            relevant_variables = {str(x) for x in relevant_variables}
        else:
            relevant_variables = set(A.gens())
        if not Q:
            Q = Quiver.trivial_quiver(A.F().gens())
        
        one = A.F()(1)
        if prefix is None: prefix = one
        if suffix is None: suffix = one
        prefix_sig = {s for (s,t) in Q.signature(prefix)}
        suffix_sig = {t for (s,t) in Q.signature(suffix)}
        mons = [(v,v,A(one)) for v in Q.G().vertices(sort=False) if v in prefix_sig]
        
        X = Q.to_dict(by_source=False)        
        prefix = A(prefix)
        suffix = A(suffix)
        
        for d in range(degbound+1):
            for s,t,m in mons:
                if s not in suffix_sig: continue
                m = prefix * m * suffix
                g = f-m
                if g.is_zero(): continue
                h = reduced_form(G,g)
                if h.is_zero(): return [g]
            if d <= degbound:
                mons = [(sx, t, m.lrmul('',T(x))) for (s,t,m) in mons for sx,x in X[s] if x in relevant_variables]
                
        return None       