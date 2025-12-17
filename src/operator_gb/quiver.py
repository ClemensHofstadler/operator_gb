# coding: utf-8
"""
Quiver
================

Class that allows to check (unifom) compatibility of polynomials with quivers.

AUTHOR:

- Clemens Hofstadler (2020-02-16)

"""

#############################################################################
#  Copyright (C) 2020 Clemens Hofstadler (clemens.hofstadler@jku.at).       #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 2, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

import itertools
from sage.all import DiGraph

from .free_algebra import MyFreeAlgebra
from .nc_polynomial import NCPolynomial

############################################################################
#  Quiver
############################################################################
class Quiver:
    def __init__(self,triples):
        G = DiGraph(multiedges=True,loops=True)
        self.__vars = []
        for (s,t,l) in triples:
            G.add_edge(str(s),str(t),str(l))
            if str(l) not in self.__vars:
                self.__vars.append(str(l))
        self.__G = G
############################################################################
    def __repr__(self):
        var = ""
        for v in self.vars():
            var += v + ", "
        var = var[:-2]
        
        return "Labelled quiver with %s vertices in the labels {%s}" % (str(self.__G.num_verts()),var)
############################################################################
    def plot(self):
        G = self.__G.plot(edge_labels=True,vertex_labels=False)
        G.show()
############################################################################
    def G(self): return self.__G
    def vars(self): return self.__vars
############################################################################
    def source(self,label):
        return [s for (s,t,l) in self.__G.edge_iterator() if l == label]
############################################################################
    def target(self,label):
        return [t for (s,t,l) in self.__G.edge_iterator() if l == label]
############################################################################
    def add_edge(self, s, t, l):
        self.__G.add_edge(str(s),str(t),str(l))
        if str(l) not in self.__vars:
            self.__vars.append(str(l))
############################################################################    
    def to_dict(self,by_source=True):
        D = dict()
        for (s,t,x) in self.__G.edge_iterator():
            a,b = (s,t) if by_source else (t,s)
            if a in D: D[a].append((b,x))
            else: D[a] = [(b,x)]
        return D
############################################################################
    def __signature_mon__(self,m):
        """
        We can assume that m is of type NCMonomial.
        """
        T = m.parent().translator()
        m = str(m)[::-1]
        #if m = 1, return {(v,v) | v\in V}
        if m == '': return {(v,v) for v in self.__G.vertex_iterator()}

        #usual case normal monomial
        begin = {(s,t) for s,t,l in self.__G.edge_iterator() if T(l) == m[0]}
        for i in range(1,len(m)):
            end = {(s,t) for s,t,l in self.__G.edge_iterator() if T(l) == m[i]}
            comb = {(s1,t2) for (s1,t1),(s2,t2) in itertools.product(begin,end) if t1 == s2}
            if len(comb) == 0:
                return set()
            begin = comb
        return begin
############################################################################
    def signature(self,f):
        if not isinstance(f,NCPolynomial):
            F = f.parent()
            A = MyFreeAlgebra(F.base_ring(),F.gens())
            f = A(f)
              
        if f.is_zero():
            return {pair for pair in itertools.product(self.__G.vertices(sort=False),repeat=2)}
        
        signatures = [self.__signature_mon__(m) for m in f.monomials()]
        return set.intersection(*signatures)    
############################################################################
    def is_compatible(self,f):
        return len(self.signature(f)) > 0
############################################################################
    def is_uniformly_compatible(self,f):
        if not isinstance(f,NCPolynomial):
            F = f.parent()
            A = MyFreeAlgebra(F.base_ring(),F.gens())
            f = A(f)
    
        if not self.is_compatible(f): return False
        
        if f.is_zero(): return True

        monomials = f.monomials()
        mon_signatures = [self.__signature_mon__(m) for m in monomials]
        return all(sig == mon_signatures[0] for sig in mon_signatures)        
############################################################################
    @staticmethod
    def trivial_quiver(X):
        return Quiver([('','',x) for x in X])