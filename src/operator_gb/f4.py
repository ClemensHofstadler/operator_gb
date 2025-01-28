# coding: utf-8
"""
F4
================

Module that implements the F4 algorithm in the free algebra that
allows to trace cofactors.

AUTHOR:

- Clemens Hofstadler (2023-02-09)

"""

#############################################################################
#  Copyright (C) 2020 Clemens Hofstadler (clemens.hofstadler@jku.at).       #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 2, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from time import time
from copy import copy

import ahocorasick

from .ambiguity import *
from .cofactor import Cofactor
from .free_algebra import MyFreeAlgebra
from .nc_monomial import NCMonomial
from .nc_polynomial import NCPolynomial
from .normal_form import interreduce
from .crit_pair import CritPair
from .linear_algebra import faugere_lachartre

class F4():

    def __init__(self,I):
        self.__parent = I.parent()
        self.__gens = I.gens()
        self.__internal_gens = I.internal_gens()
        self.__criterion = True
        self.__maxdeg = -1
        self.__lm = ahocorasick.Automaton()
        self.__suffix_trie = ahocorasick.Automaton()
        self.__verbose = 0
        self.__amb = []
        self.__G = []
        self.__constant_flag = False
############################################################################        
    def lm(self): return list(self.__lm.items()) 
    def parent(self): return self.__parent
    def gens(self): return self.__gens
    def internal_gens(self): return self.__internal_gens
    def amb(self): return self.__amb
    def G(self): return self.__G
############################################################################
    def clear(self):
        self.__lm.clear()
        self.__suffix_trie.clear()
        self.__amb = []
        self.__G = []
        self.__constant_flag = False
############################################################################
# Compute Gröbner basis
############################################################################
    def compute_basis(self,maxiter=10,maxdeg=-1,trace_cofactors=True,criterion=True,reset=True,verbose=0):
        global nr_pairs
        global zero_reductions
        
        self.__maxdeg = maxdeg
        self.__criterion = criterion
        self.__verbose = verbose
        
        nr_pairs = 0
        zero_reductions = 0
        
        if reset or not self.__G: 
            self.prepare_input(verbose,trace_cofactors=trace_cofactors)
            # constant_flag checks if GB contains 1
            if self.__constant_flag: return self.__G
            self.compute_ambiguities()
        
        self.__already_rewritten = len(self.__G)
        oldlen = len(self.__G)      
        count = 0
                       
        while count < maxiter and self.__amb:
            count += 1
            
            start = time() 
            PP = self.reduction(trace_cofactors=trace_cofactors)
            if verbose > 0:
                print("Reduction took: %.5f" % (time()-start))
                            
            self.add_polynomials(PP)  
            if self.__constant_flag: break                                             
            self.compute_ambiguities(oldlen)
            oldlen = len(self.__G)
            if verbose > 0:
                print("Iteration %d finished. G has now %d elements.\n" % (count,len(self.__G)))
                    
        if (not self.__amb or self.__constant_flag) and verbose > 0:
            print("All critial pairs could be reduced to zero.")
            
        if verbose > 1:
            print("In total, %d critical pairs were reduced." % nr_pairs)
            print("In total, %d reductions to zero occurred." % zero_reductions)
        
        if trace_cofactors:
            self.rewrite_cofactors()
            
        return self.__G

############################################################################
    def compute_ambiguities(self,oldlen=0):
        
        G = self.__G
        maxdeg = self.__maxdeg
        verbose = self.__verbose
        criterion = self.__criterion
                
        words = [str(g.lm()) for g in G]         
            
        start = time()
        prefix_trie = self.__lm
        suffix_trie = self.__suffix_trie
        for i in range(oldlen,len(words)):
            amb_i = Ambiguity.generate_with_tries(prefix_trie,suffix_trie,words,i,criterion=criterion)
            if maxdeg > 0: amb_i = [a for a in amb_i if a.degree() <= maxdeg]
            if criterion:
                self.__amb = Ambiguity.chain_criterion(self.__amb,amb_i,words[i],i)
            self.__amb += amb_i
                                    
        if verbose > 0:
            print("%d ambiguities in total (computation took %.5f)" % (len(self.__amb), time()-start))
                             
        self.__amb.sort(key = lambda p : p.degree())
############################################################################
# Reduction & Symbolic Preprocessing
############################################################################
    def reduction(self,trace_cofactors=True):
        global nr_pairs
        global zero_reductions
                
        G = self.__G
        amb = self.__amb
        verbose = self.__verbose
        
        # choose ambiguities            
        d = amb[0].degree()
        P = [a for a in amb if a.degree() == d]
        self.__amb = amb[len(P):]
        P = [CritPair(a,G[a.i()],G[a.j()]) for a in P]        
        # remove zero S-polynomials
        P = [pair for pair in P if pair.degree() > -1]
        if not P: return []
        
        if verbose > 0:
            print("%d critical pairs will be reduced." % len(P))
        
        nr_pairs += len(P) 
        
        # do symbolic preprocessing
        F = flatten([pair.fg() for pair in P])
        pivot_rows,pivot_columns,columns = self.symbolic_preprocessing(F)

        #split critical pairs into pivot and non-pivot rows
        rest_rows = [pair.g() for pair in P]
        for pair in P:
            f = pair.f()
            if f.lm() in pivot_columns:
                rest_rows.append(f)
            else:
                pivot_columns.add(f.lm())
                pivot_rows.append(f)
        
        
        #seperate pivot and non-pivot columns
        rest_columns = list(columns - pivot_columns)
        pivot_columns = list(pivot_columns)
        pivot_columns.sort(reverse=True)
        rest_columns.sort(reverse=True)
        
        size_A = len(pivot_rows)

        columns = {m:i for (i,m) in enumerate(pivot_columns + rest_columns)}
        pivot_rows.sort(key=lambda f: f.lm(),reverse=True)
        rows = pivot_rows + rest_rows
        
        if verbose > 0:
            print("Set up matrix of size (%d, %d)" % (len(rows), len(columns)))
        
        # M...the rref we are interested in
        # T,T1...the transformation matrix
        T,T1,M = faugere_lachartre(rows,columns,size_A,trace_cofactors=trace_cofactors)
                 
        PP = self.matrix_to_polies(T,T1,M,rest_columns,rows,size_A)
                
        zero_reductions += len(P) - len(PP)
            
        return PP
############################################################################
    def symbolic_preprocessing(self,F):
        done = {f.lm() for f in F}
        todo = {m for f in F for m in f.monomials()} - done
        reducers = []
        pivot_columns = set()
        while todo:
            m = todo.pop()
            done.add(m)
            agb = self.find_reducer(str(m))
            if agb:
                assert m == agb.lm()
                reducers.append(agb)
                pivot_columns.add(agb.lm())
                todo.update(set(agb.monomials()).difference(done))
        return reducers,pivot_columns,done
############################################################################
    def find_reducer(self,m):
        lm = self.__lm
        G = self.__G
                
        reducer = list(lm.iter(m))
        #reducer = next(lm.iter(m),None)
        if not reducer: return None
        #k, (i,lm) = reducer 
        k,(i,lm) = min(reducer, key=lambda r:r[1][0])
        g = G[i]
        a = m[:k-len(lm)+1]
        b = m[k+1:]
        agb = g.lrmul(a,b)
        agb.append_cofactor(Cofactor(1,a,i,b,self.__parent))
        return agb
############################################################################        
    def matrix_to_polies(self,T,T1,M,rest_columns,rows,n):
    
        #compute polynomials
        pos = M.nonzero_positions()
        if not pos: return []
        rank = pos[-1][0]+1
        coefficients = [[] for i in range(rank)]
        monomials = [[] for i in range(rank)]
        for (i,j),c in reversed(list(M.dict().items())):
            coefficients[i].append(c)
            monomials[i].append(rest_columns[j])
    
        PP = [NCPolynomial(c,m) for c,m in zip(coefficients,monomials)]
            
        # take care of cofactors
        for (i,j),tij in T1.dict().items():
            if i >= rank: continue
            PP[i].append_cofactors([c * tij for c in rows[j].cofactors()])

        for (i,j),tij in T.dict().items():
            if i >= rank: continue
            PP[i].append_cofactors([c * tij for c in rows[j+n].cofactors()])

        PP.sort(key=lambda f : f.lm())
        return PP
############################################################################        
# Auxiliary
############################################################################      
    def prepare_input(self,verbose,trace_cofactors=True):
        
        self.clear()
        # change data structure
        A = self.__parent
        G = [copy(f) for f in self.__internal_gens]
        
        # make monic and add cofactors
        for i,g in enumerate(G):
            c = g.make_monic()
            g.reset_cofactors()
            if trace_cofactors: g.append_cofactor(Cofactor(1/c,'',i,'',A))
        
        # interreduce generators
        G = interreduce(G)
        if verbose > 0:
            print("Interreduced the generators from %d elements to %d elements.\n" % (len(self.__gens),len(G)))
        
        # add interreduced generators to GB    
        self.add_polynomials(G)
        
############################################################################                    
    def add_polynomials(self,PP):
        oldlen = len(self.__G)
        for i,p in enumerate(PP):
            m = str(p.lm())
            if m:
                self.__lm.add_word(m,(i+oldlen,m))
                self.__suffix_trie.add_word(m[::-1],(oldlen+i,m))   
            else:
                self.__constant_flag = True
                break

        self.__lm.make_automaton()
        self.__G += PP
############################################################################                
    def rewrite_cofactors(self):
        
        A = self.__parent
        G = self.__G
                
        for g in G[self.__already_rewritten:]:
            cofactors_tmp = []
            for cofactor in g.cofactors():
                c,a,i,b = cofactor.caib()
                cofactors_tmp += [t.lrmul(a,b) * c for t in G[i].cofactors()]
            g.reset_cofactors()
            g.append_cofactors(cofactors_tmp)
            g.simplify_cofactors()








