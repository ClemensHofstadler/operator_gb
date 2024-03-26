# coding: utf-8
r"""
Ambiguities as defined for noncommutative Gröbner basis computations

Ambiguities characterise situations where one term can be reduced in
two different ways. The goal of noncommutative Gröbner basis computations
is to resolve all ambiguities by reducing S-polynomials formed from them. 

AUTHORS:

- Clemens Hofstadler (2023-03-01): initial version

"""

#############################################################################
#  Copyright (C) 2020 Clemens Hofstadler (clemens.hofstadler@jku.at).       #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 2, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from .auxiliary import flatten
import os

############################################################################
# ambiguities
############################################################################
class Ambiguity:
    def __init__(self,ABC,Ai,Ci,Aj,Cj,i,j):
        self.__ABC = ABC
        self.__Ai = Ai
        self.__Ci = Ci
        self.__Aj = Aj
        self.__Cj = Cj
        self.__i = i
        self.__j = j
############################################################################
    def ABC(self): return self.__ABC
    def AC(self): return self.__ABC[:self.__Ai], self.__ABC[self.__Ci:],self.__ABC[:self.__Aj], self.__ABC[self.__Cj:]
    def ij(self): return self.__i, self.__j
    def i(self): return self.__i
    def j(self): return self.__j
    def degree(self): return len(self.__ABC)
############################################################################     
    def __eq__(self,other):
       return self.__i == other.__i and \
              self.__j == other.__j and \
              self.__Ai == other.__Ai and \
              self.__Aj == other.__Aj and \
              self.__Ci == other.__Ci and \
              self.__Cj == other.__Cj and \
              self.__ABC == other.__ABC
 ############################################################################                  
    def __ne__(self,other): return not (self == other)
#############################################################################
    def __hash__(self):
       return hash( (self.__ABC,self.__Ai,self.__Ci,self.__Aj,self.__Cj,self.__i,self.__j) )
#############################################################################
    def __repr__(self):
            Ai,Ci,Aj,Cj = self.AC()
            return "(" + self.__ABC + ", " + Ai + ", " + Ci + ", " + Aj + ", " + Cj + ", (" + str(self.__i) + ", " + str(self.__j) + "))"
############################################################################
    def __truediv__(self,other):
        r"""
        Divide ``self`` by ``other`` if possible.
        
        Am ambiguity `(ABC, A_i, C_i, A_j, C_j, i ,j)` is divisible by 
        another ambiguity `(A'B'C', A_i', C_i', A_j', C_j', i', j')` if
        there exist `L` and `R` such that `A_j = L A_j' and `C_j = C_j' R`.
        
        INPUT:
        
        - ``other`` -- Ambiguity
        
        OUTPUT:
        
        - `0` if ``other`` does not divides ``self``
        
        - `1` if ``other`` divides ``self`` and `L = R = 1`
        
        - `2` if if ``other`` divides ``self`` and `LR \neq 1`, i.e, ``other``
        properly divides ``self`` 
        
        TESTS::
        
            sage: a = Ambiguity()
        
        
        """  
        
        sAj, sCj = self.__Aj, self.__Cj
        oAj, oCj = other.__Aj, other.__Cj
        sdeg = self.degree()
        odeg = other.degree()
        
        start = sAj - oAj
        end = start + odeg
        
        
        
        if start < 0 or end > sdeg: return 0
        elif self.__ABC[start:end] == other.__ABC:
            if sdeg == odeg: return 1
            else: return 2
############################################################################    
    @staticmethod
    def generate_incls(a, b):        
        amb = []
        
        i,v = a
        j,w = b
        
        if i == j: return amb
           
        k = v.find(w, 0)
        while k >= 0:
            amb.append( Ambiguity(v,0,len(v),k,k+len(w),i,j) )
            k = v.find(w,k+1)

        return amb
            
############################################################################
    @staticmethod
    def generate_with_tries(prefix_trie, suffix_trie, words, i, criterion=True):
          
        amb = []
        m = words[i]
        m_rev = m[::-1]
        len_m = len(m)

        for k in range(1,len_m):
            # overlap ambiguities with m = AB
            A = m[:-k]
            B = m[-k:]
            amb += [Ambiguity(A+BC,len(A),len(A+BC),0,len_m,j,i) for j,BC in prefix_trie.values(B) if len(BC) > k and j <= i]

            # overlap ambiguities with m = BC
            B = m_rev[-k:]
            C = m[k:]            
            amb += [Ambiguity(AB+C,0,len(AB),len(AB)-k,len(AB+C),j,i) for j,AB in suffix_trie.values(B) if len(AB) > k and j < i]
   
        # inclusion ambiguities with m = ABC       
        amb += [Ambiguity(m,k-len(B)+1,k+1,0,len(m),j,i) for k,(j,B) in prefix_trie.iter(m) if j < i]
        # inclusion ambiguities with m = B
        v = i,m
        amb += flatten([Ambiguity.generate_incls((j,w),v) for j,w in enumerate(words[:i]) if len(w) > len_m])
        
        if criterion:
            amb = Ambiguity.gebauer_moeller(amb)
                
        return amb
############################################################################
    @staticmethod
    def gebauer_moeller(amb):
        amb = sorted(amb,key=lambda a : (a.degree(), a.__i, a.__Ai))
        idx = 0
                                    
        while idx < len(amb)-1:
            a = amb[idx]
            idx += 1
            amb_tmp = amb[:idx]
            j = a.__i
            for aa in amb[idx:]:
                vw = aa / a
                if vw:
                    # divisible and vw != 1, or
                    # divisible and i > j
                    # => throw away aa
                    i = aa.__i
                    if vw > 1 or i > j: continue
                
                    # divisible, i == j, vw = 1 -> look at cofactors of ambiguities
                    # if aa_wi > a_wi -> throw away
                    if i == j and vw == 1 and aa.__Ai > a.__Ai: continue
             
                # not in one of the above cases -> have to keep aa
                amb_tmp.append(aa)
            amb = amb_tmp
                    
        return amb
    
############################################################################
    @staticmethod
    def chain_criterion(amb_old,amb_s,lm_s,s):
                
        amb = []
        set_amb_s = set(amb_s)
        len_lm_s = len(lm_s)
        for a in amb_old:
            keep = True   
            ABC = a.ABC()
            if lm_s in ABC:
                i,j = a.ij()
                Ai,Ci,Aj,Cj = a.__Ai, a.__Ci, a.__Aj, a.__Cj
                k = ABC.find(lm_s)
                len_ABC = len(ABC)
                while k > -1:
                    As,Cs = k, k+len_lm_s
                    # check if new ambiguities are trivial
                    is_trivial_i = Ai >= As + len_lm_s or Ci <= Cs - len_lm_s
                    is_trivial_j = Aj >= As + len_lm_s or Cj <= Cs - len_lm_s
                    # compute new ambiguities and shrink them
                    a_is = Ambiguity(ABC,Ai,Ci,As,Cs,i,s).shrink() if not is_trivial_i else None
                    a_js = Ambiguity(ABC,Aj,Cj,As,Cs,j,s).shrink() if not is_trivial_j else None
                    # check if new ambiguities are trivial or contained
                    if (is_trivial_i or a_is in set_amb_s) and (is_trivial_j or a_js in set_amb_s):
                        keep = False
                        break
                    k = ABC.find(lm_s,k+1)
            if keep: amb.append(a)          
        return amb 
############################################################################       
    def shrink(self):
    
        ABC = self.ABC()
        i,j = self.ij()
        Ai,Ci,Aj,Cj = self.__Ai, self.__Ci, self.__Aj, self.__Cj
        k = min(Ai,Aj)
        l = max(Ci,Cj)
        
        return Ambiguity(ABC[k:l], Ai-k, Ci-k, Aj-k, Cj-k,i,j) 