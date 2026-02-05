r"""
Implementation of noncommutative monomials

This class is not relevant for end-users and
should not be used. Instead use the call method
from ``MyFreeAlgebra`` to create noncommutative polynomials.


AUTHORS:

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


from __future__ import absolute_import

import re
from .auxiliary import *
   
############################################################################
############################################################################
# NCMonomial
############################################################################
############################################################################
class NCMonomial:
    r"""
    A basic implementation of noncommutative monomials based on strings.

    .. NOTE::

        This class should not be relevant for end-users. Do not use it.
        Instead use the call method from ``MyFreeAlgebra`` to create
        noncommutative polynomials.
        
    TETS::
     
        sage: from OperatorGB import *
        sage: F.<x,y,z> = FreeAlgebra(QQ,3)
        sage: A = MyFreeAlgebra(QQ,F.gens()) 
        sage: m = A(x*y*z*x).lm()
        sage: n = A(x*z).lm()
        sage: m == n
        False
        sage: m < n
        False
        sage: m >= n
        True
        sage: m.lrmul('',n.mon())
        x*y*z*x^2*z
    """
    
    def __init__(self,monomial,parent):
        r"""
        Constructor for noncommutative monomials.

        .. NOTE::

            End-users are not supposed to call this method.
            Instead use the call method from ``MyFreeAlgebra`` 
            to construct noncommutative polynomials.
        
        TESTS::
    
            sage: from OperatorGB import *
            sage: A = MyFreeAlgebra(QQ,['a','b','c'])
            sage: NCMonomial('',A)
            1
            sage: NCMonomial('abba',A)
            a*b^2*a
            sage: NCMonomial('aaa',A)
            a^3
        """
    
        self.__mon = monomial
        self.__parent = parent
############################################################################    
    def mon(self): return self.__mon
    def parent(self): return self.__parent
############################################################################     
    def __repr__(self):
        return self.pretty_print()
############################################################################
    def __copy__(self): return NCMonomial(self.__mon, self.__parent)
############################################################################
    def __len__(self): return len(self.__mon)
############################################################################    
    def __str__(self): return self.__mon
############################################################################    
    def __lt__(self,other): 
        assert self.parent() == other.parent()
        return self.__parent.cmp(self.__mon,other.__mon)
    def __eq__(self,other):
         return self.__mon == other.__mon
    def __ne__(self,other):
        return self.__mon != other.__mon
    def __le__(self,other): 
        return self == other or self < other
    def __gt__(self,other):
        return not (self <= other)
    def __ge__(self,other):
        return not (self < other)
############################################################################
    def __hash__(self):
         return hash(self.__mon)
############################################################################
    def __mul__(self,other):
        if self.__parent != other.__parent:
            raise ValueError("Monomials have to be defined over the same ring")
        return NCMonomial(self.__mon + other.__mon, self.__parent)
############################################################################    
    def __truediv__(self, other):
        r"""
        Divide ``self`` by ``other`` if possible, or return ``False``.
        
        INPUT:
        
        - ``other`` -- NCMonomial
        
        OUTPUT: Either strings ``a`` and ``b`` such that
        ``other.lrmul(a,b)`` is equal to ``self``, or ``False``
        if no such ``a`` and ``b`` exist.
        
        If ``other`` occurs at several positions in ``self``, then
        the leftmost position is chosen for the computation of ``a``
        and ``b``.
        
        .. NOTE::
        
            End-users are not supposed to call this method.
            Instead use the call method from ``MyFreeAlgebra`` 
            to construct noncommutative polynomials.
        
        TESTS::
    
            sage: from OperatorGB import *
            sage: F.<x,y,z> = FreeAlgebra(QQ,3)
            sage: A = MyFreeAlgebra(QQ,F.gens()) 
            sage: m = A(x*y*z*x).lm()
            sage: n = A(y*z).lm()
            sage: a,b = m / n
            sage: n.lrmul(a,b) == m 
            True
            sage: m = A(x*y*x*y).lm()
            sage: n = A(y*y).lm()
            sage: m / n
            False
            sage: n = A(x*y).lm()
            sage: a,b = m / n
            sage: a
            ''
            sage: n.lrmul(a,b) == m
            True
            sage: a,b = m / m
            sage (a,b) == ('','')
            True
               
        """
        k = self.__mon.find(other.__mon)
        if k < 0: return False
        a = self.__mon[:k]
        b = self.__mon[k+len(other):]
        return a,b
############################################################################
    def lrmul(self, l, r):
        r"""
        Multiply the monomial from the left by ``l`` and from the right 
        by ``r``.
        
        .. NOTE::
        
            End-users are not supposed to call this method.
            Instead use the call method from ``MyFreeAlgebra`` 
            to construct noncommutative polynomials.
        
        INPUT:
        
        - ``l`` -- string; the left cofactor
        - ``r`` -- string; the right cofactor
        
        OUTPUT: The monomial multiplied by the left from ``l`` and 
        by the right from ``r``.
                
        TESTS::
    
            sage: from OperatorGB import *
            sage: F.<x,y,z> = FreeAlgebra(QQ,3)
            sage: A = MyFreeAlgebra(QQ,F.gens()) 
            sage: m = A(x*y*x).lm()
            sage: m.lrmul('','') == m
            True
            sage: m.lrmul('',m.mon())
            x*y*x^2*y*x
        """
        
        out = self.__copy__()
        out.__mon = ''.join([l,self.__mon,r])
        return out
############################################################################
    def divides(self,other):
        return self.__mon in other.__mon
############################################################################
    def pretty_print(self):
        r"""
        Print the monomial in an easily readable way.
        
        EXAMPLES:
    
            sage: from OperatorGB import *
            sage: F.<x,y,z> = FreeAlgebra(QQ,3)
            sage: A = MyFreeAlgebra(QQ,F.gens())
            sage: A(x*y + 1).monomials()
            [1, x*y]
            sage: A(x^2*y*x^3 + y*y*x).monomials()
            [y^2*x, x^2*y*x^3] 
        
        TESTS::
        
            sage: F.<x,xx,xxx> = FreeAlgebra(QQ,3)
            sage: A = MyFreeAlgebra(QQ,F.gens())
            sage: A(x*xx*xx*x*x*xxx).lm() 
            x*xx^2*x^2*xxx
        """
        m = self.__mon
        T = self.__parent.translator()        
        if not m: return '1'
        
        m = "*" + "*".join([m[i:i+1] for i in range(len(m))]) + "*"
        m = T(m,to_internal=False)
        for x in self.__parent.gens():
            eq = r"((?<=\*)%s\*){2,}" % x
            l = len(x) + 1
            m = re.sub(eq, lambda y : str(x) + "^" + str((y.end() - y.start())//l) + "*", m)
        return m[1:-1]