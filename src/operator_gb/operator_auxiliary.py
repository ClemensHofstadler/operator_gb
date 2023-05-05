r"""
Useful auxiliary functions for proving operator statements

AUTHORS:

- Clemens Hofstadler (2023-03-01): initial version

"""

# ****************************************************************************
#                          Copyright (C) 2023
#      Clemens Hofstadler(clemens.hofstadler@mathematik.uni-kassel.de)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import absolute_import

from sage.all import FiniteWords, WordMorphism, prod

def pinv(a,b,a_adj,b_adj):
    r"""
    Return the Moore-Penrose equations for ``a``
    
    Return the Moore-Penrose equations for ``a`` with Moore-Penrose inverse ``b``.
    The adjoint of ``a`` is given by ``a_adj`` and the adjoint of ``b`` is given by
    ``b_adj``.
        
    Note that the addition and multiplication of all these elements have to be well-defined.
    
    INPUT:
    
    - ``a`` -- element of a FreeAlgebra
    - ``b`` -- element of a FreeAlgebra
    - ``a_adj`` -- element of a FreeAlgebra
    - ``b_adj`` -- element of a FreeAlgebra
    
    OUTPUT:
    
    The 4 Moore-Penrose equations for ``a`` plus the 2 adjoint statements.
    
    EXAMPLES::
     
        sage: from OperatorGB import *
        sage: F.<x,y,x_adj,y_adj> = FreeAlgebra(QQ,4)
        sage: pinv(x,y,x_adj,y_adj)
        [-x + x*y*x, -y + y*x*y, -x*y + y_adj*x_adj, -y*x + x_adj*y_adj, -x_adj + x_adj*y_adj*x_adj, -y_adj + y_adj*x_adj*y_adj]
        sage: F.<a,b,c,d> = FreeAlgebra(QQ,4)
        sage: pinv(a,a+b,c-d,a)
        [-a + a^3 + a*b*a, -a - b + a^3 + a^2*b + b*a^2 + b*a*b, -a^2 - a*b + a*c - a*d, -a^2 - b*a + c*a - d*a, -c + d + c*a*c - c*a*d - d*a*c + d*a*d, -a + a*c*a - a*d*a]
        
    """
    
    return [a*b*a - a, b*a*b - b, b_adj*a_adj - a*b, a_adj*b_adj - b*a]
 
def get_involution(F,symbols=None):
    # all variables
    all_symbols = F.variable_names()
    # all variables for which involution should be defined
    if not symbols: symbols = all_symbols
    
    d = dict()
    codomain = []
    for v in symbols:
        w = v[:-4] if v.endswith('_adj') else v + '_adj'
        assert w in all_symbols, f"variable {w} is missing" 
        d[v] = w
        codomain.append(w)
        
    domain = FiniteWords(symbols)
    codomain = FiniteWords(codomain)
        
    return WordMorphism(d,domain=domain,codomain=codomain)
  
def adj(f,phi=None):
    F = f.parent()
    # prepare involution endomorphism
    if not phi:
        symbols = list(map(str,f.variables()))
        phi = get_involution(F,symbols)
    
    new_monomials = []
    for m in f.support():
        w = phi(m.to_word()).reversal()
        m = F(str(w.to_monoid_element()))
        new_monomials.append(m)
    
    return sum(map(prod,zip(f.coefficients(),new_monomials)))
        
def add_adj(F):
    out = F
    # prepare involution endomorphism
    for f in F:
        g = adj(f)
        if not (g in out or -g in out): out.append(g)
    
    return out  