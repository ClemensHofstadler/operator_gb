r"""
Useful auxiliary functions

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

from sage.all import ZZ,prod
import re
from itertools import chain

############################################################################
def flatten(l):
    r"""
    Flatten a list of lists.
    
    No checks are performed whether ``l`` is really a list of lists.

    INPUT: 
    
    - ``l`` -- a list of lists
    
    OUTPUT: The flattened list
            
    EXAMPLES::
     
        sage: from OperatorGB import *
        sage: flatten([[1,2],[3,4]])
        [1,2,3,4]
        sage: flatten([[1,2],[],[3,4]])
        [1,2,3,4]
        sage: flatten([[],[]])
        []
    """
    return list(chain(*l))
############################################################################  
def my_pretty_print(string,A):
    r"""
    Pretty-print a string representing a monomial.

    """
    m = string        
    if not m: return '1'
    
    T = A.translator()
    
    X = set(m)
    m = "*".join([m[i:i+1] for i in range(len(m))]) + "*"
    for x in X:
        eq = r"(%s\*){2,}" % x
        m = re.sub(eq, lambda y : str(x) + "^" + str((y.end() - y.start())//2) + "*", m)
    return T(m[:-1],to_internal=False)
############################################################################
def simplify_str(string):
    r"""
    Simplify a string representing a monomial
    
    Remove `*` characters and expand powers `x^n` as `x\dots x`.

    INPUT: 
    
    - ``string`` -- string
    
    OUTPUT: ``string`` with all `*` characters removed and all powers
     `x^n` expanded as `x\dots x`. The string '1' gets simplified to ''.
            
    TESTS::
     
        sage: from OperatorGB import *
        sage: simplify_string('1')
        ''
        sage: simplify_string('a*b*c')
        'abc'
        sage: simplify_string('a^2*b^3*c')
        'aabbbc'
    """
    m = ''
    if string == '1': return m
    else:
        for v in str(string).split('*'):
            j = v.find('^')
            if j != -1:
                m += v[:j] * ZZ(v[j+1:])
            else:
                m += v
        return m
############################################################################
def expand_cofactors(cofactors,gens):
    r"""
    Expand a list of cofactors
    
    Expand a list of cofactors w.r.t. to the polynomials in ``gens``
    and returns the resulting polynomial as an element in the FreeAlgebra.

    INPUT: 
    
    - ``cofactors`` -- a list of cofactors
    
    - ``gens`` -- a list of noncommutative polynomials; elements in the classical
    FreeAlgebra (not element in a MyFreeAlgebra)
    
    OUTPUT: If ``cofactors`` contains elements of the form `(c_j,a_j,i_j,b_j)`
    with coefficient ``c_j``, cofactors ``a_j``, ``b_j`` and non-negative 
    integer ``j_i`` and if ``gens`` contains polynomials `f_i`, then the
    output is the polynomial `\sum_j c_j a_j f_{i_j} b_j`.
     
    TESTS::
     
    """
    F = gens[0].parent()
    out = [c.to_free_algebra(F,gens[c.i()]) for c in cofactors] 
    return sum(map(prod,out))
############################################################################
def pretty_print_proof(proof,assumptions):
    if not proof: return "0 = 0"
    lhs = str(expand_cofactors(proof,assumptions))
    rhs = ""
    F = proof[0]
    for cofactor in proof:
        rhs += cofactor.pretty_print(assumptions[cofactor.i()])
        
    rhs = rhs[1:]
    if rhs[0] == "+": rhs = rhs[2:]
    return lhs + " = " + rhs