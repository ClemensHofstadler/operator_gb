# coding: utf-8
"""
OperatorGB
================

Module to compute normal forms of noncommutative polynomials

AUTHOR:

- Clemens Hofstadler (2023-02-12)

"""

#############################################################################
#  Copyright (C) 2020 Clemens Hofstadler (clemens.hofstadler@jku.at).       #
#                                                                           #
#  Distributed under the terms of the GNU General Public License (GPL)      #
#  either version 2, or (at your option) any later version                  #
#                                                                           #
#  http://www.gnu.org/licenses/                                             #
#############################################################################

from copy import deepcopy

from sage.all import matrix,QQ
from sage.matrix.matrix_rational_sparse import Matrix_rational_sparse

from .cofactor import Cofactor
from .free_algebra import MyFreeAlgebra
from .nc_monomial import NCMonomial
from .nc_polynomial import NCPolynomial
from .rational_linear_algebra import set_up_matrix, rational_augment
from .modular_linear_algebra import modn_augment


############################################################################
# Interreduction
############################################################################
def interreduce(G,one_sided=None):
    
    if not G: return G
    if not isinstance(G[0],NCPolynomial):
        raise ValueError("Input has to be of type NCPolynomial")
        
    A = G[0].parent()
    s = len(G)
    H = deepcopy(G)
    
    changes = True
    first_change = 0

    while changes:
        changes = False
        for i in range(first_change,s):
            
            if H[i] is None: continue
            normal_form = __reduced_form__(H,i,A,one_sided=one_sided)
        
            # normal form = g_i
            if normal_form is None: continue
            # normal form = 0
            elif normal_form.is_zero(): H[i].set_to_zero()
            # normal form != 0
            elif normal_form:
                H[i] = normal_form
                if not changes:
                    first_change = i
                    changes = True

    return [g for g in H if not g.is_zero()]
############################################################################
# Normal form computation
############################################################################
def reduced_form(G, f, trace_cofactors=True):
    
    # check input
    if not G: return f
    if not (isinstance(f,NCPolynomial) and isinstance(G[0],NCPolynomial)):
        raise ValueError("Input has to be of type NCPolynomial")
    if f.parent() != G[0].parent():
        raise ValueError("Different parent structures")
    A = f.parent()
        
    if trace_cofactors and not G[0].cofactors():
         for i,g in enumerate(G): g.append_cofactor(Cofactor(1,'',i,'',A))
            
    h = __reduced_form__(G + [f], len(G), f.parent(), intern=False, trace_cofactors=trace_cofactors)
    # h is such that f + h.cofactors() = h
    # make h such that f - h.cofactors() = h
    for c in h.cofactors(): c.multiply_by(-1)
    
    if not (trace_cofactors and isinstance(f,NCPolynomial)):
        h = A(h)
        
    return h
############################################################################
def __reduced_form__(G,i,A,one_sided=None,intern=True,trace_cofactors=True):
        
    f = G[i]
    rows,columns = symbolic_preprocessing_reduction(f,G,i,one_sided=one_sided)
    columns = list(columns)
    columns.sort(reverse=True)
    nr_columns = len(columns)
    
    rows.insert(0,f)
    M = set_up_matrix(rows,{m:i for i,m in enumerate(columns)})
    if isinstance(M,Matrix_rational_sparse):
        M = rational_augment(M)
    else:
        M = modn_augment(M)
    M.echelonize()
    
    row_idx = M.nonzero_positions_in_column(nr_columns)[-1]
    pos = M.nonzero_positions_in_row(row_idx)
    
    coefficients = []
    monomials = []
    c = M[row_idx, nr_columns]
    for j in reversed(pos):
        if j >= nr_columns: continue
        coefficients.append(1/c * M[row_idx,j])
        monomials.append(columns[j])
        
    if monomials:
        normal_form = NCPolynomial(coefficients,monomials)
    else:
        normal_form = NCPolynomial.zero(A) 

    # no reduction
    if intern and normal_form == f:
        return None
    
    # some reduction - take care of cofactors
    if trace_cofactors:
        # treat cofactors of the reducee separately
        t = M[row_idx,nr_columns]
        cofactors = [cof * t for cof in f.cofactors()]
        normal_form.append_cofactors(cofactors)
        # treat cofactors of the reducers
        for j in pos:
            if j <= nr_columns: continue
            t = M[row_idx,j]
            for cofactor in rows[j-nr_columns].cofactors():
                # the coefficient will always be 1, since these are just
                # dummy cofactors to keep track of a,i,b
                _,a,i,b = cofactor.caib()
                cofactors = [cof.lrmul(a,b) * t for cof in G[i].cofactors()]
                normal_form.append_cofactors(cofactors) 
            normal_form.simplify_cofactors()
    
    return normal_form
############################################################################
def symbolic_preprocessing_reduction(f,G,idx,one_sided=None):
    A = f.parent()
    todo = set(f.monomials())
    lm = [(i,g.lm()) for i,g in enumerate(G) if i != idx and not g.is_zero()]
    done = set()
    reducers = []
    while todo:
        t = todo.pop()
        done.add(t)
        agb = find_reducer_reduction(A,G,lm,str(t),one_sided=one_sided)
        if agb:
            reducers.append(agb)
            todo.update(set(agb.monomials()).difference(done))
    return reducers, done
############################################################################
def find_reducer_reduction(A,G,lm,t,one_sided=None):
    
    a = None
    
    if not one_sided:
        i,m = next(((i,m) for i,m in lm if str(m) in t), (None,None))
        if m is not None:
            str_m = str(m)
            if str_m == '': a,b = t,''
            else: a,b = t.split(str_m,1)
    elif one_sided == 'right':
        i,m = next(((i,m) for i,m in lm if t.startswith(str(m))), (None,None))
        if m is not None:
            a,b = '',t[len(m):]
    elif one_sided == 'left':
        raise NotImplementedError
    else:
        raise ValueError("Wrong value for argument one_sided")
    
    if a is None: return None
    
    g = G[i]
    agb = g.lrmul(a,b)
    agb.append_cofactor(Cofactor(1,a,i,b,A))
    return agb

