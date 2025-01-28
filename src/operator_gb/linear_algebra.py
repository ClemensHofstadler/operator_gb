from __future__ import absolute_import

#from sage.all import *

from sage.matrix.matrix_rational_sparse import Matrix_rational_sparse

from .rational_linear_algebra import *
from .modular_linear_algebra import *



from copy import copy

  
def faugere_lachartre(rows,columns,m,trace_cofactors=True):
    """
    m = len(pivot_rows)
    """
    
    M = set_up_matrix(rows,columns)

    AB,CD = split_along_rows(M,m)
    A,B = split_along_columns(AB,m)
    C,D = split_along_columns(CD,m)
        
    #compute D - CA^{-1}B
    if isinstance(A,Matrix_rational_sparse):
        CA = rational_trsm(A,C)    
        CAB = rational_mat_mul(CA,B)
    else:
        CA = modn_trsm(A,C)
        CAB = modn_mat_mul(CA,B)
    #assert CAB == CA * B  
    diff_in_place(D,CAB)
    d_cols = D.ncols()
        
    # compute RREF(D-CA^{-1}B)
    if trace_cofactors: 
        if isinstance(D,Matrix_rational_sparse):
            D = rational_augment(D)  
        else:
            D = modn_augment(D)  
    D.echelonize()

    #get transformation matrix
    M,T = split_along_columns(D,d_cols)
    
    if isinstance(D,Matrix_rational_sparse):
        T1 = rational_mat_mul(-T,CA)
    else:
        T1 = modn_mat_mul(-T,CA)
    
    return T,T1,M
    
