from __future__ import absolute_import

#from sage.all import *

from .rational_linear_algebra import *
from time import time

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
    CA = trsm(A,C)
    #assert CA == C*A.inverse()       
    CAB = mat_mul(CA,B)
    #assert CAB == CA * B  
    diff_in_place(D,CAB)
    d_cols = D.ncols()
        
    # compute RREF(D-CA^{-1}B)
    if trace_cofactors: D = augment(D)    
    D.echelonize()

    #get transformation matrix
    M,T = split_along_columns(D,d_cols)

    T1 = mat_mul(-T,CA)
    #assert T1 == - T*CA
        
    return T,T1,M
    
