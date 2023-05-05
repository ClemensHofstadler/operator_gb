# coding: utf-8
"""
MonomialOrder
================

Module to set different monomial orderings for Greobner basis computations
in the free algebra.

AUTHOR:

- Clemens Hofstadler (2023-02-09)

"""

############################################################################
# monomial order
############################################################################
def deglex(a,b,X,blocks=None):
    
    la = len(a)
    lb = len(b)
    
    # compare degree
    if la != lb: return la < lb

    # compare lex
    for x,y in zip(a,b):
        if x != y:
            return X[x] < X[y]
    
    # this is a strict order
    return False
############################################################################
def multilex(a, b, X, blocks):

    for block in blocks:
        Xa = sum(x in block for x in a)
        Xb = sum(x in block for x in b)
        if Xa != Xb: return Xa < Xb

    return deglex(a,b,X)