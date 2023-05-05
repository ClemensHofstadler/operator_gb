# coding: utf-8
"""
Certify
================

Module to automatically prove operator identities.

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
from time import time
from random import choice
import warnings

from sage.all import *

from .free_algebra import MyFreeAlgebra
from .nc_polynomial import NCPolynomial
from .nc_ideal import NCIdeal
from .auxiliary import *
from .normal_form import interreduce, reduced_form
from .quiver import Quiver

############################################################################
# Certify
############################################################################
def certify(assumptions,claims,quiver=None,maxiter=10,maxdeg=-1,order=None,verbose=1):
    
    F = assumptions[0].parent()
    X = F.gens() if not order else order
    A = MyFreeAlgebra(F.base_ring(),X)
    
    # prepare input
    original_claims_was_list = True
    if type(claims) is not list:
        claims = [A(claims)]
        original_claims_was_list = False
    else:
        claims = [A(claim) for claim in claims]  
    assumptions = [A(f) for f in assumptions]
        
    # check compatibility
    if quiver is not None:
        check_compatibility(assumptions,claims,quiver)
    
    # make ideal
    I = NCIdeal(assumptions,order=order)
    
    # main loop Groebner basis computation & reduction
    if verbose > 0:
        print("Computing a (partial) Groebner basis and reducing the claims...\n")
    idx = 0
    reduced_claims = [False] * len(claims)
    
    while idx < maxiter and not all(reduced_claims):
        idx += 1
        if verbose > 0 and idx % 5 == 0:
            print("Starting iteration %d..." % idx)
        G = I.groebner_basis(maxiter=1,maxdeg=maxdeg,reset=False,verbose=verbose-1)

        #try to reduce the claim
        for i,f in enumerate(claims):
            if not reduced_claims[i]:
                h = reduced_form(G,f)
                if h.is_zero():
                    reduced_claims[i] = h.cofactors()
        
    # negative outcome
    if not all(reduced_claims):
        print("Failed! Not all ideal memberships could be verified.")
        return False

    # positive outcome
    if verbose > 0:
        print("Done! Ideal membership of all claims could be verified!")
    
    # check coefficients
    for claim in reduced_claims:
        if not all(c.c() in ZZ for c in claim):
            warnings.warn("There appear non-integer coefficients in the cofactor representations!")

    if original_claims_was_list:
        return reduced_claims
    else:
        return reduced_claims[0]
############################################################################
def check_compatibility(assumptions,claims,Q):
    for f in assumptions:
        if not Q.is_uniformly_compatible(f):
            raise ValueError("The assumption " + f.__repr__() + " is not compatible with the quiver")
    for f in claims:
        if not Q.is_compatible(f):
                raise ValueError("The claim " + f.__repr__() + " is not compatible with the quiver")