# distutils: language=3

from __future__ import absolute_import

from cysignals.signals cimport sig_on, sig_off
from cysignals.memory cimport sig_malloc, sig_free

from sage.all import QQ, matrix
from sage.rings.rational cimport Rational

from sage.modules.vector_rational_sparse cimport *
from sage.libs.gmp.mpq cimport *


############################################################################
############################################################################
# Linear algebra
############################################################################
############################################################################

###################################################
# general auxiliary
###################################################
cpdef set_up_matrix(rows, columns):
    cdef Matrix_sparse A
    cdef Py_ssize_t nr, nc, i, j
    
    K = rows[0].parent().base_ring()
   
    nr = len(rows)
    nc = len(columns)
    A = matrix(K,nr,nc,sparse=True)
                    
    for i,f in enumerate(rows):
        for c,m in zip(f.coefficients(), f.monomials()):
            j = columns[m]
            A.set_unsafe(i,j,c)
    
    return A
############################################################################    
cpdef split_along_rows(Matrix_sparse A, m):
    
    B = A.matrix_from_rows_and_columns(range(m),range(A._ncols))
    C = A.matrix_from_rows_and_columns(range(m,A._nrows),range(A._ncols))    
  
    return B,C
############################################################################        
cpdef split_along_columns(Matrix_sparse A, n):
    
    B = A.matrix_from_rows_and_columns(range(A._nrows),range(n))
    C = A.matrix_from_rows_and_columns(range(A._nrows),range(n,A._ncols)) 
            
    return B,C
###################################################
# rational auxiliary
###################################################
cpdef rational_augment(Matrix_rational_sparse A):
    cdef Matrix_rational_sparse B
    cdef Py_ssize_t nr, nc, c, r, k
    cdef mpq_vector v,w
    cdef mpq_t one, a
    mpq_init(one)
    mpq_init(a)
    mpq_set_ui(one, 1, 1)
    
    nr = A._nrows
    nc = A._ncols
    B = A.new_matrix(nr, nr + nc)
                
    sig_on()
    for r from 0 <= r < nr:
       mpq_vector_set_entry(&B._matrix[r], nc + r, one)     
       v = A._matrix[r]
       w = B._matrix[r]

       for k from 0 <= k < v.num_nonzero: 
           c = v.positions[k]
           mpq_vector_get_entry(a, &v, c)
           mpq_vector_set_entry(&B._matrix[r], c, a)
    sig_off()
    
    mpq_clear(one)
    mpq_clear(a)
    
    return B
###################################################
# general arithmetic
###################################################    
cpdef diff_in_place(Matrix_sparse A, Matrix_sparse B):
    cdef Py_ssize_t i,j
    
    for (i,j) in B.nonzero_positions(copy=False,column_order=False):
        aij = A.get_unsafe(i,j)
        bij = B.get_unsafe(i,j)
        v = aij - bij
        A.set_unsafe(i,j,v)
###################################################
# rational arithmetic
###################################################
cpdef rational_mat_mul(Matrix_rational_sparse self, Matrix_rational_sparse right):
    cdef Matrix_rational_sparse ans

    cdef mpq_vector* v
    cdef mpq_vector* w

    # Build a table that gives the nonzero positions in each column of right
    cdef list nonzero_positions_in_columns = [set() for _ in range(right._ncols)]
    cdef Py_ssize_t i, j, k
    
    for i from 0 <= i < right._nrows:
        v = &(right._matrix[i])
        for j from 0 <= j < v.num_nonzero:
            (<set> nonzero_positions_in_columns[v.positions[j]]).add(i)

    ans = self.new_matrix(self._nrows, right._ncols)

    # Now do the multiplication, getting each row completely before filling it in.
    cdef set relevant_columns
    cdef mpq_t x, y, s
    cdef Py_ssize_t pos
    mpq_init(x)
    mpq_init(y)
    mpq_init(s)
    
    for i from 0 <= i < self._nrows:
        v = &(self._matrix[i])
        if not v.num_nonzero: continue
                
        # find all relevant columns
        relevant_columns = set()
        for k from 0 <= k < v.num_nonzero:
            w = &(right._matrix[v.positions[k]])
            (<set> relevant_columns).update({w.positions[p] for p in range(w.num_nonzero)})
                
        # compute dot product with relevant columns
        for j in list(relevant_columns):    
            mpq_set_si(s, 0, 1)
            for k from 0 <= k < v.num_nonzero:
                if v.positions[k] in nonzero_positions_in_columns[j]:
                    mpq_vector_get_entry(y, &right._matrix[v.positions[k]], j)
                    mpq_mul(x, v.entries[k], y)
                    mpq_add(s, s, x)
            mpq_vector_set_entry(&ans._matrix[i], j, s)

    mpq_clear(x)
    mpq_clear(y)
    mpq_clear(s)
    return ans
############################################################################
cpdef rational_trsm(Matrix_rational_sparse A, Matrix_rational_sparse B):
    r"""
    Given A,B with A upper triangular, compute BA^{-1}
    """
    cdef Py_ssize_t i, k, nc
    cdef mpq_vector tmp
    cdef mpq_t a_minus,minus_one
    cdef Matrix_rational_sparse AT,BT

    mpq_init(a_minus)
    mpq_init(minus_one)
    mpq_set_si(minus_one,-1,1)
    
    AT = A.transpose()
    BT = B.transpose()
    
    for i from 0 < i < BT._nrows:
        a = &(AT._matrix[i])
        if a.num_nonzero == 1: continue
        for k from 0 <= k < a.num_nonzero-1:
            mpq_neg(a_minus, a.entries[k])
            add_mpq_vector_init(&tmp, &BT._matrix[i], &BT._matrix[a.positions[k]], a_minus)
            mpq_vector_clear(&BT._matrix[i])
            BT._matrix[i] = tmp

    mpq_clear(a_minus)
    mpq_clear(minus_one)
    
    return BT.transpose()
