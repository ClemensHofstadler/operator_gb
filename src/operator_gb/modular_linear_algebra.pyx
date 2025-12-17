# distutils: language=3

from __future__ import absolute_import

from cysignals.signals cimport sig_on, sig_off
from cysignals.memory cimport sig_malloc, sig_free

from sage.all import matrix
from sage.modules.vector_modn_sparse cimport *


############################################################################
############################################################################
# Linear algebra
############################################################################
############################################################################

cpdef modn_augment(Matrix_modn_sparse A):
    cdef Matrix_modn_sparse B
    cdef Py_ssize_t nr, nc, c, r, k
    cdef int_fast64_t a
    cdef c_vector_modint* v
    
    nr = A._nrows
    nc = A._ncols
    B = A.new_matrix(nr, nr + nc)
                
    sig_on()
    for r from 0 <= r < nr:
       set_entry(&B.rows[r], nc + r, 1)     
       v = &(A.rows[r])

       for k from 0 <= k < v.num_nonzero: 
           c = v.positions[k]
           a = get_entry(v, c)
           set_entry(&B.rows[r], c, a)
    sig_off()
    
    return B
###################################################
# modular arithmetic
###################################################
cpdef modn_mat_mul(Matrix_modn_sparse self, Matrix_modn_sparse right):
    cdef Matrix_modn_sparse ans

    cdef c_vector_modint* v
    cdef c_vector_modint* w
    
    cdef int_fast64_t mod = self.p

    # Build a table that gives the nonzero positions in each column of right
    cdef list nonzero_positions_in_columns = [set() for _ in range(right._ncols)]
    cdef Py_ssize_t i, j, k
    
    for i from 0 <= i < right._nrows:
        v = &(right.rows[i])
        for j from 0 <= j < v.num_nonzero:
            (<set> nonzero_positions_in_columns[v.positions[j]]).add(i)

    ans = self.new_matrix(self._nrows, right._ncols)

    # Now do the multiplication, getting each row complete before filling it in.
    cdef set relevant_columns
    cdef int_fast64_t x, y, s
    cdef Py_ssize_t pos
    
    for i from 0 <= i < self._nrows:
        v = &(self.rows[i])
        if not v.num_nonzero: continue
                
        # find all relevant columns
        relevant_columns = set()
        for k from 0 <= k < v.num_nonzero:
            w = &(right.rows[v.positions[k]])
            (<set> relevant_columns).update({w.positions[p] for p in range(w.num_nonzero)})
                
        # compute dot product with relevant columns
        for j in list(relevant_columns):    
            s = 0
            for k from 0 <= k < v.num_nonzero:
                if v.positions[k] in nonzero_positions_in_columns[j]:
                    y = get_entry(&right.rows[v.positions[k]], j)
                    x = v.entries[k] * y
                    s += x
                    s = s % mod
                    set_entry(&ans.rows[i], j, s)

    return ans
############################################################################
cpdef modn_trsm(Matrix_modn_sparse A, Matrix_modn_sparse B):
    r"""
    Given A,B with A upper triangular, compute BA^{-1}
    """
    cdef Py_ssize_t i, k, nc
    cdef c_vector_modint* a
    cdef c_vector_modint tmp
    cdef Matrix_modn_sparse AT,BT
    cdef int_fast64_t a_minus, mod
    
    mod = A.p
    
    AT = A.transpose()
    BT = B.transpose()
    
    for i from 0 < i < BT._nrows:
        a = &(AT.rows[i])
        if a.num_nonzero == 1: continue
        for k from 0 <= k < a.num_nonzero-1:
            a_minus = mod - a.entries[k]
            add_c_vector_modint_init(&tmp, &BT.rows[i], &BT.rows[a.positions[k]], a_minus)
            clear_c_vector_modint(&BT.rows[i])
            BT.rows[i] = tmp

    return BT.transpose()
