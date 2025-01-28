from __future__ import absolute_import

from sage.matrix.matrix_sparse cimport Matrix_sparse
from sage.matrix.matrix_rational_sparse cimport Matrix_rational_sparse

###################################################
# general auxiliary
###################################################

cpdef set_up_matrix(rows, columns)

cpdef split_along_columns(Matrix_sparse A, n)

cpdef split_along_rows(Matrix_sparse A, m)

###################################################
# rational auxiliary
###################################################

cpdef rational_augment(Matrix_rational_sparse)

###################################################
# general arithmetic
###################################################

cpdef diff_in_place(Matrix_sparse A, Matrix_sparse B)

###################################################
# rational arithmetic
###################################################

cpdef rational_mat_mul(Matrix_rational_sparse self, Matrix_rational_sparse right)

cpdef rational_trsm(Matrix_rational_sparse A, Matrix_rational_sparse B)
