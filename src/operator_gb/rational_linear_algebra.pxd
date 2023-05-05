from __future__ import absolute_import

from sage.matrix.matrix_rational_sparse cimport Matrix_rational_sparse

###################################################
# rational auxiliary
###################################################

cpdef set_up_matrix(rows, columns)

cpdef split_along_columns(Matrix_rational_sparse A, n)

cpdef split_along_rows(Matrix_rational_sparse A, m)

cpdef augment(Matrix_rational_sparse)

###################################################
# rational arithmetic
###################################################

cpdef mat_mul(Matrix_rational_sparse self, Matrix_rational_sparse right)

cpdef trsm(Matrix_rational_sparse A, Matrix_rational_sparse B)

cpdef diff_in_place(Matrix_rational_sparse A, Matrix_rational_sparse B)

###################################################
# rational echelon computations
###################################################

# cpdef _echelon_in_place(Matrix_rational_sparse self, bint transformation=*)

