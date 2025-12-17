from __future__ import absolute_import

from sage.matrix.matrix_modn_sparse cimport Matrix_modn_sparse

###################################################
# modular auxiliary
###################################################

cpdef modn_augment(Matrix_modn_sparse)

###################################################
# rational arithmetic
###################################################

cpdef modn_mat_mul(Matrix_modn_sparse self, Matrix_modn_sparse right)

cpdef modn_trsm(Matrix_modn_sparse A, Matrix_modn_sparse B)

