#ifndef _SVTSPARSEARRAY_CLASS_H_
#define _SVTSPARSEARRAY_CLASS_H_

#include <Rdefines.h>

SEXP C_get_SVT_SparseArray_nzdata_length(
	SEXP x_dim,
	SEXP x_SVT
);

SEXP C_set_SVT_SparseArray_type(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP new_type
);

SEXP C_from_SVT_SparseArray_to_Rarray(
	SEXP x_dim,
	SEXP x_dimnames,
	SEXP x_type,
	SEXP x_SVT
);

SEXP C_build_SVT_from_Rarray(
	SEXP x,
	SEXP new_type
);

SEXP C_from_SVT_SparseArray_to_CsparseMatrix(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT
);

SEXP C_build_SVT_from_dgCMatrix(
	SEXP x,
	SEXP new_type
);

SEXP C_from_SVT_SparseArray_to_COO_SparseArray(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT
);

SEXP C_build_SVT_from_COO_SparseArray(
	SEXP x_dim,
	SEXP x_nzcoo,
	SEXP x_nzdata,
	SEXP new_type
);

SEXP C_transpose_SVT_SparseArray(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT
);

#endif  /* _SVTSPARSEARRAY_CLASS_H_ */

