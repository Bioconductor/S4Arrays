#ifndef _SVTSPARSEARRAY_CLASS_H_
#define _SVTSPARSEARRAY_CLASS_H_

#include <Rdefines.h>

SEXP C_get_SVT_SparseArray_nzdata_length(
	SEXP x_dim,
	SEXP x_svtree
);

SEXP C_from_SVT_SparseArray_to_COO_SparseArray(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_svtree
);

SEXP C_from_COO_SparseArray_to_SVT_SparseArray(
	SEXP x_dim,
	SEXP x_nzindex,
        SEXP x_nzdata
);

SEXP C_make_SVT_SparseArray_from_dgCMatrix(
	SEXP x,
	SEXP as_integer
);

SEXP C_from_SVT_SparseArray_to_CsparseMatrix(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_svtree
);

SEXP C_from_SVT_SparseArray_to_array(
	SEXP x_dim,
	SEXP x_dimnames,
	SEXP x_type,
	SEXP x_svtree
);

SEXP C_from_array_to_SVT_SparseArray(
	SEXP x
);

#endif  /* _SVTSPARSEARRAY_CLASS_H_ */

