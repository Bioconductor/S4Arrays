#ifndef _SVTSPARSEARRAY_CLASS_H_
#define _SVTSPARSEARRAY_CLASS_H_

#include <Rdefines.h>

SEXPTYPE _get_Rtype_from_Rstring(SEXP type);

SEXP _new_leaf_vector(
	SEXP lv_offs,
	SEXP lv_vals
);

static inline int _split_leaf_vector(SEXP lv, SEXP *lv_offs, SEXP *lv_vals)
{
	R_xlen_t lv_offs_len;

	/* Sanity checks (should never fail). */
	if (!isVectorList(lv))  // IS_LIST() is broken
		return -1;
	if (LENGTH(lv) != 2)
		return -1;
	*lv_offs = VECTOR_ELT(lv, 0);
	*lv_vals = VECTOR_ELT(lv, 1);
	if (!IS_INTEGER(*lv_offs))
		return -1;
	lv_offs_len = XLENGTH(*lv_offs);
	if (lv_offs_len > INT_MAX)
		return -1;
	if (XLENGTH(*lv_vals) != lv_offs_len)
		return -1;
	return (int) lv_offs_len;
}

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
	SEXP ans_type
);

SEXP C_from_SVT_SparseMatrix_to_CsparseMatrix(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT
);

SEXP C_build_SVT_from_CsparseMatrix(
	SEXP x,
	SEXP ans_type
);

SEXP C_from_SVT_SparseArray_to_COO_SparseArray(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT
);

SEXP C_build_SVT_from_COO_SparseArray(
	SEXP x_dim,
	SEXP x_nzcoo,
	SEXP x_nzvals,
	SEXP ans_type
);

SEXP C_transpose_SVT_SparseMatrix(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT
);

#endif  /* _SVTSPARSEARRAY_CLASS_H_ */

