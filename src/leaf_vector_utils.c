/****************************************************************************
 *                   Basic manipulation of "leaf vectors"                   *
 ****************************************************************************/
#include "leaf_vector_utils.h"

#include "S4Vectors_interface.h"

#include "Rvector_utils.h"
#include "coerceVector2.h"

#include <limits.h>  /* for INT_MAX */
#include <string.h>  /* for memcpy() */


/****************************************************************************
 * collect_offsets_of_nonzero_Rsubvec_elts()
 */

static int collect_offsets_of_nonzero_int_elts(
		const int *in, int in_len, int *offsets)
{
	int *off_p = offsets;
	for (int offset = 0; offset < in_len; offset++, in++)
		if (*in != 0)
			*(off_p++) = offset;
	return (int) (off_p - offsets);
}

static int collect_offsets_of_nonzero_double_elts(
		const double *in, int in_len, int *offsets)
{
	int *off_p = offsets;
	for (int offset = 0; offset < in_len; offset++, in++)
		if (*in != 0.0)
			*(off_p++) = offset;
	return (int) (off_p - offsets);
}

static int collect_offsets_of_nonzero_Rcomplex_elts(
		const Rcomplex *in, int in_len, int *offsets)
{
	int *off_p = offsets;
	for (int offset = 0; offset < in_len; offset++, in++)
		if (in->r != 0.0 || in->i != 0.0)
			*(off_p++) = offset;
	return (int) (off_p - offsets);
}

static int collect_offsets_of_nonzero_Rbyte_elts(
		const Rbyte *in, int in_len, int *offsets)
{
	int *off_p = offsets;
	for (int offset = 0; offset < in_len; offset++, in++)
		if (*in != 0)
			*(off_p++) = offset;
	return (int) (off_p - offsets);
}

typedef int (*RVectorEltIsZero_FUNType)(SEXP x, R_xlen_t i);

/* Compare with 'integer(1)'. */
static inline int INTEGER_elt_is_zero(SEXP x, R_xlen_t i)
{
	return INTEGER(x)[i] == 0;
}

/* Compare with 'numeric(1)'. */
static inline int NUMERIC_elt_is_zero(SEXP x, R_xlen_t i)
{
	return REAL(x)[i] == 0.0;
}

/* Compare with 'complex(1)'. */
static inline int COMPLEX_elt_is_zero(SEXP x, R_xlen_t i)
{
	Rcomplex x_elt;

	x_elt = COMPLEX(x)[i];
	return x_elt.r == 0.0 && x_elt.i == 0.0;
}

/* Compare with 'raw(1)'. */
static inline int RAW_elt_is_zero(SEXP x, R_xlen_t i)
{
	return RAW(x)[i] == 0;
}

/* Compare with 'character(1)'. */
static inline int CHARACTER_elt_is_zero(SEXP x, R_xlen_t i)
{
	SEXP x_elt;

	x_elt = STRING_ELT(x, i);
	return x_elt != NA_STRING && XLENGTH(x_elt) == 0;
}

/* Compare with 'vector("list", length=1)[[1]]' (which is a NULL). */
static inline int LIST_elt_is_zero(SEXP x, R_xlen_t i)
{
	return VECTOR_ELT(x, i) == R_NilValue;
}

static RVectorEltIsZero_FUNType select_Rvector_elt_is_zero_FUN(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case LGLSXP: case INTSXP: return INTEGER_elt_is_zero;
	    case REALSXP:             return NUMERIC_elt_is_zero;
	    case CPLXSXP:             return COMPLEX_elt_is_zero;
	    case RAWSXP:              return RAW_elt_is_zero;
	    case STRSXP:              return CHARACTER_elt_is_zero;
	    case VECSXP:              return LIST_elt_is_zero;
	}
	return NULL;
}

static int collect_offsets_of_nonzero_Rsubvec_elts(
		SEXP Rvector, R_xlen_t subvec_offset, int subvec_len,
		int *offsets)
{
	SEXPTYPE Rtype;
	RVectorEltIsZero_FUNType Rvector_elt_is_zero_FUN;

	Rtype = TYPEOF(Rvector);

	/* Optimized for LGLSXP, INTSXP, REALSXP, CPLXSXP, and RAWSXP. */
	switch (Rtype) {
	    case LGLSXP: case INTSXP:
		return collect_offsets_of_nonzero_int_elts(
				INTEGER(Rvector) + subvec_offset,
				subvec_len, offsets);
	    case REALSXP:
		return collect_offsets_of_nonzero_double_elts(
				REAL(Rvector) + subvec_offset,
				subvec_len, offsets);
	    case CPLXSXP:
		return collect_offsets_of_nonzero_Rcomplex_elts(
				COMPLEX(Rvector) + subvec_offset,
				subvec_len, offsets);
	    case RAWSXP:
		return collect_offsets_of_nonzero_Rbyte_elts(
				RAW(Rvector) + subvec_offset,
				subvec_len, offsets);
	}

	/* STRSXP and VECSXP cases. */
	Rvector_elt_is_zero_FUN = select_Rvector_elt_is_zero_FUN(Rtype);
	if (Rvector_elt_is_zero_FUN == NULL)
		error("S4Arrays internal error in "
		      "collect_offsets_of_nonzero_Rsubvec_elts():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));

	int *off_p = offsets;
	for (int offset = 0; offset < subvec_len; offset++, subvec_offset++)
		if (!Rvector_elt_is_zero_FUN(Rvector, subvec_offset))
			*(off_p++) = offset;
	return (int) (off_p - offsets);
}


/****************************************************************************
 * _new_leaf_vector()
 * _alloc_leaf_vector()
 * _alloc_and_split_leaf_vector()
 */

SEXP _new_leaf_vector(SEXP lv_offs, SEXP lv_vals)
{
	const char *msg;
	R_xlen_t lv_offs_len;
	SEXP ans;

	/* Sanity checks (should never fail). */
	msg = "S4Arrays internal error in _new_leaf_vector():\n"
	      "    invalid 'lv_offs' and/or 'lv_vals' arguments";
	if (!IS_INTEGER(lv_offs))
		error(msg);
	lv_offs_len = XLENGTH(lv_offs);
	if (lv_offs_len != XLENGTH(lv_vals) ||
	    lv_offs_len < 1 || lv_offs_len > INT_MAX)
		error(msg);

	ans = PROTECT(NEW_LIST(2));
	SET_VECTOR_ELT(ans, 0, lv_offs);
	SET_VECTOR_ELT(ans, 1, lv_vals);
	UNPROTECT(1);
	return ans;
}

SEXP _alloc_leaf_vector(int lv_len, SEXPTYPE Rtype)
{
	SEXP lv_offs, lv_vals, ans;

	lv_offs = PROTECT(NEW_INTEGER(lv_len));
	lv_vals = PROTECT(allocVector(Rtype, lv_len));
	ans = _new_leaf_vector(lv_offs, lv_vals);
	UNPROTECT(2);
	return ans;
}

SEXP _alloc_and_split_leaf_vector(int lv_len, SEXPTYPE Rtype,
		SEXP *lv_offs, SEXP *lv_vals)
{
	SEXP ans;

	ans = PROTECT(_alloc_leaf_vector(lv_len, Rtype));
	_split_leaf_vector(ans, lv_offs, lv_vals);
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * _make_leaf_vector_from_Rsubvec()
 *
 * Returns R_NilValue or a "leaf vector".
 */

SEXP _make_leaf_vector_from_Rsubvec(
		SEXP Rvector, R_xlen_t subvec_offset, int subvec_len,
		int *offs_buf)
{
	int ans_len;
	SEXP ans_offs, ans_vals, ans;

	ans_len = collect_offsets_of_nonzero_Rsubvec_elts(
			Rvector, subvec_offset, subvec_len,
			offs_buf);
	if (ans_len == 0)
		return R_NilValue;

	ans_offs = PROTECT(NEW_INTEGER(ans_len));
	memcpy(INTEGER(ans_offs), offs_buf, sizeof(int) * ans_len);

	if (ans_len == subvec_len &&
	    subvec_offset == 0 && XLENGTH(Rvector) == subvec_len)
	{
		/* The full 'Rvector' contains no zeros and can be reused
		   as-is without the need to copy its nonzero values to a
		   new SEXP. */
		ans = _new_leaf_vector(ans_offs, Rvector);
		UNPROTECT(1);
		return ans;
	}

	ans_vals = PROTECT(allocVector(TYPEOF(Rvector), ans_len));
	_copy_selected_Rsubvec_elts(Rvector, subvec_offset, offs_buf, ans_vals);
	ans = _new_leaf_vector(ans_offs, ans_vals);
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * _remove_zeros_from_leaf_vector()
 *
 * Returns R_NilValue or a "leaf vector" of the same length or shorter than
 * the input "leaf vector".
 */

SEXP _remove_zeros_from_leaf_vector(SEXP lv, int *offs_buf)
{
	SEXP lv_offs, lv_vals, ans, ans_offs;
	int lv_len, ans_len;

	lv_len = _split_leaf_vector(lv, &lv_offs, &lv_vals);
	ans = _make_leaf_vector_from_Rsubvec(lv_vals, 0, lv_len, offs_buf);
	if (ans == R_NilValue)
		return ans;

	PROTECT(ans);
	ans_offs = VECTOR_ELT(ans, 0);
	ans_len = LENGTH(ans_offs);
	if (ans_len == lv_len) {
		/* 'lv' contains no zeros. */
		UNPROTECT(1);
		return lv;
	}
	_copy_selected_ints(INTEGER(lv_offs), INTEGER(ans_offs), ans_len,
			    INTEGER(ans_offs));
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * _coerce_leaf_vector()
 *
 * Returns R_NilValue or a "leaf vector" of the same length or shorter than
 * the input "leaf vector".
 */

SEXP _coerce_leaf_vector(SEXP lv, SEXPTYPE new_Rtype, int *warn, int *offs_buf)
{
	SEXP lv_offs, lv_vals, ans_vals, ans;
	SEXPTYPE old_Rtype;
	int can_add_zeros;

	_split_leaf_vector(lv, &lv_offs, &lv_vals);
	ans_vals = PROTECT(_coerceVector2(lv_vals, new_Rtype, warn));
	ans = PROTECT(_new_leaf_vector(lv_offs, ans_vals));
	/* The above coercion can introduce zeros in 'ans_vals' e.g. when
	   going from double/complex to int/raw. We need to remove them. */
	old_Rtype = TYPEOF(lv_vals);
	can_add_zeros = new_Rtype == RAWSXP
		|| (old_Rtype == REALSXP && new_Rtype == INTSXP)
		|| (old_Rtype == CPLXSXP && (new_Rtype == INTSXP ||
					     new_Rtype == REALSXP))
		|| old_Rtype == STRSXP
		|| old_Rtype == VECSXP;
	if (can_add_zeros)
		ans = _remove_zeros_from_leaf_vector(ans, offs_buf);
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * _merge_leaf_vectors()
 *
 * Returns a "leaf vector" whose length is guaranteed to not exceed
 * min(length(lv1) + length(lv2), INT_MAX).
 */

SEXP _merge_leaf_vectors(SEXP lv1, SEXP lv2)
{
	SEXP lv1_offs, lv1_vals, lv2_offs, lv2_vals, ans, ans_offs, ans_vals;
	int lv1_len, lv2_len, ans_len, k1, k2, k, n;
	SEXPTYPE Rtype;
	CopyRVectorElt_FUNType copy_Rvector_elt_FUN;
	CopyRVectorElts_FUNType copy_Rvector_elts_FUN;
	const int *offs1_p, *offs2_p;
	int *ans_offs_p;

	lv1_len = _split_leaf_vector(lv1, &lv1_offs, &lv1_vals);
	lv2_len = _split_leaf_vector(lv2, &lv2_offs, &lv2_vals);
	Rtype = TYPEOF(lv1_vals);
	if (TYPEOF(lv2_vals) != Rtype)
		error("S4Arrays internal error in _merge_leaf_vectors():\n"
		      "    leaf vectors to merge have different types");
	copy_Rvector_elt_FUN = _select_copy_Rvector_elt_FUN(Rtype);
	if (copy_Rvector_elt_FUN == NULL)
		error("S4Arrays internal error in _merge_leaf_vectors():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));
	copy_Rvector_elts_FUN = _select_copy_Rvector_elts_FUN(Rtype);
	if (copy_Rvector_elts_FUN == NULL)
		error("S4Arrays internal error in _merge_leaf_vectors():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));

	offs1_p = INTEGER(lv1_offs);
	offs2_p = INTEGER(lv2_offs);
	ans_len = k1 = k2 = 0;
	while (k1 < lv1_len && k2 < lv2_len) {
		if (*offs1_p < *offs2_p) {
			offs1_p++;
			k1++;
		} else if (*offs1_p > *offs2_p) {
			offs2_p++;
			k2++;
		} else {
			/* *offs1_p == *offs2_p */
			offs1_p++;
			offs2_p++;
			k1++;
			k2++;
		}
		ans_len++;
	}
	if (k1 < lv1_len) {
		ans_len += lv1_len - k1;
	} else if (k2 < lv2_len) {
		ans_len += lv2_len - k2;
	}
        ans = PROTECT(_alloc_and_split_leaf_vector(ans_len, Rtype,
						   &ans_offs, &ans_vals));
	offs1_p = INTEGER(lv1_offs);
	offs2_p = INTEGER(lv2_offs);
	ans_offs_p = INTEGER(ans_offs);
	k = k1 = k2 = 0;
	while (k1 < lv1_len && k2 < lv2_len) {
		if (*offs1_p < *offs2_p) {
			*ans_offs_p = *offs1_p;
			copy_Rvector_elt_FUN(lv1_vals, (R_xlen_t) k1,
					     ans_vals, (R_xlen_t) k);
			offs1_p++;
			k1++;
		} else if (*offs1_p > *offs2_p) {
			*ans_offs_p = *offs2_p;
			copy_Rvector_elt_FUN(lv2_vals, (R_xlen_t) k2,
					     ans_vals, (R_xlen_t) k);
			offs2_p++;
			k2++;
		} else {
			/* *offs1_p == *offs2_p */
			*ans_offs_p = *offs2_p;
			copy_Rvector_elt_FUN(lv2_vals, (R_xlen_t) k2,
					     ans_vals, (R_xlen_t) k);
			offs1_p++;
			offs2_p++;
			k1++;
			k2++;
		}
		ans_offs_p++;
		k++;
	}
	if (k1 < lv1_len) {
		n = lv1_len - k1;
		memcpy(ans_offs_p, offs1_p, sizeof(int) * n);
		copy_Rvector_elts_FUN(lv1_vals, (R_xlen_t) k1,
				      ans_vals, (R_xlen_t) k,
				      (R_xlen_t) n);
	} else if (k2 < lv2_len) {
		n = lv2_len - k2;
		memcpy(ans_offs_p, offs2_p, sizeof(int) * n);
		copy_Rvector_elts_FUN(lv2_vals, (R_xlen_t) k2,
				      ans_vals, (R_xlen_t) k,
				      (R_xlen_t) n);
	}
	UNPROTECT(1);
	return ans;
}

