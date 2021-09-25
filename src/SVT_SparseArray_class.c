/****************************************************************************
 *               Core manipulation of SVT_SparseArray objects               *
 ****************************************************************************/
#include "SVT_SparseArray_class.h"

#include "S4Vectors_interface.h"

#include "Rvector_utils.h"
#include "coerceVector2.h"

#include <limits.h>  /* for INT_MAX */
#include <string.h>  /* for memset(), memcpy() */
#include <stdlib.h>  /* for qsort() */


/* All the atomic types + "list". */
static const SEXPTYPE supported_SVT_Rtypes[] = {
	LGLSXP,   // "logical"
	INTSXP,   // "integer"
	REALSXP,  // "double"
	CPLXSXP,  // "complex"
	RAWSXP,   // "raw"
	STRSXP,   // "character"

	VECSXP    // "list"
};

/* Also checks the supplied 'type'. */
SEXPTYPE _get_Rtype_from_Rstring(SEXP type)
{
	SEXP type0;
	SEXPTYPE Rtype;
	int ntypes, i;

	if (!IS_CHARACTER(type) || LENGTH(type) != 1)
		return 0;
	type0 = STRING_ELT(type, 0);
	if (type0 == NA_STRING)
		return 0;
	Rtype = str2type(CHAR(type0));
	ntypes = sizeof(supported_SVT_Rtypes) / sizeof(SEXPTYPE);
	for (i = 0; i < ntypes; i++)
		if (Rtype == supported_SVT_Rtypes[i])
			return Rtype;
	return 0;
}


/****************************************************************************
 * Low-level utils
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

/* General purpose copy function.
   We only support the 7 SEXP types listed in 'supported_SVT_Rtypes' above. */
static inline int copy_Rvector_elts(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset,
		R_xlen_t nelt)
{
	SEXPTYPE Rtype;
	CopyRVectorElts_FUNType copy_Rvector_elts_FUN;

	Rtype = TYPEOF(in);
	copy_Rvector_elts_FUN = _select_copy_Rvector_elts_FUN(Rtype);
	if (copy_Rvector_elts_FUN == NULL)
		return -1;
	if (TYPEOF(out) != Rtype)
		return -1;
	if (in_offset  + nelt > XLENGTH(in))
		return -1;
	if (out_offset + nelt > XLENGTH(out))
		return -1;
	copy_Rvector_elts_FUN(in, in_offset, out, out_offset, nelt);
	return 0;
}


/****************************************************************************
 * Basic manipulation of "leaf vectors"
 *
 * A "leaf vector" is a vector of offset/value pairs sorted by strictly
 * ascending offset. It is represented by a list of 2 parallel vectors:
 * an integer vector of offsets (i.e. 0-based positions) and a vector
 * (atomic or list) of nonzero values.
 * The length of a leaf vector is always >= 1 and <= INT_MAX.
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

static SEXP alloc_leaf_vector(int lv_len, SEXPTYPE Rtype)
{
	SEXP lv_offs, lv_vals, ans;

	lv_offs = PROTECT(NEW_INTEGER(lv_len));
	lv_vals = PROTECT(allocVector(Rtype, lv_len));
	ans = _new_leaf_vector(lv_offs, lv_vals);
	UNPROTECT(2);
	return ans;
}

static SEXP alloc_and_split_leaf_vector(int lv_len, SEXPTYPE Rtype,
		SEXP *lv_offs, SEXP *lv_vals)
{
	SEXP ans;

	ans = PROTECT(alloc_leaf_vector(lv_len, Rtype));
	_split_leaf_vector(ans, lv_offs, lv_vals);
	UNPROTECT(1);
	return ans;
}

/* Always used in a context where 'offs' and 'vals' are not empty.
   Will sort the supplied offset/value pairs by strictly ascending offset
   before storing them in the "leaf vector".
   The presence of duplicates in 'offs' will trigger an error. */
static SEXP make_leaf_vector(SEXP offs, SEXP vals,
		int *order_buf, unsigned short int *rxbuf1, int *rxbuf2)
{
	int ans_len, k, ret;
	const int *offs_p;
	SEXP ans, ans_offs, ans_vals;

	ans_len = LENGTH(offs);
	if (LENGTH(vals) != ans_len)
		error("S4Arrays internal error in make_leaf_vector():\n"
		      "    LENGTH(offs) != LENGTH(vals)");
	if (ans_len == 1)
		return _new_leaf_vector(offs, vals);

	/* Sort the offset/value pairs by strictly ascending offset. */

	for (k = 0; k < ans_len; k++)
		order_buf[k] = k;
	offs_p = INTEGER(offs);

	ret = sort_ints(order_buf, ans_len, offs_p, 0, 1, rxbuf1, rxbuf2);
	/* Note that ckecking the value returned by sort_ints() is not really
	   necessary here because sort_ints() should never fail when 'rxbuf1'
	   and 'rxbuf2' are supplied (see implementation of _sort_ints() in
	   S4Vectors/src/sort_utils.c for the details). We perform this check
	   nonetheless just to be on the safe side in case the implementation
	   of sort_ints() changes in the future. */
	if (ret < 0)
		error("S4Arrays internal error in make_leaf_vector():\n"
		      "    sort_ints() returned an error");

	/* Check for duplicated offsets. */
	for (k = 1; k < ans_len; k++)
		if (offs_p[order_buf[k]] == offs_p[order_buf[k - 1]])
			error("coercion from COO_SparseArray to "
			      "SVT_SparseArray is not supported when\n"
			      "  the \"nzcoo\" slot of the object to "
			      "coerce contains duplicated coordinates");

	ans = PROTECT(alloc_and_split_leaf_vector(ans_len, TYPEOF(vals),
						  &ans_offs, &ans_vals));
	_copy_selected_Rsubvec_elts(offs, 0, order_buf, ans_offs);
	_copy_selected_Rsubvec_elts(vals, 0, order_buf, ans_vals);
	UNPROTECT(1);
	return ans;
}

/* Always used in a context where 'pos' and 'vals' are not empty. */
static SEXP make_leaf_vector_from_one_based_positions(
		const int *pos, SEXP vals, int maxpos,
		int *order_buf, unsigned short int *rxbuf1, int *rxbuf2)
{
	int ans_len, k, p, *offs_p;
	SEXP offs, ans;

	ans_len = LENGTH(vals);
	offs = PROTECT(NEW_INTEGER(ans_len));
	for (k = 0, offs_p = INTEGER(offs); k < ans_len; k++, offs_p++) {
		p = pos[k];
		if (p == NA_INTEGER || p < 1 || p > maxpos) {
			UNPROTECT(1);
			error("\"nzcoo\" slot of COO_SparseArray object "
			      "to coerce contains invalid coordinates");
		}
		*offs_p = p - 1;
	}
	ans = make_leaf_vector(offs, vals,
			       order_buf, rxbuf1, rxbuf2);
	UNPROTECT(1);
	return ans;
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

/* Returns R_NilValue or a "leaf vector". */
static SEXP make_leaf_vector_from_Rsubvec(
		SEXP Rvector, R_xlen_t subvec_offset, int subvec_len,
		int *offs_buf)
{
	int lv_len;
	SEXP lv, lv_offs, lv_vals;

	lv_len = collect_offsets_of_nonzero_Rsubvec_elts(
			Rvector, subvec_offset, subvec_len,
			offs_buf);

	if (lv_len == 0)
		return R_NilValue;

	lv = PROTECT(alloc_and_split_leaf_vector(lv_len, TYPEOF(Rvector),
						 &lv_offs, &lv_vals));
	memcpy(INTEGER(lv_offs), offs_buf, sizeof(int) * lv_len);
	_copy_selected_Rsubvec_elts(Rvector, subvec_offset, offs_buf, lv_vals);
	UNPROTECT(1);
	return lv;
}

/* Returns R_NilValue or a "leaf vector" of the same length or shorter than
   the input "leaf vector". */
static SEXP remove_zeros_from_leaf_vector(SEXP lv, int *offs_buf)
{
	int lv_len, new_lv_len, k, *p;
	SEXP lv_offs, lv_vals, new_lv, new_lv_offs;

	lv_len = _split_leaf_vector(lv, &lv_offs, &lv_vals);
	new_lv = make_leaf_vector_from_Rsubvec(lv_vals, 0, lv_len, offs_buf);
	if (new_lv == R_NilValue)
		return new_lv;

	PROTECT(new_lv);
	new_lv_offs = VECTOR_ELT(new_lv, 0);
	new_lv_len = LENGTH(new_lv_offs);
	for (k = 0, p = INTEGER(new_lv_offs); k < new_lv_len; k++, p++)
		*p = INTEGER(lv_offs)[*p];
	UNPROTECT(1);
	return new_lv;
}

/* Returns R_NilValue or a "leaf vector" of the same length or shorter than
   the input "leaf vector". */
static SEXP coerce_leaf_vector(SEXP lv, SEXPTYPE new_Rtype,
		int *warn, int *offs_buf)
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
		ans = remove_zeros_from_leaf_vector(ans, offs_buf);
	UNPROTECT(2);
	return ans;
}


/****************************************************************************
 * Basic manipulation of "proto leaf vectors"
 *
 * Used as temporary objects to collect offset/value pairs.
 * Two differences with "leaf vectors":
 *   1. Offset/value pairs can be appended to a "proto leaf vector" as they
 *      are collected.
 *   2. Offset/value pairs don't need to be sorted by strictly ascending
 *      offset.
 */

static SEXP alloc_proto_leaf_vector(int plv_len, SEXPTYPE Rtype)
{
	SEXP plv_offs, plv_vals, plv_nelt, ans;

	plv_offs = PROTECT(NEW_INTEGER(plv_len));
	plv_vals = PROTECT(allocVector(Rtype, plv_len));
	plv_nelt = PROTECT(NEW_INTEGER(1));
	INTEGER(plv_nelt)[0] = 0;

	ans = PROTECT(NEW_LIST(3));
	SET_VECTOR_ELT(ans, 0, plv_offs);
	SET_VECTOR_ELT(ans, 1, plv_vals);
	SET_VECTOR_ELT(ans, 2, plv_nelt);
	UNPROTECT(4);
	return ans;
}

static SEXP alloc_list_of_proto_leaf_vectors(
		const int *plv_lens, int plv_lens_len,
		SEXPTYPE Rtype)
{
	SEXP plvs, plv;
	int i, plv_len;

	plvs = PROTECT(NEW_LIST(plv_lens_len));
	for (i = 0; i < plv_lens_len; i++) {
		plv_len = plv_lens[i];
		if (plv_len != 0) {
			plv = PROTECT( alloc_proto_leaf_vector(plv_len, Rtype));
			SET_VECTOR_ELT(plvs, i, plv);
			UNPROTECT(1);
		}
	}
	UNPROTECT(1);
	return plvs;
}

/* 'plv' must be a "proto leaf vector".
   Returns a value >= 0 if success (1 if the "proto leaf vector" is full,
   0 otherwise), and < 0 if error. */
static inline int append_off_val_pair_to_proto_leaf_vector(SEXP plv,
		int off, SEXP nzvals, int nzvals_offset,
		CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	SEXP plv_offs, plv_vals, plv_nelt;
	int plv_len, nelt, *plv_offs_p;

	plv_offs = VECTOR_ELT(plv, 0);
	plv_vals = VECTOR_ELT(plv, 1);
	plv_nelt = VECTOR_ELT(plv, 2);

	plv_len = LENGTH(plv_offs);
	nelt = INTEGER(plv_nelt)[0];
	if (nelt >= plv_len)
		return -1;
	plv_offs_p = INTEGER(plv_offs);
	plv_offs_p[nelt] = off;
	copy_Rvector_elt_FUN(nzvals, (R_xlen_t) nzvals_offset,
			     plv_vals, (R_xlen_t) nelt);
	return ++INTEGER(plv_nelt)[0] == plv_len;
}

/* Always used in a context where 'plv' is not empty. */
static SEXP make_leaf_vector_from_full_proto_leaf_vector(SEXP plv,
		int *order_buf, unsigned short int *rxbuf1, int *rxbuf2)
{
	SEXP plv_offs, plv_vals;

	plv_offs = VECTOR_ELT(plv, 0);
	plv_vals = VECTOR_ELT(plv, 1);
	return make_leaf_vector(plv_offs, plv_vals,
				order_buf, rxbuf1, rxbuf2);
}


/****************************************************************************
 * C_get_SVT_SparseArray_nzdata_length()
 */

/* Recursive. */
static R_xlen_t REC_sum_leaf_vector_lengths(SEXP SVT, int ndim)
{
	R_xlen_t ans;
	int SVT_len, i;
	SEXP subSVT;

	if (SVT == R_NilValue)
		return 0;

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		return XLENGTH(VECTOR_ELT(SVT, 0));
	}

	/* 'SVT' is a regular node (list). */
	ans = 0;
	SVT_len = LENGTH(SVT);
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		ans += REC_sum_leaf_vector_lengths(subSVT, ndim - 1);
	}
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_get_SVT_SparseArray_nzdata_length(SEXP x_dim, SEXP x_SVT)
{
	R_xlen_t nzdata_len;

	nzdata_len = REC_sum_leaf_vector_lengths(x_SVT, LENGTH(x_dim));
	if (nzdata_len > INT_MAX)
		return ScalarReal((double) nzdata_len);
	return ScalarInteger((int) nzdata_len);
}


/****************************************************************************
 * type() setter
 */

/* Recursive. */
static int REC_set_SVT_type(SEXP SVT, const int *dim, int ndim,
		SEXPTYPE new_Rtype, int *warn, int *offs_buf)
{
	SEXP new_lv, subSVT;
	int SVT_len, is_empty, i, ret;

	if (SVT == R_NilValue)
		return 1;

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		new_lv = coerce_leaf_vector(SVT, new_Rtype, warn, offs_buf);
		if (new_lv == R_NilValue)
			return 1;
		PROTECT(new_lv);
		SET_VECTOR_ELT(SVT, 0, VECTOR_ELT(new_lv, 0));
		SET_VECTOR_ELT(SVT, 1, VECTOR_ELT(new_lv, 1));
		UNPROTECT(1);
		return 0;
	}

	/* 'SVT' is a regular node (list). */
	SVT_len = LENGTH(SVT);

	/* Sanity check (should never fail). */
	if (SVT_len != dim[ndim - 1])
		return -1;

	is_empty = 1;
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		ret = REC_set_SVT_type(subSVT, dim, ndim - 1,
				       new_Rtype, warn, offs_buf);
		if (ret < 0)
			return -1;
		if (ret == 1) {
			SET_VECTOR_ELT(SVT, i, R_NilValue);
		} else {
			is_empty = 0;
		}
	}
	return is_empty;
}

/* --- .Call ENTRY POINT --- */
SEXP C_set_SVT_SparseArray_type(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		SEXP new_type)
{
	SEXPTYPE new_Rtype;
	int warn, *offs_buf, ret;
	SEXP ans;

	new_Rtype = _get_Rtype_from_Rstring(new_type);
	if (new_Rtype == 0)
		error("invalid supplied type");

	if (x_SVT == R_NilValue)
		return x_SVT;

	warn = 0;
	offs_buf = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));
	ans = PROTECT(duplicate(x_SVT));
	ret = REC_set_SVT_type(ans, INTEGER(x_dim), LENGTH(x_dim),
			       new_Rtype, &warn, offs_buf);
	if (ret < 0) {
		UNPROTECT(1);
		error("S4Arrays internal error in "
		      "C_set_SVT_SparseArray_type():\n"
		      "    REC_set_SVT_type() returned an error");
	}
	if (warn)
		_CoercionWarning(warn);
	UNPROTECT(1);
	return ret == 1 ? R_NilValue : ans;
}


/****************************************************************************
 * Going from SVT_SparseArray to ordinary array
 */

/* Recursive. */
static int REC_dump_SVT_to_Rsubarray(SEXP SVT,
		const int *dim, int ndim,
		SEXP Rarray, R_xlen_t subarr_offset, R_xlen_t subarr_len,
		CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	int lv_len, k, SVT_len, i, ret;
	SEXP lv_offs, lv_vals, subSVT;
	const int *p;
	R_xlen_t offset;

	if (SVT == R_NilValue)
		return 0;

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		lv_len = _split_leaf_vector(SVT, &lv_offs, &lv_vals);
		if (lv_len < 0)
			return -1;
		for (k = 0, p = INTEGER(lv_offs); k < lv_len; k++, p++) {
			offset = subarr_offset + *p;
			copy_Rvector_elt_FUN(lv_vals, (R_xlen_t) k,
					     Rarray, offset);
		}
		return 0;
	}

	/* 'SVT' is a regular node (list). */
	SVT_len = LENGTH(SVT);
	if (SVT_len != dim[ndim - 1])
		return -1;

	subarr_len /= SVT_len;
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		ret = REC_dump_SVT_to_Rsubarray(subSVT,
				dim, ndim - 1,
				Rarray, subarr_offset, subarr_len,
				copy_Rvector_elt_FUN);
		if (ret < 0)
			return -1;
		subarr_offset += subarr_len;
	}
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_from_SVT_SparseArray_to_Rarray(SEXP x_dim, SEXP x_dimnames,
		SEXP x_type, SEXP x_SVT)
{
	SEXPTYPE Rtype;
	CopyRVectorElt_FUNType copy_Rvector_elt_FUN;
	SEXP ans;
	int ret;

	Rtype = _get_Rtype_from_Rstring(x_type);
	copy_Rvector_elt_FUN = _select_copy_Rvector_elt_FUN(Rtype);
	if (copy_Rvector_elt_FUN == NULL)
		error("S4Arrays internal error in "
		      "C_from_SVT_SparseArray_to_Rarray():\n"
		      "    SVT_SparseArray object has invalid type");

	ans = PROTECT(_new_Rarray(Rtype, x_dim, x_dimnames));
	ret = REC_dump_SVT_to_Rsubarray(x_SVT,
				INTEGER(x_dim), LENGTH(x_dim),
				ans, 0, XLENGTH(ans),
				copy_Rvector_elt_FUN);
	UNPROTECT(1);
	if (ret < 0)
		error("S4Arrays internal error "
		      "in C_from_SVT_SparseArray_to_Rarray():\n"
		      "    invalid SVT_SparseArray object");
	return ans;
}


/****************************************************************************
 * Going from ordinary array to SVT_SparseArray
 */

/* Recursive. */
static SEXP REC_build_SVT_from_Rsubarray(
		SEXP Rarray, R_xlen_t subarr_offset, R_xlen_t subarr_len,
		const int *dim, int ndim,
		SEXPTYPE ans_Rtype, int *warn, int *offs_buf)
{
	SEXP ans, ans_elt;
	int SVT_len, is_empty, i;

	if (ndim == 1) {
		/* Sanity check (should never fail). */
		if (dim[0] != subarr_len)
			error("S4Arrays internal error "
			      "in REC_build_SVT_from_Rsubarray():\n"
			      "    dim[0] != subarr_len");
		ans = make_leaf_vector_from_Rsubvec(
					Rarray, subarr_offset, dim[0],
					offs_buf);
		if (ans_Rtype == TYPEOF(Rarray) || ans == R_NilValue)
			return ans;
		PROTECT(ans);
		ans = coerce_leaf_vector(ans, ans_Rtype, warn, offs_buf);
		UNPROTECT(1);
		return ans;
	}

	SVT_len = dim[ndim - 1];  /* cannot be 0 so safe to divide below */
	subarr_len /= SVT_len;
	ans = PROTECT(NEW_LIST(SVT_len));
	is_empty = 1;
	for (i = 0; i < SVT_len; i++) {
		ans_elt = REC_build_SVT_from_Rsubarray(
					Rarray, subarr_offset, subarr_len,
					dim, ndim - 1,
					ans_Rtype, warn, offs_buf);
		if (ans_elt != R_NilValue) {
			PROTECT(ans_elt);
			SET_VECTOR_ELT(ans, i, ans_elt);
			UNPROTECT(1);
			is_empty = 0;
		}
		subarr_offset += subarr_len;
	}
	UNPROTECT(1);
	return is_empty ? R_NilValue : ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_build_SVT_from_Rarray(SEXP x, SEXP ans_type)
{
	SEXPTYPE ans_Rtype;
	RVectorEltIsZero_FUNType Rvector_elt_is_zero_FUN;
	int x_ndim, warn, *offs_buf;
	R_xlen_t x_len;
	SEXP x_dim, ans;

	ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	if (ans_Rtype == 0)
		error("invalid requested type");

	Rvector_elt_is_zero_FUN = select_Rvector_elt_is_zero_FUN(TYPEOF(x));
	if (Rvector_elt_is_zero_FUN == NULL)
		error("input array has invalid type");

	x_len = XLENGTH(x);
	if (x_len == 0)  /* means that 'any(dim(x) == 0)' is TRUE */
		return R_NilValue;

	x_dim = GET_DIM(x);  /* does not contain zeros */
	x_ndim = LENGTH(x_dim);
	warn = 0;
	offs_buf = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));
	ans = REC_build_SVT_from_Rsubarray(x, 0, x_len,
					   INTEGER(x_dim), x_ndim,
					   ans_Rtype, &warn, offs_buf);
	if (warn) {
		if (ans != R_NilValue)
			PROTECT(ans);
		_CoercionWarning(warn);
		if (ans != R_NilValue)
			UNPROTECT(1);
	}
	return ans;
}


/****************************************************************************
 * Going from SVT_SparseMatrix to [d|l]gCMatrix
 */

/* Return nb of nonzero values in column. */
static int dump_col_to_CsparseMatrix_slots(SEXP SVT, int col_idx,
		SEXP ans_i, SEXP ans_x, int offset)
{
	SEXP subSVT, lv_offs, lv_vals;
	int lv_len, ret;

	subSVT = VECTOR_ELT(SVT, col_idx);
	if (subSVT == R_NilValue)
		return 0;

	/* 'subSVT' is a "leaf vector". */
	lv_len = _split_leaf_vector(subSVT, &lv_offs, &lv_vals);
	if (lv_len < 0)
		return -1;

	/* Copy 0-based row indices from 'lv_offs' to 'ans_i'. */
	_copy_INTEGER_elts(lv_offs, (R_xlen_t) 0,
			ans_i, (R_xlen_t) offset,
			XLENGTH(lv_offs));

	ret = copy_Rvector_elts(lv_vals, (R_xlen_t) 0,
			ans_x, (R_xlen_t) offset,
			XLENGTH(lv_vals));
	if (ret < 0)
		return -1;

	return lv_len;
}

static int dump_SVT_to_CsparseMatrix_slots(SEXP x_SVT, int x_ncol,
		SEXP ans_p, SEXP ans_i, SEXP ans_x)
{
	int offset, j, nzdata_len;

	INTEGER(ans_p)[0] = 0;
	offset = 0;
	for (j = 0; j < x_ncol; j++) {
		nzdata_len = dump_col_to_CsparseMatrix_slots(x_SVT, j,
						ans_i, ans_x, offset);
		if (nzdata_len < 0)
			return -1;
		offset += nzdata_len;
		INTEGER(ans_p)[j + 1] = offset;
	}
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_from_SVT_SparseMatrix_to_CsparseMatrix(SEXP x_dim,
		SEXP x_type, SEXP x_SVT)
{
	R_xlen_t nzdata_len;
	SEXPTYPE x_Rtype;
	int x_ncol, ret;
	SEXP ans_p, ans_i, ans_x, ans;

	if (LENGTH(x_dim) != 2)
		error("object to coerce to [d|l]gCMatrix "
		      "must have exactly 2 dimensions");

	nzdata_len = REC_sum_leaf_vector_lengths(x_SVT, 2);
	if (nzdata_len > INT_MAX)
		error("SVT_SparseMatrix object contains too many nonzero "
		      "values to be turned into a dgCMatrix or lgCMatrix "
		      "object");

	x_Rtype = _get_Rtype_from_Rstring(x_type);
	if (x_Rtype == 0)
		error("S4Arrays internal error in "
		      "C_from_SVT_SparseMatrix_to_CsparseMatrix():\n"
		      "    SVT_SparseMatrix object has invalid type");

	x_ncol = INTEGER(x_dim)[1];

	ans_i = PROTECT(NEW_INTEGER(nzdata_len));
	ans_x = PROTECT(allocVector(x_Rtype, nzdata_len));
	if (nzdata_len == 0) {
		ans_p = PROTECT(_new_Rvector(INTSXP, (R_xlen_t) x_ncol + 1));
	} else {
		ans_p = PROTECT(NEW_INTEGER(x_ncol + 1));
		ret = dump_SVT_to_CsparseMatrix_slots(x_SVT, x_ncol,
						      ans_p, ans_i, ans_x);
		if (ret < 0) {
			UNPROTECT(3);
			error("S4Arrays internal error in "
			      "C_from_SVT_SparseMatrix_to_CsparseMatrix():\n"
			      "    invalid SVT_SparseMatrix object");
		}
	}

	ans = PROTECT(NEW_LIST(3));
	SET_VECTOR_ELT(ans, 0, ans_p);
	SET_VECTOR_ELT(ans, 1, ans_i);
	SET_VECTOR_ELT(ans, 2, ans_x);
	UNPROTECT(4);
	return ans;
}


/****************************************************************************
 * Going from [d|l]gCMatrix to SVT_SparseMatrix
 */

/* Returns R_NilValue or a "leaf vector" of length <= nzcount. */
static SEXP build_leaf_vector_from_CsparseMatrix_col(SEXP x_i, SEXP x_x,
		int offset, int nzcount,
		SEXPTYPE ans_Rtype, int *warn, int *offs_buf,
		CopyRVectorElts_FUNType copy_Rvector_elts_FUN)
{
	SEXP lv_offs, lv_vals, ans;

	lv_offs = PROTECT(NEW_INTEGER(nzcount));
	/* Copy 0-based row indices from 'x_i' to 'lv_offs'. */
	_copy_INTEGER_elts(x_i, (R_xlen_t) offset,
			lv_offs, (R_xlen_t) 0,
			XLENGTH(lv_offs));
	lv_vals = PROTECT(allocVector(TYPEOF(x_x), nzcount));
	copy_Rvector_elts_FUN(x_x, (R_xlen_t) offset,
			lv_vals, (R_xlen_t) 0,
			XLENGTH(lv_vals));
	ans = _new_leaf_vector(lv_offs, lv_vals);
	if (ans_Rtype == TYPEOF(x_x)) {
		UNPROTECT(2);
		return ans;
	}
	PROTECT(ans);
	ans = coerce_leaf_vector(ans, ans_Rtype, warn, offs_buf);
	UNPROTECT(3);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_build_SVT_from_CsparseMatrix(SEXP x, SEXP ans_type)
{
	SEXP x_x, x_Dim, x_p, x_i, ans, ans_elt;
	CopyRVectorElts_FUNType copy_Rvector_elts_FUN;
	SEXPTYPE ans_Rtype;
	int x_nrow, x_ncol, j, offset, nzcount, warn, *offs_buf, is_empty;

	x_x = GET_SLOT(x, install("x"));
	copy_Rvector_elts_FUN = _select_copy_Rvector_elts_FUN(TYPEOF(x_x));
	if (copy_Rvector_elts_FUN == NULL)
		error("unsupported CsparseMatrix type");

	ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	if (ans_Rtype == 0)
		error("invalid requested type");

	x_Dim = GET_SLOT(x, install("Dim"));
	x_nrow = INTEGER(x_Dim)[0];
	x_ncol = INTEGER(x_Dim)[1];
	x_p = GET_SLOT(x, install("p"));

	if (INTEGER(x_p)[x_ncol] == 0)
		return R_NilValue;

	x_i = GET_SLOT(x, install("i"));

	warn = 0;
	offs_buf = (int *) R_alloc(x_nrow, sizeof(int));
	ans = PROTECT(NEW_LIST(x_ncol));
	is_empty = 1;
	for (j = 0; j < x_ncol; j++) {
		offset = INTEGER(x_p)[j];
		nzcount = INTEGER(x_p)[j + 1] - offset;
		if (nzcount != 0) {
			ans_elt = build_leaf_vector_from_CsparseMatrix_col(
						x_i, x_x,
						offset, nzcount,
						ans_Rtype, &warn, offs_buf,
						copy_Rvector_elts_FUN);
			if (ans_elt != R_NilValue) {
				PROTECT(ans_elt);
				SET_VECTOR_ELT(ans, j, ans_elt);
				UNPROTECT(1);
				is_empty = 0;
			}
		}
	}
	if (warn)
		_CoercionWarning(warn);
	UNPROTECT(1);
	return is_empty ? R_NilValue : ans;
}


/****************************************************************************
 * Going from SVT_SparseArray to COO_SparseArray
 */

static SEXP alloc_nzvals(R_xlen_t nzdata_len, SEXP type)
{
	SEXPTYPE Rtype;

	Rtype = _get_Rtype_from_Rstring(type);
	if (Rtype == 0)
		error("S4Arrays internal error in "
		      "alloc_nzvals():\n"
		      "    SVT_SparseArray object has invalid type");
	return allocVector(Rtype, nzdata_len);
}

/* Recursive. */
static int REC_extract_nzcoo_and_nzvals_from_SVT(SEXP SVT,
		SEXP nzvals, int *nzdata_offset,
		int *nzcoo, int nzcoo_nrow, int nzcoo_ncol,
		int *rowbuf, int rowbuf_offset)
{
	int SVT_len, i, ret, lv_len, k, *p, j;
	SEXP subSVT, lv_offs, lv_vals;

	if (SVT == R_NilValue)
		return 0;

	if (rowbuf_offset > 0) {
		if (!isVectorList(SVT))  // IS_LIST() is broken
			return -1;
		SVT_len = LENGTH(SVT);
		for (i = 0; i < SVT_len; i++) {
			subSVT = VECTOR_ELT(SVT, i);
			rowbuf[rowbuf_offset] = i + 1;
			ret = REC_extract_nzcoo_and_nzvals_from_SVT(
					subSVT,
					nzvals, nzdata_offset,
					nzcoo, nzcoo_nrow, nzcoo_ncol,
					rowbuf, rowbuf_offset - 1);
			if (ret < 0)
				return -1;
		}
		return 0;
	}

	/* 'SVT' is a "leaf vector". */
	lv_len = _split_leaf_vector(SVT, &lv_offs, &lv_vals);
	if (lv_len < 0)
		return -1;

	ret = copy_Rvector_elts(lv_vals, (R_xlen_t) 0,
				nzvals, (R_xlen_t) *nzdata_offset,
				XLENGTH(lv_vals));
	if (ret < 0)
		return -1;

	for (k = 0; k < lv_len; k++) {
		rowbuf[0] = INTEGER(lv_offs)[k] + 1;

		/* Copy 'rowbuf' to 'nzcoo'. */
		p = nzcoo + *nzdata_offset;
		for (j = 0; j < nzcoo_ncol; j++) {
			*p = rowbuf[j];
			p += nzcoo_nrow;
		}

		(*nzdata_offset)++;
	}
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_from_SVT_SparseArray_to_COO_SparseArray(SEXP x_dim,
		SEXP x_type, SEXP x_SVT)
{
	R_xlen_t nzdata_len;
	int nzcoo_nrow, nzcoo_ncol, *rowbuf, nzdata_offset, ret;
	SEXP nzcoo, nzvals, ans;

	nzdata_len = REC_sum_leaf_vector_lengths(x_SVT, LENGTH(x_dim));
	if (nzdata_len > INT_MAX)
		error("SVT_SparseArray object contains too many nonzero "
		      "values to be turned into a COO_SparseArray object");

	nzvals = PROTECT(alloc_nzvals(nzdata_len, x_type));

	nzcoo_nrow = (int) nzdata_len;
	nzcoo_ncol = LENGTH(x_dim);
	rowbuf = (int *) R_alloc(nzcoo_ncol, sizeof(int));
	nzcoo = PROTECT(allocMatrix(INTSXP, nzcoo_nrow, nzcoo_ncol));

	nzdata_offset = 0;
	ret = REC_extract_nzcoo_and_nzvals_from_SVT(x_SVT,
			nzvals, &nzdata_offset,
			INTEGER(nzcoo), nzcoo_nrow, nzcoo_ncol,
			rowbuf, nzcoo_ncol - 1);
	if (ret < 0) {
		UNPROTECT(2);
		error("S4Arrays internal error "
		      "in C_from_SVT_SparseArray_to_COO_SparseArray():\n"
		      "    invalid SVT_SparseArray object");
	}

	/* Sanity check (should never fail). */
	if (nzdata_offset != nzcoo_nrow) {
		UNPROTECT(2);
		error("S4Arrays internal error "
		      "in C_from_SVT_SparseArray_to_COO_SparseArray():\n"
		      "    *out_offset != nzcoo_nrow");
	}

	ans = PROTECT(NEW_LIST(2));
	SET_VECTOR_ELT(ans, 0, nzcoo);
	SET_VECTOR_ELT(ans, 1, nzvals);
	UNPROTECT(3);
	return ans;
}


/****************************************************************************
 * Going from COO_SparseArray to SVT_SparseArray
 */

static SEXP check_Mindex_dim(SEXP Mindex, R_xlen_t vals_len, int ndim,
		const char *what1, const char *what2, const char *what3)
{
	SEXP Mindex_dim;

	Mindex_dim = GET_DIM(Mindex);
	if (Mindex_dim == R_NilValue || LENGTH(Mindex_dim) != 2)
		error("'%s' must be a matrix", what1);
	if (!IS_INTEGER(Mindex))
		error("'%s' must be an integer matrix", what1);
	if (INTEGER(Mindex_dim)[0] != vals_len)
		error("nrow(%s) != %s", what1, what2);
	if (INTEGER(Mindex_dim)[1] != ndim)
		error("ncol(%s) != %s", what1, what3);
	return Mindex_dim;
}

static int grow_SVT(SEXP SVT,
		const int *dim, int ndim,
		const int *nzcoo, int nzdata_len, int nzdata_offset)
{
	const int *p;
	int along, i;
	SEXP subSVT;

	p = nzcoo + nzdata_offset;
	if (*p == NA_INTEGER || *p < 1 || *p > dim[0])
		return -1;

	if (ndim >= 3) {
		p += (size_t) nzdata_len * ndim;
		for (along = ndim - 2; along >= 1; along--) {
			p -= nzdata_len;
			if (*p == NA_INTEGER || *p < 1 || *p > dim[along])
				return -1;
			i = *p - 1;
			subSVT = VECTOR_ELT(SVT, i);
			if (along == 1)
				break;
			/* 'subSVT' is NULL or a list. */
			if (subSVT == R_NilValue) {
				subSVT = PROTECT(NEW_LIST(dim[along]));
				SET_VECTOR_ELT(SVT, i, subSVT);
				UNPROTECT(1);
			}
			SVT = subSVT;
		}
		/* 'subSVT' is NULL or an integer vector of counts. */
		if (subSVT == R_NilValue) {
			subSVT = PROTECT(
				_new_Rvector(INTSXP, (R_xlen_t) dim[along])
			);
			SET_VECTOR_ELT(SVT, i, subSVT);
			UNPROTECT(1);
		}
		SVT = subSVT;
	}

	p = nzcoo + nzdata_offset + nzdata_len;
	if (*p == NA_INTEGER || *p < 1 || *p > LENGTH(SVT))
		return -1;
	INTEGER(SVT)[*p - 1]++;
	return 0;
}

static int store_nzoff_and_nzval_in_SVT(
		const int *nzcoo, int nzdata_len, int nzcoo_ncol,
		SEXP nzvals, int nzdata_offset,
		SEXP SVT,
		CopyRVectorElt_FUNType copy_Rvector_elt_FUN,
		int *order_buf, unsigned short int *rxbuf1, int *rxbuf2)
{
	const int *p;
	int along, i, ret;
	SEXP subSVT;

	if (nzcoo_ncol >= 3) {
		p = nzcoo + nzdata_offset + (size_t) nzdata_len * nzcoo_ncol;
		for (along = nzcoo_ncol - 2; along >= 1; along--) {
			p -= nzdata_len;
			i = *p - 1;
			subSVT = VECTOR_ELT(SVT, i);
			if (along == 1)
				break;
			SVT = subSVT;
		}
		/* 'subSVT' is an integer vector of counts or a list. */
		if (IS_INTEGER(subSVT)) {
			subSVT = PROTECT(
				alloc_list_of_proto_leaf_vectors(
						INTEGER(subSVT),
						LENGTH(subSVT),
						TYPEOF(nzvals))
			);
			SET_VECTOR_ELT(SVT, i, subSVT);
			UNPROTECT(1);
		}
		SVT = subSVT;
	}

	p = nzcoo + nzdata_offset + nzdata_len;
	i = *p - 1;
	subSVT = VECTOR_ELT(SVT, i);

	/* 'subSVT' is a "proto leaf vector". */
	ret = append_off_val_pair_to_proto_leaf_vector(subSVT,
					nzcoo[nzdata_offset] - 1,
					nzvals, nzdata_offset,
					copy_Rvector_elt_FUN);
	if (ret < 0)
		return ret;
	if (ret == 1) {
		/* "Proto leaf vector" 'subSVT' is now full.
		   Turn it into a regular (i.e. non-proto) "leaf vector". */
		subSVT = PROTECT(
			make_leaf_vector_from_full_proto_leaf_vector(subSVT,
						order_buf, rxbuf1, rxbuf2)
		);
		SET_VECTOR_ELT(SVT, i, subSVT);
		UNPROTECT(1);
	}
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_build_SVT_from_COO_SparseArray(
		SEXP x_dim, SEXP x_nzcoo, SEXP x_nzvals,
		SEXP ans_type)
{
	//SEXPTYPE ans_Rtype;
	CopyRVectorElt_FUNType copy_Rvector_elt_FUN;
	int x_ndim, nzdata_len, ans_len, i, ret;
	SEXP ans;
	int *order_buf, *rxbuf2;
	unsigned short int *rxbuf1;

	/* The 'ans_type' argument is currently ignored.
	   For now we take care of honoring the user-requested type at the R
	   level and it's not clear that there's much to gain by handling this
	   at the C level. */
	//ans_Rtype = _get_Rtype_from_Rstring(ans_type);
	//if (ans_Rtype == 0)
	//	error("invalid requested type");

	copy_Rvector_elt_FUN = _select_copy_Rvector_elt_FUN(TYPEOF(x_nzvals));
	if (copy_Rvector_elt_FUN == NULL)
		error("'x@nzvals' has invalid type");

	x_ndim = LENGTH(x_dim);
	nzdata_len = LENGTH(x_nzvals);

	/* Check 'x_nzcoo'. */
	check_Mindex_dim(x_nzcoo, (R_xlen_t) nzdata_len, x_ndim,
			 "x@nzcoo", "length(x@nzvals)", "length(x@dim)");

	if (nzdata_len == 0)
		return R_NilValue;

	order_buf = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));
	rxbuf1 = (unsigned short int *)
			R_alloc(INTEGER(x_dim)[0], sizeof(unsigned short int));
	rxbuf2 = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));

	if (x_ndim == 1)
		return make_leaf_vector_from_one_based_positions(
				INTEGER(x_nzcoo), x_nzvals, INTEGER(x_dim)[0],
				order_buf, rxbuf1, rxbuf2);

	ans_len = INTEGER(x_dim)[x_ndim - 1];

	/* 1st pass: Grow the branches of the tree but don't add any
	   leaf vectors to it, only compute their lengths. */
	if (x_ndim == 2) {
		ans = PROTECT(_new_Rvector(INTSXP, (R_xlen_t) ans_len));
	} else {
		ans = PROTECT(NEW_LIST(ans_len));
	}
	for (i = 0; i < nzdata_len; i++) {
		ret = grow_SVT(ans,
			       INTEGER(x_dim), x_ndim,
			       INTEGER(x_nzcoo), nzdata_len, i);
		if (ret < 0) {
			UNPROTECT(1);
			error("\"nzcoo\" slot of COO_SparseArray object "
			      "to coerce contains invalid coordinates");
		}
	}

	/* 2nd pass: Add the leaf vectors to the tree. */
	if (x_ndim == 2)
		ans = PROTECT(
			alloc_list_of_proto_leaf_vectors(
					INTEGER(ans), ans_len,
					TYPEOF(x_nzvals))
		);
	for (i = 0; i < nzdata_len; i++) {
		ret = store_nzoff_and_nzval_in_SVT(
				INTEGER(x_nzcoo), nzdata_len, x_ndim,
				x_nzvals, i,
				ans,
				copy_Rvector_elt_FUN,
				order_buf, rxbuf1, rxbuf2);
		if (ret < 0) {
			UNPROTECT(1);
			error("S4Arrays internal error in "
			      "C_build_SVT_from_COO_SparseArray():\n"
			      "    store_nzoff_and_nzval_in_SVT() "
			      "returned an unexpected error");
		}
	}

	if (x_ndim == 2)
		UNPROTECT(1);
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * Basic manipulation of "extended leaves"
 *
 * An "extended leaf" is used to temporarily attach a subset of the incoming
 * data (represented by 'Mindex' and 'vals') to a "bottom leaf" (a "bottom
 * leaf" being a leaf located at the deepest possible depth in the SVT, i.e.
 * at depth N - 1 where N is the number of dimensions of the sparse array).
 * An "extended leaf" is **either**:
 *   - A Global Offset Buffer (GOB). Global offsets are offsets w.r.t. the
 *     incoming data i.e. w.r.t. 'Mindex' and 'vals'. A GOB is represented
 *     by an IntAE buffer placed behind an external pointer.
 *   - An "extended leaf vector" i.e. a "leaf vector with a GOB on it".
 *     This is represented as a list of length 3: the 2 list elements of a
 *     regular "leaf vector" (lv_offs + lv_vals) + a GOB.
 */

static SEXP new_GOB()
{
	IntAE *go_buf;

	go_buf = new_IntAE(1, 0, 0);
	return R_MakeExternalPtr(go_buf, R_NilValue, R_NilValue);
}

static SEXP new_extended_leaf_vector(SEXP lv)
{
	SEXP lv_offs, lv_vals, gob, ans;
	int lv_len;

	lv_len = _split_leaf_vector(lv, &lv_offs, &lv_vals);
	if (lv_len < 0)
		error("S4Arrays internal error in "
		      "new_extended_leaf_vector():\n"
		      "    unexpected error");
	gob = PROTECT(new_GOB());
	ans = PROTECT(NEW_LIST(3));
	SET_VECTOR_ELT(ans, 0, lv_offs);
	SET_VECTOR_ELT(ans, 1, lv_vals);
	SET_VECTOR_ELT(ans, 2, gob);
	UNPROTECT(2);
	return ans;
}

/* As a side effect the function also puts a new GOB on 'bottom_leaf' if
   it doesn't have one yet. More precisely:
   - If 'bottom_leaf' is NULL, it gets replaced with a GOB.
   - If 'bottom_leaf' is a "leaf vector", it gets replaced with an "extended
     leaf vector".
*/
static int get_GOB(SEXP leaf_parent, int i, SEXP bottom_leaf, SEXP *gob)
{
	if (bottom_leaf == R_NilValue) {
		*gob = PROTECT(new_GOB());
		SET_VECTOR_ELT(leaf_parent, i, *gob);
		UNPROTECT(1);
		return 0;
	}
	if (TYPEOF(bottom_leaf) == EXTPTRSXP) {
		/* 'bottom_leaf' is a GOB. */
		*gob = bottom_leaf;
		return 0;
	}
	if (!isVectorList(bottom_leaf))
		error("S4Arrays internal error in get_GOB():\n"
		      "    unexpected error");
	/* 'bottom_leaf' is a "leaf vector" or an "extended leaf vector". */
	if (LENGTH(bottom_leaf) == 2) {
		/* 'bottom_leaf' is a "leaf vector". */
		bottom_leaf = PROTECT(new_extended_leaf_vector(bottom_leaf));
		SET_VECTOR_ELT(leaf_parent, i, bottom_leaf);
		UNPROTECT(1);
	} else if (LENGTH(bottom_leaf) != 3) {
		error("S4Arrays internal error in get_GOB():\n"
		      "    unexpected error");
	}
	*gob = VECTOR_ELT(bottom_leaf, 2);
	return 0;
}

static void append_global_offset_to_GOB(SEXP gob, int global_offset)
{
	IntAE *go_buf;

	go_buf = (IntAE *) R_ExternalPtrAddr(gob);
	IntAE_insert_at(go_buf, go_buf->_nelt, global_offset);
	return;
}

static const int *tmp_Mindex;  /* needed by compar_global_offsets() below */

static inline int compar_global_offsets(const void *p1, const void *p2)
{
	int global_offset1, global_offset2, ret;

	global_offset1 = *((const int *) p1);
	global_offset2 = *((const int *) p2);
	ret = tmp_Mindex[global_offset1] - tmp_Mindex[global_offset2];
	if (ret != 0)
		return ret;
	/* Break tie by global offset value so the ordering is "stable". */
	return global_offset1 - global_offset2;
}

static void qsort_GOB(SEXP gob, const int *Mindex)
{
	IntAE *go_buf;

	tmp_Mindex = Mindex;
	go_buf = (IntAE *) R_ExternalPtrAddr(gob);
	qsort(go_buf->elts, go_buf->_nelt, sizeof(int), compar_global_offsets);
	return;
}

static int uniq_GOB(SEXP gob, const int *Mindex)
{
	IntAE *go_buf;
	int *p1, k2;
	const int *p2;

	go_buf = (IntAE *) R_ExternalPtrAddr(gob);
	if (go_buf->_nelt <= 1)
		return go_buf->_nelt;
	p1 = go_buf->elts;
	for (k2 = 1, p2 = p1 + 1; k2 < go_buf->_nelt; k2++, p2++) {
		if (Mindex[*p1] != Mindex[*p2])
			p1++;
		*p1 = *p2;
	}
	return go_buf->_nelt = p1 - go_buf->elts + 1;
}

/* IMPORTANT NOTE: Offset/value pairs where the value is zero are NOT dropped
   at the moment! This means that the function always returns a "leaf vector"
   of the same length as the GOB (which should never be 0).
   FIXME: Offset/value pairs where the value is zero should be dropped.
   Once this is implemented the function **will** sometimes return an
   R_NilValue instead of a "leaf vector". */
static SEXP make_leaf_vector_from_GOB(SEXP gob, SEXP Mindex, SEXP vals, int d)
{
	IntAE *go_buf;
	SEXP ans, ans_offs, ans_vals;
	int *ans_offs_p, k, m;
	const int *Mindex_p, *selection;

	go_buf = (IntAE *) R_ExternalPtrAddr(gob);
	qsort_GOB(gob, INTEGER(Mindex));  /* preserves GOB length */
	uniq_GOB(gob, INTEGER(Mindex));  /* preserves GOB length */
	//if (go_buf->_nelt == 0)  /* can't happen at the moment */
	//	return R_NilValue;
	ans = PROTECT(alloc_and_split_leaf_vector(go_buf->_nelt, TYPEOF(vals),
						  &ans_offs, &ans_vals));
	Mindex_p = INTEGER(Mindex);
	selection = go_buf->elts;
	ans_offs_p = INTEGER(ans_offs);
	for (k = 0; k < go_buf->_nelt; k++, selection++, ans_offs_p++) {
		m = Mindex_p[*selection];
		//if (m == NA_INTEGER || m < 1 || m > d)
		//	error("'Mindex' contains invalid coordinates");
		*ans_offs_p = m - 1;
	}
	_copy_selected_Rsubvec_elts(vals, 0, go_buf->elts, ans_vals);
	UNPROTECT(1);
	return ans;
}

/* Returns R_NilValue or a "leaf vector". */
static SEXP merge_GOB_into_leaf_vector(SEXP xlv, SEXP Mindex, SEXP vals)
{
	error("merge_GOB_into_leaf_vector() is not ready yet");
	return R_NilValue;
}


/****************************************************************************
 * C_subassign_SVT_by_Mindex() and C_subassign_SVT_by_Lindex()
 */

static SEXP shallow_list_copy(SEXP x)
{
	int x_len, i;
	SEXP ans;

	if (!isVectorList(x))  // IS_LIST() is broken
		error("S4Arrays internal error in shallow_list_copy():\n"
		      "    'x' is not a list");
	x_len = LENGTH(x);
	ans = PROTECT(NEW_LIST(x_len));
	for (i = 0; i < x_len; i++)
		SET_VECTOR_ELT(ans, i, VECTOR_ELT(x, i));
	UNPROTECT(1);
	return ans;
}

static int go_to_SVT_bottom(SEXP SVT, SEXP SVT0,
		const int *dim, int ndim,
		const int *M, int vals_len,
		SEXP *leaf_parent, int *idx, SEXP *bottom_leaf)
{
	const int *m_p;
	SEXP subSVT0, subSVT;
	int along, d, m, i;

	m_p = M + (size_t) vals_len * ndim;
	subSVT0 = R_NilValue;
	for (along = ndim - 1; along >= 1; along--) {
		d = dim[along];
		m_p -= vals_len;
		m = *m_p;
		//if (m == NA_INTEGER || m < 1 || m > d)
		//	error("'Mindex' contains invalid coordinates");
		i = m - 1;
		subSVT = VECTOR_ELT(SVT, i);
		if (along == 1)
			break;
		if (SVT0 != R_NilValue)
			subSVT0 = VECTOR_ELT(SVT0, i);
		/* 'subSVT' is NULL or a list. */
		if (subSVT == R_NilValue) {
			subSVT = PROTECT(NEW_LIST(d));
			SET_VECTOR_ELT(SVT, i, subSVT);
			UNPROTECT(1);
		} else if (isVectorList(subSVT) && LENGTH(subSVT) == d) {
			/* Shallow copy **only** if 'subSVT' == corresponding
			   node in original 'SVT0'. */
			if (subSVT == subSVT0) {
				subSVT = PROTECT(shallow_list_copy(subSVT));
				SET_VECTOR_ELT(SVT, i, subSVT);
				UNPROTECT(1);
			}
		} else {
			error("S4Arrays internal error in "
			      "go_to_SVT_bottom():\n"
			      "    unexpected error");
		}
		SVT = subSVT;
		if (SVT0 != R_NilValue)
			SVT0 = subSVT0;
	}
	*leaf_parent = SVT;
	*idx = i;
	*bottom_leaf = subSVT;
	return 0;
}

static int attach_incoming_vals_to_SVT(SEXP SVT, SEXP SVT0,
		const int *dim, int ndim,
		const int *Mindex, SEXP vals)
{
	int vals_len, global_offset, i, ret;
	SEXP leaf_parent, bottom_leaf, gob;

	/* We know 'vals' cannot be a long vector because 'XLENGTH(vals)'
	   went thru check_Mindex_dim(). */
	vals_len = LENGTH(vals);
	for (global_offset = 0; global_offset < vals_len; global_offset++) {
		ret = go_to_SVT_bottom(SVT, SVT0, dim, ndim,
				Mindex + global_offset, vals_len,
				&leaf_parent, &i, &bottom_leaf);
		if (ret < 0)
			return -1;
		ret = get_GOB(leaf_parent, i, bottom_leaf, &gob);
		if (ret < 0)
			return -1;
		append_global_offset_to_GOB(gob, global_offset);
	}
	return 0;
}

/* Recursive. */
static SEXP REC_merge_attached_vals_to_SVT(SEXP SVT,
		const int *dim, int ndim, SEXP Mindex, SEXP vals)
{
	int SVT_len, is_empty, i;
	SEXP subSVT;

	if (SVT == R_NilValue)
		return R_NilValue;

	if (ndim == 1) {
		/* 'SVT' is a bottom leaf (GOB, "leaf vector", or
		   "extended leaf vector"). */
		if (TYPEOF(SVT) == EXTPTRSXP) {
			/* 'SVT' is a GOB. */
			return make_leaf_vector_from_GOB(SVT, Mindex, vals,
							 dim[0]);
		}
		SVT_len = LENGTH(SVT);
		if (SVT_len == 2) {
			/* 'SVT' is a "leaf vector". */
			return SVT;
		}
		if (SVT_len == 3) {
			/* 'SVT' is an "extended leaf vector". */
			return merge_GOB_into_leaf_vector(SVT, Mindex, vals);
		}
		error("S4Arrays internal error in "
		      "REC_merge_attached_vals_to_SVT():\n"
		      "    unexpected error");
	}

	/* 'SVT' is a regular node (list). */
	SVT_len = LENGTH(SVT);  /* should be equal to 'd = dim[ndim - 1]' */
	is_empty = 1;
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		subSVT = REC_merge_attached_vals_to_SVT(subSVT,
					dim, ndim - 1, Mindex, vals);
		if (subSVT != R_NilValue) {
			PROTECT(subSVT);
			SET_VECTOR_ELT(SVT, i, subSVT);
			UNPROTECT(1);
			is_empty = 0;
		} else {
			SET_VECTOR_ELT(SVT, i, subSVT);
		}
	}
	return is_empty ? R_NilValue : SVT;
}

/* --- .Call ENTRY POINT --- */
SEXP C_subassign_SVT_by_Mindex(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		SEXP Mindex, SEXP vals)
{
	SEXPTYPE Rtype;
	int x_ndim, d, ret;
	R_xlen_t vals_len;
	SEXP ans;

	Rtype = _get_Rtype_from_Rstring(x_type);
	if (Rtype == 0)
		error("S4Arrays internal error in "
		      "C_subassign_SVT_by_Mindex():\n"
		      "    SVT_SparseMatrix object has invalid type");

	x_ndim = LENGTH(x_dim);
	vals_len = XLENGTH(vals);
	check_Mindex_dim(Mindex, vals_len, x_ndim,
			 "Mindex", "length(vals)", "length(dim(x))");
	if (vals_len == 0)
		return x_SVT;  /* no-op */

	if (x_ndim == 1) {
		/* Will need special treatment like in
		   C_build_SVT_from_COO_SparseArray(). */
		error("'x_ndim == 1' case is not supported yet");
	}

	/* 1st pass */
	d = INTEGER(x_dim)[x_ndim - 1];
	if (x_SVT == R_NilValue) {
		ans = PROTECT(NEW_LIST(d));
	} else if (isVectorList(x_SVT) && LENGTH(x_SVT) == d) {
		ans = PROTECT(shallow_list_copy(x_SVT));
	} else {
		error("unexpected error");
	}
	ret = attach_incoming_vals_to_SVT(ans, x_SVT,
			INTEGER(x_dim), LENGTH(x_dim),
			INTEGER(Mindex), vals);
	if (ret < 0) {
		UNPROTECT(1);
		error("S4Arrays internal error in "
		      "C_subassign_SVT_by_Mindex():\n"
		      "    attach_incoming_vals_to_SVT() returned an error");
	}

	/* 2nd pass */
	ans = REC_merge_attached_vals_to_SVT(ans,
			INTEGER(x_dim), LENGTH(x_dim), Mindex, vals);
	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_subassign_SVT_by_Lindex(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		SEXP Lindex, SEXP vals)
{
	if (!IS_INTEGER(Lindex) || !IS_NUMERIC(Lindex))
		error("'Lindex' must be an integer or numeric vector");
	return R_NilValue;
}


/****************************************************************************
 * Transposition
 */

static void count_nonzero_vals_per_row(SEXP SVT, int nrow, int ncol,
		int *nzcounts)
{
	int j, lv_len, k;
	SEXP subSVT, lv_offs, lv_vals;
	const int *p;

	memset(nzcounts, 0, sizeof(int) * nrow);
	for (j = 0; j < ncol; j++) {
		subSVT = VECTOR_ELT(SVT, j);
		if (subSVT == R_NilValue)
			continue;
		/* 'subSVT' is a "leaf vector". */
		lv_len = _split_leaf_vector(subSVT, &lv_offs, &lv_vals);
		if (lv_len < 0)
			error("S4Arrays internal error in "
			      "count_nonzero_vals_per_row():\n"
			      "    invalid SVT_SparseMatrix object");
		for (k = 0, p = INTEGER(lv_offs); k < lv_len; k++, p++)
			nzcounts[*p]++;
	}
	return;
}

static void **set_quick_out_vals_p(SEXP out_SVT, SEXPTYPE Rtype)
{
	int out_SVT_len, i;
	SEXP lv;

	out_SVT_len = LENGTH(out_SVT);
	switch (Rtype) {
	    case LGLSXP: case INTSXP: {
		int **vals_p, **p;
		vals_p = (int **) R_alloc(out_SVT_len, sizeof(int *));
		for (i = 0, p = vals_p; i < out_SVT_len; i++, p++) {
			lv = VECTOR_ELT(out_SVT, i);
			if (lv != R_NilValue)
				*p = INTEGER(VECTOR_ELT(lv, 1));
		}
		return (void **) vals_p;
	    }
	    case REALSXP: {
		double **vals_p, **p;
		vals_p = (double **) R_alloc(out_SVT_len, sizeof(double *));
		for (i = 0, p = vals_p; i < out_SVT_len; i++, p++) {
			lv = VECTOR_ELT(out_SVT, i);
			if (lv != R_NilValue)
				*p = REAL(VECTOR_ELT(lv, 1));
		}
		return (void **) vals_p;
	    }
	    case CPLXSXP: {
		Rcomplex **vals_p, **p;
		vals_p = (Rcomplex **) R_alloc(out_SVT_len, sizeof(Rcomplex *));
		for (i = 0, p = vals_p; i < out_SVT_len; i++, p++) {
			lv = VECTOR_ELT(out_SVT, i);
			if (lv != R_NilValue)
				*p = COMPLEX(VECTOR_ELT(lv, 1));
		}
		return (void **) vals_p;
	    }
	    case RAWSXP: {
		Rbyte **vals_p, **p;
		vals_p = (Rbyte **) R_alloc(out_SVT_len, sizeof(Rbyte *));
		for (i = 0, p = vals_p; i < out_SVT_len; i++, p++) {
			lv = VECTOR_ELT(out_SVT, i);
			if (lv != R_NilValue)
				*p = RAW(VECTOR_ELT(lv, 1));
		}
		return (void **) vals_p;
	    }
	}
	/* STRSXP and VECSXP cases. */
	return NULL;
}

typedef void (*TransposeCol_FUNType)(int col_idx,
		const int *offs, SEXP lv_vals,
		int **quick_out_offs_p, void **quick_out_vals_p,
		SEXP out_SVT, int *nzcounts);

/* Ignores 'out_SVT' and 'nzcounts'. */
static void transpose_INTEGER_col(int col_idx,
		const int *offs, SEXP lv_vals,
		int **quick_out_offs_p, void **quick_out_vals_p,
		SEXP out_SVT, int *nzcounts)
{
	int **vals_p;
	int lv_len, k, row_idx;
	const int *v;

	vals_p = (int **) quick_out_vals_p;
	lv_len = LENGTH(lv_vals);
	for (k = 0, v = INTEGER(lv_vals); k < lv_len; k++, v++) {
		row_idx = *offs;
		*(quick_out_offs_p[row_idx]++) = col_idx;
		*(vals_p[row_idx]++) = *v;
		offs++;
	}
	return;
}

/* Ignores 'out_SVT' and 'nzcounts'. */
static void transpose_NUMERIC_col(int col_idx,
		const int *offs, SEXP lv_vals,
		int **quick_out_offs_p, void **quick_out_vals_p,
		SEXP out_SVT, int *nzcounts)
{
	double **vals_p;
	int lv_len, k, row_idx;
	const double *v;

	vals_p = (double **) quick_out_vals_p;
	lv_len = LENGTH(lv_vals);
	for (k = 0, v = REAL(lv_vals); k < lv_len; k++, v++) {
		row_idx = *offs;
		*(quick_out_offs_p[row_idx]++) = col_idx;
		*(vals_p[row_idx]++) = *v;
		offs++;
	}
	return;
}

/* Ignores 'out_SVT' and 'nzcounts'. */
static void transpose_COMPLEX_col(int col_idx,
		const int *offs, SEXP lv_vals,
		int **quick_out_offs_p, void **quick_out_vals_p,
		SEXP out_SVT, int *nzcounts)
{
	Rcomplex **vals_p;
	int lv_len, k, row_idx;
	const Rcomplex *v;

	vals_p = (Rcomplex **) quick_out_vals_p;
	lv_len = LENGTH(lv_vals);
	for (k = 0, v = COMPLEX(lv_vals); k < lv_len; k++, v++) {
		row_idx = *offs;
		*(quick_out_offs_p[row_idx]++) = col_idx;
		*(vals_p[row_idx]++) = *v;
		offs++;
	}
	return;
}

/* Ignores 'out_SVT' and 'nzcounts'. */
static void transpose_RAW_col(int col_idx,
		const int *offs, SEXP lv_vals,
		int **quick_out_offs_p, void **quick_out_vals_p,
		SEXP out_SVT, int *nzcounts)
{
	Rbyte **vals_p;
	int lv_len, k, row_idx;
	const Rbyte *v;

	vals_p = (Rbyte **) quick_out_vals_p;
	lv_len = LENGTH(lv_vals);
	for (k = 0, v = RAW(lv_vals); k < lv_len; k++, v++) {
		row_idx = *offs;
		*(quick_out_offs_p[row_idx]++) = col_idx;
		*(vals_p[row_idx]++) = *v;
		offs++;
	}
	return;
}

/* Ignores 'quick_out_vals_p'. */
static void transpose_CHARACTER_col(int col_idx,
		const int *offs, SEXP lv_vals,
		int **quick_out_offs_p, void **quick_out_vals_p,
		SEXP out_SVT, int *nzcounts)
{
	int lv_len, k, row_idx;
	SEXP out_lv;

	lv_len = LENGTH(lv_vals);
	for (k = 0; k < lv_len; k++) {
		row_idx = *offs;
		*(quick_out_offs_p[row_idx]++) = col_idx;
		out_lv = VECTOR_ELT(out_SVT, row_idx);
		_copy_CHARACTER_elt(lv_vals, (R_xlen_t) k,
			VECTOR_ELT(out_lv, 1), (R_xlen_t) nzcounts[row_idx]++);
		offs++;
	}
	return;
}

/* Ignores 'quick_out_vals_p'. */
static void transpose_LIST_col(int col_idx,
		const int *offs, SEXP lv_vals,
		int **quick_out_offs_p, void **quick_out_vals_p,
		SEXP out_SVT, int *nzcounts)
{
	int lv_len, k, row_idx;
	SEXP out_lv;

	lv_len = LENGTH(lv_vals);
	for (k = 0; k < lv_len; k++) {
		row_idx = *offs;
		*(quick_out_offs_p[row_idx]++) = col_idx;
		out_lv = VECTOR_ELT(out_SVT, row_idx);
		_copy_LIST_elt(lv_vals, (R_xlen_t) k,
			VECTOR_ELT(out_lv, 1), (R_xlen_t) nzcounts[row_idx]++);
		offs++;
	}
	return;
}

static TransposeCol_FUNType select_transpose_col_FUN(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case LGLSXP: case INTSXP: return transpose_INTEGER_col;
	    case REALSXP:             return transpose_NUMERIC_col;
	    case CPLXSXP:             return transpose_COMPLEX_col;
	    case RAWSXP:              return transpose_RAW_col;
	    case STRSXP:              return transpose_CHARACTER_col;
	    case VECSXP:              return transpose_LIST_col;
	}
	return NULL;
}

static SEXP transpose_SVT(SEXP SVT, SEXPTYPE Rtype, int nrow, int ncol,
		int *nzcounts)
{
	TransposeCol_FUNType transpose_col_FUN;
	SEXP ans, ans_elt, subSVT, lv_offs, lv_vals;
	int **quick_out_offs_p;
	void **quick_out_vals_p;
	int i, j, lv_len;

	transpose_col_FUN = select_transpose_col_FUN(Rtype);
	if (transpose_col_FUN == NULL)
		error("S4Arrays internal error in "
		      "transpose_SVT():\n"
		      "    SVT_SparseMatrix object has invalid type");

	ans = PROTECT(NEW_LIST(nrow));
	quick_out_offs_p = (int **) R_alloc(nrow, sizeof(int *));
	for (i = 0; i < nrow; i++) {
		lv_len = nzcounts[i];
		if (lv_len != 0) {
			ans_elt = PROTECT(alloc_leaf_vector(lv_len, Rtype));
			SET_VECTOR_ELT(ans, i, ans_elt);
			UNPROTECT(1);
			quick_out_offs_p[i] = INTEGER(VECTOR_ELT(ans_elt, 0));
		}
	}
	quick_out_vals_p = set_quick_out_vals_p(ans, Rtype);

	memset(nzcounts, 0, sizeof(int) * nrow);
	for (j = 0; j < ncol; j++) {
		subSVT = VECTOR_ELT(SVT, j);
		if (subSVT == R_NilValue)
			continue;
		/* 'subSVT' is a "leaf vector". */
		lv_len = _split_leaf_vector(subSVT, &lv_offs, &lv_vals);
		if (lv_len < 0) {
			UNPROTECT(1);
			error("S4Arrays internal error in "
			      "transpose_SVT():\n"
			      "    invalid SVT_SparseMatrix object");
		}
		transpose_col_FUN(j,
			INTEGER(lv_offs), lv_vals,
			quick_out_offs_p, quick_out_vals_p,
			ans, nzcounts);
	}
	UNPROTECT(1);
	return ans;
}

SEXP C_transpose_SVT_SparseMatrix(SEXP x_dim, SEXP x_type, SEXP x_SVT)
{
	SEXPTYPE Rtype;
	int x_nrow, x_ncol, *nzcounts;

	Rtype = _get_Rtype_from_Rstring(x_type);
	if (Rtype == 0)
		error("S4Arrays internal error in "
		      "C_transpose_SVT_SparseMatrix():\n"
		      "    SVT_SparseMatrix object has invalid type");

	if (LENGTH(x_dim) != 2)
		error("object to transpose must have exactly 2 dimensions");

	if (x_SVT == R_NilValue)
		return x_SVT;

	x_nrow = INTEGER(x_dim)[0];
	x_ncol = INTEGER(x_dim)[1];
	nzcounts = (int *) R_alloc(x_nrow, sizeof(int));

	/* 1st pass: Count the number of nonzero values per row in the
	   input object. */
	count_nonzero_vals_per_row(x_SVT, x_nrow, x_ncol, nzcounts);

	/* 2nd pass: Build the transposed SVT. */
	return transpose_SVT(x_SVT, Rtype, x_nrow, x_ncol, nzcounts);
}

