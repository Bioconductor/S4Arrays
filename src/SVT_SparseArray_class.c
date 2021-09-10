/****************************************************************************
 *            Low-level manipulation of SVT_SparseArray objects             *
 ****************************************************************************/
#include "SVT_SparseArray_class.h"

#include "coerceVector2.h"

#include <limits.h>  /* for INT_MAX */
#include <string.h>  /* for memcpy(), memset() */


/* All the atomic types + "list". */
static const SEXPTYPE supported_SVT_Rtypes[] = {
	LGLSXP,   // "logical"
	INTSXP,   // "integer"
	REALSXP,  // "double"
	CPLXSXP,  // "complex"
	STRSXP,   // "character"
	RAWSXP,   // "raw"

	VECSXP    // "list"
};


/****************************************************************************
 * Low-level utils
 */

/* Returns 1-based positions with respect to 'in'. */
static int collect_pos_of_nonzero_int_elts(const int *in, int in_len,
		int *posbuf)
{
	int *p = posbuf;
	for (int pos = 1; pos <= in_len; pos++, in++)
		if (*in != 0)
			*(p++) = pos;  /* 1-based */
	return (int) (p - posbuf);
}

/* Returns 1-based positions with respect to 'in'. */
static int collect_pos_of_nonzero_double_elts(const double *in, int in_len,
		int *posbuf)
{
	int *p = posbuf;
	for (int pos = 1; pos <= in_len; pos++, in++)
		if (*in != 0.0)
			*(p++) = pos;  /* 1-based */
	return (int) (p - posbuf);
}

/* Returns 1-based positions with respect to 'in'. */
static int collect_pos_of_nonzero_Rcomplex_elts(const Rcomplex *in, int in_len,
		int *posbuf)
{
	int *p = posbuf;
	for (int pos = 1; pos <= in_len; pos++, in++)
		if (in->r != 0.0 || in->i != 0.0)
			*(p++) = pos;  /* 1-based */
	return (int) (p - posbuf);
}

/* Returns 1-based positions with respect to 'in'. */
static int collect_pos_of_nonzero_Rbyte_elts(const Rbyte *in, int in_len,
		int *posbuf)
{
	int *p = posbuf;
	for (int pos = 1; pos <= in_len; pos++, in++)
		if (*in != 0)
			*(p++) = pos;  /* 1-based */
	return (int) (p - posbuf);
}

static void copy_selected_int_elts(const int *in,
		const int *in_offsets, int n, int *out)
{
	for (int k = 0; k < n; k++, out++)
		*out = *(in + in_offsets[k]);
	return;
}

static void copy_selected_double_elts(const double *in,
		const int *in_offsets, int n, double *out)
{
	for (int k = 0; k < n; k++, out++)
		*out = *(in + in_offsets[k]);
	return;
}

static void copy_selected_Rcomplex_elts(const Rcomplex *in,
		const int *in_offsets, int n, Rcomplex *out)
{
	for (int k = 0; k < n; k++, out++)
		*out = *(in + in_offsets[k]);
	return;
}

static void copy_selected_Rbyte_elts(const Rbyte *in,
		const int *in_offsets, int n, Rbyte *out)
{
	for (int k = 0; k < n; k++, out++)
		*out = *(in + in_offsets[k]);
	return;
}

static inline size_t get_Rtype_size(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case LGLSXP: case INTSXP: return sizeof(int);
	    case REALSXP:             return sizeof(double);
	    case CPLXSXP:             return sizeof(Rcomplex);
	    case RAWSXP:              return sizeof(Rbyte);
	}
	return 0;
}

/* Like allocVector() but with initialization of the vector elements. */
static SEXP new_Rvector(SEXPTYPE Rtype, R_xlen_t len)
{
	SEXP ans;
	size_t Rtype_size;

	ans = PROTECT(allocVector(Rtype, len));
	/* allocVector() does NOT initialize the vector elements, except
	   for a list or a character vector. */
	if (Rtype != VECSXP && Rtype != STRSXP) {
		Rtype_size = get_Rtype_size(Rtype);
		if (Rtype_size == 0) {
			UNPROTECT(1);
			error("S4Arrays internal error in new_Rvector():\n"
			      "  unsupported 'Rtype'");
		}
		memset(DATAPTR(ans), 0, Rtype_size * XLENGTH(ans));
	}
	UNPROTECT(1);
	return ans;
}

/* Like allocArray() but with initialization of the array elements and
   addition of the dimnames. */
static SEXP new_Rarray(SEXPTYPE Rtype, SEXP dim, SEXP dimnames)
{
	SEXP ans;
	size_t Rtype_size;

	ans = PROTECT(allocArray(Rtype, dim));
	/* allocArray() is just a thin wrapper for allocVector() and the
	   latter does NOT initialize the vector elements, except for a
	   list or a character vector. */
	if (Rtype != VECSXP && Rtype != STRSXP) {
		Rtype_size = get_Rtype_size(Rtype);
		if (Rtype_size == 0) {
			UNPROTECT(1);
			error("S4Arrays internal error in new_Rarray():\n"
			      "  unsupported 'Rtype'");
		}
		memset(DATAPTR(ans), 0, Rtype_size * XLENGTH(ans));
	}
	SET_DIMNAMES(ans, dimnames);
	UNPROTECT(1);
	return ans;
}

typedef int (*RVectorEltIsZero_FUNType)(SEXP x, R_xlen_t i);

typedef void (*CopyRVectorElt_FUNType)(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset);

typedef void (*CopyRVectorElts_FUNType)(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset,
		R_xlen_t nelt);

/* Compare with 'integer(1)'. */
static inline int INTEGER_elt_is_zero(SEXP x, R_xlen_t i)
{
	return INTEGER(x)[i] == 0;
}

static inline void copy_INTEGER_elt(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset)
{
	INTEGER(out)[out_offset] = INTEGER(in)[in_offset];
	return;
}

static inline void copy_INTEGER_elts(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset,
		R_xlen_t nelt)
{
	void *dest, *src;

	dest = INTEGER(out) + out_offset;
	src  = INTEGER(in)  + in_offset;
	memcpy(dest, src, sizeof(int) * nelt);
	return;
}

/* Compare with 'numeric(1)'. */
static inline int NUMERIC_elt_is_zero(SEXP x, R_xlen_t i)
{
	return REAL(x)[i] == 0.0;
}

static inline void copy_NUMERIC_elt(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset)
{
	REAL(out)[out_offset] = REAL(in)[in_offset];
	return;
}

static inline void copy_NUMERIC_elts(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset,
		R_xlen_t nelt)
{
	void *dest, *src;

	dest = REAL(out) + out_offset;
	src  = REAL(in)  + in_offset;
	memcpy(dest, src, sizeof(double) * nelt);
	return;
}

/* Compare with 'complex(1)'. */
static inline int COMPLEX_elt_is_zero(SEXP x, R_xlen_t i)
{
	Rcomplex x_elt;

	x_elt = COMPLEX(x)[i];
	return x_elt.r == 0.0 && x_elt.i == 0.0;
}

static inline void copy_COMPLEX_elt(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset)
{
	COMPLEX(out)[out_offset] = COMPLEX(in)[in_offset];
	return;
}

static inline void copy_COMPLEX_elts(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset,
		R_xlen_t nelt)
{
	void *dest, *src;

	dest = COMPLEX(out) + out_offset;
	src  = COMPLEX(in)  + in_offset;
	memcpy(dest, src, sizeof(Rcomplex) * nelt);
	return;
}

/* Compare with 'raw(1)'. */
static inline int RAW_elt_is_zero(SEXP x, R_xlen_t i)
{
	return RAW(x)[i] == 0;
}

static inline void copy_RAW_elt(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset)
{
	RAW(out)[out_offset] = RAW(in)[in_offset];
	return;
}

static inline void copy_RAW_elts(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset,
		R_xlen_t nelt)
{
	void *dest, *src;

	dest = RAW(out) + out_offset;
	src  = RAW(in)  + in_offset;
	memcpy(dest, src, sizeof(Rbyte) * nelt);
	return;
}

/* Compare with 'character(1)'. */
static inline int CHARACTER_elt_is_zero(SEXP x, R_xlen_t i)
{
	SEXP x_elt;

	x_elt = STRING_ELT(x, i);
	return x_elt != NA_STRING && XLENGTH(x_elt) == 0;
}

static inline void copy_CHARACTER_elt(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset)
{
	SET_STRING_ELT(out, out_offset, STRING_ELT(in, in_offset));
	return;
}

static inline void copy_CHARACTER_elts(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset,
		R_xlen_t nelt)
{
	R_xlen_t k;

	for (k = 0; k < nelt; k++)
		copy_CHARACTER_elt(in, in_offset + k, out, out_offset + k);
	return;
}

/* Compare with 'vector("list", length=1)[[1]]' (which is a NULL). */
static inline int LIST_elt_is_zero(SEXP x, R_xlen_t i)
{
	return isNull(VECTOR_ELT(x, i));
}

static inline void copy_LIST_elt(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset)
{
	SET_VECTOR_ELT(out, out_offset, VECTOR_ELT(in, in_offset));
	return;
}

static inline void copy_LIST_elts(
		SEXP in,  R_xlen_t in_offset,
		SEXP out, R_xlen_t out_offset,
		R_xlen_t nelt)
{
	R_xlen_t k;

	for (k = 0; k < nelt; k++)
		copy_LIST_elt(in, in_offset + k, out, out_offset + k);
	return;
}

static RVectorEltIsZero_FUNType select_Rvector_elt_is_zero_FUN(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case LGLSXP: case INTSXP: return INTEGER_elt_is_zero;
	    case REALSXP:             return NUMERIC_elt_is_zero;
	    case CPLXSXP:             return COMPLEX_elt_is_zero;
	    case RAWSXP:              return RAW_elt_is_zero;
	    case VECSXP:              return LIST_elt_is_zero;
	    case STRSXP:              return CHARACTER_elt_is_zero;
	}
	return NULL;
}

static CopyRVectorElt_FUNType select_copy_Rvector_elt_FUN(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case LGLSXP: case INTSXP: return copy_INTEGER_elt;
	    case REALSXP:             return copy_NUMERIC_elt;
	    case CPLXSXP:             return copy_COMPLEX_elt;
	    case RAWSXP:              return copy_RAW_elt;
	    case VECSXP:              return copy_LIST_elt;
	    case STRSXP:              return copy_CHARACTER_elt;
	}
	return NULL;
}

static CopyRVectorElts_FUNType select_copy_Rvector_elts_FUN(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case LGLSXP: case INTSXP: return copy_INTEGER_elts;
	    case REALSXP:             return copy_NUMERIC_elts;
	    case CPLXSXP:             return copy_COMPLEX_elts;
	    case RAWSXP:              return copy_RAW_elts;
	    case VECSXP:              return copy_LIST_elts;
	    case STRSXP:              return copy_CHARACTER_elts;
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
	copy_Rvector_elts_FUN = select_copy_Rvector_elts_FUN(Rtype);
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

/* Also checks the supplied 'type'. */
static SEXPTYPE get_Rtype_from_Rstring(SEXP type)
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
 * Basic manipulation of "leaf vectors"
 *
 * A "leaf vector" is a sparse vector represented by a list of 2 parallel
 * vectors: an integer vector of positions (1-based) and a vector (atomic
 * or list) of nonzero values.
 * The length of a leaf vector is always <= INT_MAX.
 */

static SEXP new_leaf_vector(SEXP lv_pos, SEXP lv_vals)
{
	const char *msg;
	R_xlen_t lv_pos_len;
	SEXP ans;

	/* Sanity checks (should never fail). */
	msg = "S4Arrays internal error in new_leaf_vector():\n"
	      "  invalid 'lv_pos' and/or 'lv_vals' arguments";
	if (!IS_INTEGER(lv_pos))
		error(msg);
	lv_pos_len = XLENGTH(lv_pos);
	if (lv_pos_len > INT_MAX || lv_pos_len != XLENGTH(lv_vals))
		error(msg);

	ans = PROTECT(NEW_LIST(2));
	SET_VECTOR_ELT(ans, 0, lv_pos);
	SET_VECTOR_ELT(ans, 1, lv_vals);
	UNPROTECT(1);
	return ans;
}

static SEXP alloc_leaf_vector(int lv_len, SEXPTYPE Rtype)
{
	SEXP lv_pos, lv_vals, ans;

	lv_pos  = PROTECT(NEW_INTEGER(lv_len));
	lv_vals = PROTECT(allocVector(Rtype, lv_len));
	ans = new_leaf_vector(lv_pos, lv_vals);
	UNPROTECT(2);
	return ans;
}

static inline int split_leaf_vector(SEXP lv, SEXP *lv_pos, SEXP *lv_vals)
{
	R_xlen_t lv_pos_len;

	/* Sanity checks (should never fail). */
	if (!isVectorList(lv))  // IS_LIST() is broken
		return -1;
	if (LENGTH(lv) != 2)
		return -1;
	*lv_pos = VECTOR_ELT(lv, 0);
	*lv_vals = VECTOR_ELT(lv, 1);
	if (!IS_INTEGER(*lv_pos))
		return -1;
	lv_pos_len = XLENGTH(*lv_pos);
	if (lv_pos_len > INT_MAX)
		return -1;
	if (XLENGTH(*lv_vals) != lv_pos_len)
		return -1;
	return (int) lv_pos_len;
}

static SEXP alloc_and_split_leaf_vector(int lv_len, SEXPTYPE Rtype,
		SEXP *lv_pos, SEXP *lv_vals)
{
	SEXP ans;

	ans = PROTECT(alloc_leaf_vector(lv_len, Rtype));
	split_leaf_vector(ans, lv_pos, lv_vals);
	UNPROTECT(1);
	return ans;
}

/* 'pos' must contain 1-based positions. */
static SEXP make_leaf_vector(const int *pos, SEXP lv_vals, int maxpos)
{
	int lv_len, k, p;
	SEXP lv_pos, ans;

	lv_len = LENGTH(lv_vals);
	lv_pos = PROTECT(NEW_INTEGER(lv_len));
	for (k = 0; k < lv_len; k++) {
		p = pos[k];
		if (p < 1 || p > maxpos) {
			UNPROTECT(1);
			error("the supplied matrix contains "
			      "out-of-bound values");
		}
		INTEGER(lv_pos)[k] = p;
	}
	ans = new_leaf_vector(lv_pos, lv_vals);
	UNPROTECT(1);
	return ans;
}

static int collect_pos_of_nonzero_Rsubvec_elts(
		SEXP Rvector, R_xlen_t subvec_offset, int subvec_len,
		int *posbuf,
		RVectorEltIsZero_FUNType Rvector_elt_is_zero_FUN)
{
	/* Optimized for LGLSXP, INTSXP, REALSXP, CPLXSXP, and RAWSXP. */
	switch (TYPEOF(Rvector)) {
	    case LGLSXP: case INTSXP:
		return collect_pos_of_nonzero_int_elts(
				INTEGER(Rvector) + subvec_offset,
				subvec_len, posbuf);
	    case REALSXP:
		return collect_pos_of_nonzero_double_elts(
				REAL(Rvector) + subvec_offset,
				subvec_len, posbuf);
	    case CPLXSXP:
		return collect_pos_of_nonzero_Rcomplex_elts(
				COMPLEX(Rvector) + subvec_offset,
				subvec_len, posbuf);
	    case RAWSXP:
		return collect_pos_of_nonzero_Rbyte_elts(
				RAW(Rvector) + subvec_offset,
				subvec_len, posbuf);
	}
	/* STRSXP and VECSXP cases. */
	R_xlen_t offset = subvec_offset;
	int *p = posbuf;
	for (int pos = 1; pos <= subvec_len; pos++, offset++)
		if (!Rvector_elt_is_zero_FUN(Rvector, offset))
			*(p++) = pos;  /* 1-based */
	return (int) (p - posbuf);
}

static void copy_selected_Rsubvec_elts(
		SEXP Rvector, R_xlen_t subvec_offset,
		const int *pos, SEXP lv_vals,
		CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	int lv_len;

	lv_len = LENGTH(lv_vals);
	subvec_offset--;  /* so we can use the 1-based positions in 'pos'
			     as offsets below */

	/* Optimized for LGLSXP, INTSXP, REALSXP, CPLXSXP, and RAWSXP. */
	switch (TYPEOF(Rvector)) {
	    case LGLSXP: case INTSXP:
		copy_selected_int_elts(INTEGER(Rvector) + subvec_offset,
				pos, lv_len, INTEGER(lv_vals));
		return;
	    case REALSXP:
		copy_selected_double_elts(REAL(Rvector) + subvec_offset,
				pos, lv_len, REAL(lv_vals));
		return;
	    case CPLXSXP:
		copy_selected_Rcomplex_elts(COMPLEX(Rvector) + subvec_offset,
				pos, lv_len, COMPLEX(lv_vals));
		return;
	    case RAWSXP:
		copy_selected_Rbyte_elts(RAW(Rvector) + subvec_offset,
				pos, lv_len, RAW(lv_vals));
		return;
	}
	/* STRSXP and VECSXP cases. */
	for (int k = 0; k < lv_len; k++) {
		R_xlen_t offset = subvec_offset + pos[k];
		copy_Rvector_elt_FUN(Rvector, offset, lv_vals, k);
	}
	return;
}

/* Returns R_NilValue or a "leaf vector". */
static SEXP make_leaf_vector_from_Rsubvec(
		SEXP Rvector, R_xlen_t subvec_offset, int subvec_len,
		int *posbuf,
		RVectorEltIsZero_FUNType Rvector_elt_is_zero_FUN,
		CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	int lv_len;
	SEXP lv, lv_pos, lv_vals;

	lv_len = collect_pos_of_nonzero_Rsubvec_elts(
			Rvector, subvec_offset, subvec_len,
			posbuf,
			Rvector_elt_is_zero_FUN);

	if (lv_len == 0)
		return R_NilValue;

	lv = PROTECT(alloc_and_split_leaf_vector(lv_len, TYPEOF(Rvector),
						 &lv_pos, &lv_vals));
	/* Fill 'lv_pos'. */
	memcpy(INTEGER(lv_pos), posbuf, sizeof(int) * lv_len);
	/* Fill 'lv_vals'. */
	copy_selected_Rsubvec_elts(Rvector, subvec_offset, posbuf, lv_vals,
				   copy_Rvector_elt_FUN);
	UNPROTECT(1);
	return lv;
}


/****************************************************************************
 * Basic manipulation of "appendable leaf vectors"
 */

static SEXP alloc_appendable_leaf_vector(int alv_len, SEXPTYPE Rtype)
{
	SEXP alv_pos, alv_vals, alv_nelt, ans;

	alv_pos  = PROTECT(NEW_INTEGER(alv_len));
	alv_vals = PROTECT(allocVector(Rtype, alv_len));
	alv_nelt = PROTECT(NEW_INTEGER(1));
	INTEGER(alv_nelt)[0] = 0;

	ans = PROTECT(NEW_LIST(3));
	SET_VECTOR_ELT(ans, 0, alv_pos);
	SET_VECTOR_ELT(ans, 1, alv_vals);
	SET_VECTOR_ELT(ans, 2, alv_nelt);
	UNPROTECT(4);
	return ans;
}

static SEXP alloc_list_of_appendable_leaf_vectors(
		const int *alv_lens, int alv_lens_len,
		SEXPTYPE Rtype)
{
	SEXP alvs, alv;
	int k, alv_len;

	alvs = PROTECT(NEW_LIST(alv_lens_len));
	for (k = 0; k < alv_lens_len; k++) {
		alv_len = alv_lens[k];
		if (alv_len != 0) {
			alv = PROTECT(
				alloc_appendable_leaf_vector(alv_len, Rtype)
			);
			SET_VECTOR_ELT(alvs, k, alv);
			UNPROTECT(1);
		}
	}
	UNPROTECT(1);
	return alvs;
}

/* 'alv' must be an "appendable leaf vector.
   'pos' must be a 1-based position. */
static inline int append_pos_val_pair_to_leaf_vector(SEXP alv,
		int pos, SEXP nzdata, int nzdata_offset,
		CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	SEXP alv_pos, alv_vals, alv_nelt;
	int alv_len, *alv_nelt_p;

	alv_pos  = VECTOR_ELT(alv, 0);
	alv_vals = VECTOR_ELT(alv, 1);
	alv_nelt = VECTOR_ELT(alv, 2);
	alv_len  = LENGTH(alv_pos);
	alv_nelt_p = INTEGER(alv_nelt);
	if (*alv_nelt_p >= alv_len)
		return -1;
	INTEGER(alv_pos)[*alv_nelt_p] = pos;
	copy_Rvector_elt_FUN(nzdata, (R_xlen_t) nzdata_offset,
			     alv_vals, (R_xlen_t) *alv_nelt_p);
	(*alv_nelt_p)++;
	return *alv_nelt_p == alv_len;
}


/****************************************************************************
 * C_get_SVT_SparseArray_nzdata_length()
 */

/* Recursive. */
static R_xlen_t REC_sum_leaf_vector_lengths(SEXP SVT, int ndim)
{
	R_xlen_t ans;
	int SVT_len, k;
	SEXP subSVT;

	if (isNull(SVT))
		return 0;

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		return XLENGTH(VECTOR_ELT(SVT, 0));
	}

	/* 'SVT' is a regular node (list). */
	ans = 0;
	SVT_len = LENGTH(SVT);
	for (k = 0; k < SVT_len; k++) {
		subSVT = VECTOR_ELT(SVT, k);
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
		SEXPTYPE new_Rtype, int *warn)
{
	SEXP lv_vals, subSVT;
	int SVT_len, k, ret;

	if (isNull(SVT))
		return 0;

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		lv_vals = VECTOR_ELT(SVT, 1);
		lv_vals = PROTECT(_coerceVector2(lv_vals, new_Rtype, warn));
		SET_VECTOR_ELT(SVT, 1, lv_vals);
		UNPROTECT(1);
		return 0;
	}

	/* 'SVT' is a regular node (list). */
	SVT_len = LENGTH(SVT);
	if (SVT_len != dim[ndim - 1])
		return -1;
	for (k = 0; k < SVT_len; k++) {
		subSVT = VECTOR_ELT(SVT, k);
		ret = REC_set_SVT_type(subSVT, dim, ndim - 1, new_Rtype, warn);
		if (ret < 0)
			return -1;
	}
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_set_SVT_SparseArray_type(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		SEXP new_type)
{
	SEXPTYPE new_Rtype;
	SEXP ans;
	int warn, ret;

	new_Rtype = get_Rtype_from_Rstring(new_type);
	if (new_Rtype == 0)
		error("invalid supplied type");
	warn = 0;
	ans = PROTECT(duplicate(x_SVT));
	ret = REC_set_SVT_type(ans, INTEGER(x_dim), LENGTH(x_dim),
			       new_Rtype, &warn);
	if (ret < 0) {
		UNPROTECT(1);
		error("S4Arrays internal error in "
		      "C_set_SVT_SparseArray_type():\n"
		      "  REC_set_SVT_type() returned an error");
	}
	if (warn)
		_CoercionWarning(warn);
	UNPROTECT(1);
	return ans;
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
	int lv_len, k, SVT_len, ret;
	SEXP lv_pos, lv_vals, subSVT;
	R_xlen_t offset;

	if (isNull(SVT))
		return 0;

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		lv_len = split_leaf_vector(SVT, &lv_pos, &lv_vals);
		if (lv_len < 0)
			return -1;
		for (k = 0; k < lv_len; k++) {
			offset = subarr_offset + INTEGER(lv_pos)[k] - 1;
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
	for (k = 0; k < SVT_len; k++) {
		subSVT = VECTOR_ELT(SVT, k);
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

	Rtype = get_Rtype_from_Rstring(x_type);
	if (Rtype == 0)
		error("S4Arrays internal error in "
		      "C_from_SVT_SparseArray_to_Rarray():\n"
		      "  SVT_SparseArray object has invalid type");
	copy_Rvector_elt_FUN = select_copy_Rvector_elt_FUN(Rtype);
	ans = PROTECT(new_Rarray(Rtype, x_dim, x_dimnames));
	ret = REC_dump_SVT_to_Rsubarray(x_SVT,
				INTEGER(x_dim), LENGTH(x_dim),
				ans, 0, XLENGTH(ans),
				copy_Rvector_elt_FUN);
	UNPROTECT(1);
	if (ret < 0)
		error("S4Arrays internal error "
		      "in C_from_SVT_SparseArray_to_Rarray():\n"
		      "  invalid SVT_SparseArray object");
	return ans;
}


/****************************************************************************
 * Going from ordinary array to SVT_SparseArray
 */

/* Recursive. */
static SEXP REC_build_SVT_from_Rsubarray(
		SEXP Rarray, R_xlen_t subarr_offset, R_xlen_t subarr_len,
		const int *dim, int ndim,
		int *posbuf,
		RVectorEltIsZero_FUNType Rvector_elt_is_zero_FUN,
		CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	int SVT_len, k, empty;
	SEXP ans, ans_elt;

	if (ndim == 1) { /* Sanity check (should never fail). */
		if (dim[0] != subarr_len)
			error("S4Arrays internal error "
			      "in REC_build_SVT_from_Rsubarray():\n"
			      "  dim[0] != subarr_len");
		return make_leaf_vector_from_Rsubvec(
					Rarray, subarr_offset, dim[0],
					posbuf,
					Rvector_elt_is_zero_FUN,
					copy_Rvector_elt_FUN);
	}

	SVT_len = dim[ndim - 1];  /* cannot be 0 so safe to divide below */
	subarr_len /= SVT_len;
	ans = PROTECT(NEW_LIST(SVT_len));
	empty = 1;
	for (k = 0; k < SVT_len; k++) {
		ans_elt = REC_build_SVT_from_Rsubarray(
					Rarray, subarr_offset, subarr_len,
					dim, ndim - 1,
					posbuf,
					Rvector_elt_is_zero_FUN,
					copy_Rvector_elt_FUN);
		if (!isNull(ans_elt)) {
			PROTECT(ans_elt);
			SET_VECTOR_ELT(ans, k, ans_elt);
			UNPROTECT(1);
			empty = 0;
		}
		subarr_offset += subarr_len;
	}
	UNPROTECT(1);
	return empty ? R_NilValue : ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_build_SVT_from_Rarray(SEXP x, SEXP as_integer)
{
	int as_int, x_ndim, *posbuf;
	RVectorEltIsZero_FUNType Rvector_elt_is_zero_FUN;
	CopyRVectorElt_FUNType copy_Rvector_elt_FUN;
	R_xlen_t x_len;
	SEXP x_dim;

	as_int = LOGICAL(as_integer)[0];
	if (as_int)
		error("'as.integer=TRUE' is not supported yet");

	Rvector_elt_is_zero_FUN = select_Rvector_elt_is_zero_FUN(TYPEOF(x));
	if (Rvector_elt_is_zero_FUN == NULL)
		error("input array has invalid type");

	copy_Rvector_elt_FUN = select_copy_Rvector_elt_FUN(TYPEOF(x));

	x_len = XLENGTH(x);
	if (x_len == 0)  /* means that 'any(dim(x) == 0)' is TRUE */
		return R_NilValue;

	x_dim = GET_DIM(x);  /* does not contain zeros */
	x_ndim = LENGTH(x_dim);
	posbuf = (int *) R_alloc(INTEGER(x_dim)[0], sizeof(int));
	return REC_build_SVT_from_Rsubarray(x, 0, x_len,
					    INTEGER(x_dim), x_ndim,
					    posbuf,
					    Rvector_elt_is_zero_FUN,
					    copy_Rvector_elt_FUN);
}


/****************************************************************************
 * Going from SVT_SparseArray to [d|l]gCMatrix
 */

static int dump_SVT_to_CsparseMatrix_slots(SEXP x_SVT, int x_ncol,
		SEXP ans_p, SEXP ans_i, SEXP ans_x)
{
	int offset, j, lv_len, ret, k;
	SEXP subSVT, lv_pos, lv_vals;

	INTEGER(ans_p)[0] = 0;
	offset = 0;
	for (j = 0; j < x_ncol; j++) {
		subSVT = VECTOR_ELT(x_SVT, j);
		if (!isNull(subSVT)) {
			/* 'subSVT' is a "leaf vector". */
			lv_len = split_leaf_vector(subSVT, &lv_pos, &lv_vals);
			if (lv_len < 0)
				return -1;
			ret = copy_Rvector_elts(lv_vals, (R_xlen_t) 0,
						ans_x, (R_xlen_t) offset,
						XLENGTH(lv_vals));
			if (ret < 0)
				return -1;
			/* Row indices in 'ans_i' must be 0-based. */
			for (k = 0; k < lv_len; k++) {
				INTEGER(ans_i)[offset] = INTEGER(lv_pos)[k] - 1;
				offset++;
			}
		}
		INTEGER(ans_p)[j + 1] = offset;
	}
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_from_SVT_SparseArray_to_CsparseMatrix(SEXP x_dim,
		SEXP x_type, SEXP x_SVT)
{
	R_xlen_t nzdata_len;
	SEXPTYPE Rtype;
	int x_ncol, ret;
	SEXP ans_p, ans_i, ans_x, ans;

	if (LENGTH(x_dim) != 2)
		error("object to coerce to dgCMatrix "
		      "must have exactly 2 dimensions");

	nzdata_len = REC_sum_leaf_vector_lengths(x_SVT, 2);
	if (nzdata_len > INT_MAX)
		error("SVT_SparseArray object contains too many nonzero "
		      "values to be turned into a COO_SparseArray object");

	Rtype = get_Rtype_from_Rstring(x_type);
	if (Rtype == 0)
		error("S4Arrays internal error in "
		      "C_from_SVT_SparseArray_to_CsparseMatrix():\n"
		      "  SVT_SparseArray object has invalid type");
	x_ncol = INTEGER(x_dim)[1];

	ans_i = PROTECT(NEW_INTEGER(nzdata_len));
	ans_x = PROTECT(allocVector(Rtype, nzdata_len));
	if (nzdata_len == 0) {
		ans_p = PROTECT(new_Rvector(INTSXP, (R_xlen_t) x_ncol + 1));
	} else {
		ans_p = PROTECT(NEW_INTEGER(x_ncol + 1));
		ret = dump_SVT_to_CsparseMatrix_slots(x_SVT, x_ncol,
							 ans_p, ans_i, ans_x);
		if (ret < 0) {
			UNPROTECT(3);
			error("S4Arrays internal error "
			      "in C_from_SVT_SparseArray_to_CsparseMatrix():\n"
			      "  invalid SVT_SparseArray object");
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
 * Going from dgCMatrix to SVT_SparseArray
 */

static SEXP build_leaf_vector_from_dgCMatrix_col(SEXP x_i, SEXP x_x,
						 int offset, int lv_len)
{
	SEXP lv_pos, lv_vals, ans;
	int k;

	lv_pos  = PROTECT(NEW_INTEGER(lv_len));
	lv_vals = PROTECT(NEW_NUMERIC(lv_len));
	for (k = 0; k < lv_len; k++) {
		INTEGER(lv_pos)[k] = INTEGER(x_i)[offset] + 1;  /* 1-based */
		REAL(lv_vals)[k]   = REAL(x_x)[offset];
		offset++;
	}
	ans = new_leaf_vector(lv_pos, lv_vals);
	UNPROTECT(2);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_build_SVT_from_dgCMatrix(SEXP x, SEXP as_integer)
{
	SEXP x_Dim, x_p, x_i, x_x, ans, lv;
	int as_int, x_ncol, j, offset, lv_len;

	as_int = LOGICAL(as_integer)[0];
	if (as_int)
		error("'as.integer=TRUE' is not supported yet");

	x_Dim = GET_SLOT(x, install("Dim"));
	x_ncol = INTEGER(x_Dim)[1];
	x_p = GET_SLOT(x, install("p"));

	if (INTEGER(x_p)[x_ncol] == 0)
		return R_NilValue;

	x_i = GET_SLOT(x, install("i"));
	x_x = GET_SLOT(x, install("x"));

	ans = PROTECT(NEW_LIST(x_ncol));
	for (j = 0; j < x_ncol; j++) {
		offset = INTEGER(x_p)[j];
		lv_len = INTEGER(x_p)[j + 1] - offset;
		if (lv_len != 0) {
			lv = PROTECT(
				build_leaf_vector_from_dgCMatrix_col(x_i, x_x,
							offset, lv_len)
			);
			SET_VECTOR_ELT(ans, j, lv);
			UNPROTECT(1);
		}
	}
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * Going from SVT_SparseArray to COO_SparseArray
 */

static SEXP alloc_nzdata(R_xlen_t nzdata_len, SEXP type)
{
	SEXPTYPE Rtype;

	Rtype = get_Rtype_from_Rstring(type);
	if (Rtype == 0)
		error("S4Arrays internal error in "
		      "alloc_nzdata():\n"
		      "  SVT_SparseArray object has invalid type");
	return allocVector(Rtype, nzdata_len);
}

/* Recursive. */
static int REC_extract_nzindex_and_nzdata_from_SVT(SEXP SVT,
		SEXP nzdata, int *nzdata_offset,
		int *nzindex, int nzindex_nrow, int nzindex_ncol,
		int *rowbuf, int rowbuf_offset)
{
	int SVT_len, k, ret, lv_len, *p, j;
	SEXP subSVT, lv_pos, lv_vals;

	if (isNull(SVT))
		return 0;

	if (rowbuf_offset > 0) {
		if (!isVectorList(SVT))  // IS_LIST() is broken
			return -1;
		SVT_len = LENGTH(SVT);
		for (k = 0; k < SVT_len; k++) {
			subSVT = VECTOR_ELT(SVT, k);
			rowbuf[rowbuf_offset] = k + 1;
			ret = REC_extract_nzindex_and_nzdata_from_SVT(
					subSVT,
					nzdata, nzdata_offset,
					nzindex, nzindex_nrow, nzindex_ncol,
					rowbuf, rowbuf_offset - 1);
			if (ret < 0)
				return -1;
		}
		return 0;
	}

	/* 'SVT' is a "leaf vector". */
	lv_len = split_leaf_vector(SVT, &lv_pos, &lv_vals);
	if (lv_len < 0)
		return -1;

	ret = copy_Rvector_elts(lv_vals, (R_xlen_t) 0,
				nzdata, (R_xlen_t) *nzdata_offset,
				XLENGTH(lv_vals));
	if (ret < 0)
		return -1;

	for (k = 0; k < lv_len; k++) {
		rowbuf[0] = INTEGER(lv_pos)[k];

		/* Copy 'rowbuf' to 'nzindex'. */
		p = nzindex + *nzdata_offset;
		for (j = 0; j < nzindex_ncol; j++) {
			*p = rowbuf[j];
			p += nzindex_nrow;
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
	int nzindex_nrow, nzindex_ncol, *rowbuf, nzdata_offset, ret;
	SEXP nzindex, nzdata, ans;

	nzdata_len = REC_sum_leaf_vector_lengths(x_SVT, LENGTH(x_dim));
	if (nzdata_len > INT_MAX)
		error("SVT_SparseArray object contains too many nonzero "
		      "values to be turned into a COO_SparseArray object");

	nzdata = PROTECT(alloc_nzdata(nzdata_len, x_type));

	nzindex_nrow = (int) nzdata_len;
	nzindex_ncol = LENGTH(x_dim);
	rowbuf = (int *) R_alloc(nzindex_ncol, sizeof(int));
	nzindex = PROTECT(allocMatrix(INTSXP, nzindex_nrow, nzindex_ncol));

	nzdata_offset = 0;
	ret = REC_extract_nzindex_and_nzdata_from_SVT(x_SVT,
			nzdata, &nzdata_offset,
			INTEGER(nzindex), nzindex_nrow, nzindex_ncol,
			rowbuf, nzindex_ncol - 1);
	if (ret < 0) {
		UNPROTECT(2);
		error("S4Arrays internal error "
		      "in C_from_SVT_SparseArray_to_COO_SparseArray():\n"
		      "  invalid SVT_SparseArray object");
	}

	/* Sanity check (should never fail). */
	if (nzdata_offset != nzindex_nrow) {
		UNPROTECT(2);
		error("S4Arrays internal error "
		      "in C_from_SVT_SparseArray_to_COO_SparseArray():\n"
		      "  *out_offset != nzindex_nrow");
	}

	ans = PROTECT(NEW_LIST(2));
	SET_VECTOR_ELT(ans, 0, nzindex);
	SET_VECTOR_ELT(ans, 1, nzdata);
	UNPROTECT(3);
	return ans;
}


/****************************************************************************
 * Going from COO_SparseArray to SVT_SparseArray
 */

static int grow_SVT(SEXP SVT,
		const int *dim, int ndim,
		const int *nzindex, int nzdata_len, int nzdata_offset)
{
	const int *p;
	int j, k;
	SEXP subSVT;

	p = nzindex + nzdata_offset;
	if (*p < 1  || *p > dim[0])
		return -1;

	if (ndim >= 3) {
		p += (size_t) nzdata_len * ndim;
		for (j = ndim - 2; j >= 1; j--) {
			p -= nzdata_len;
			k = *p - 1;
			if (k < 0 || k >= LENGTH(SVT))
				return -1;
			subSVT = VECTOR_ELT(SVT, k);
			if (j == 1)
				break;
			/* 'subSVT' is NULL or a list. */
			if (isNull(subSVT)) {
				subSVT = PROTECT(NEW_LIST(dim[j]));
				SET_VECTOR_ELT(SVT, k, subSVT);
				UNPROTECT(1);
			}
			SVT = subSVT;
		}
		/* 'subSVT' is NULL or an integer vector of counts. */
		if (isNull(subSVT)) {
			subSVT = PROTECT(
				new_Rvector(INTSXP, (R_xlen_t) dim[j])
			);
			SET_VECTOR_ELT(SVT, k, subSVT);
			UNPROTECT(1);
		}
		SVT = subSVT;
	}

	p = nzindex + nzdata_offset + nzdata_len;
	k = *p - 1;
	if (k < 0 || k >= LENGTH(SVT))
		return -1;
	INTEGER(SVT)[k]++;
	return 0;
}

static int store_nzpos_and_nzval_in_SVT(
		const int *nzindex, int nzdata_len, int nzindex_ncol,
		SEXP nzdata, int nzdata_offset,
		SEXP SVT,
		CopyRVectorElt_FUNType copy_Rvector_elt_FUN)
{
	const int *p;
	int j, k, ret;
	SEXP subSVT;

	if (nzindex_ncol >= 3) {
		p = nzindex + nzdata_offset +
			      (size_t) nzdata_len * nzindex_ncol;
		for (j = nzindex_ncol - 2; j >= 1; j--) {
			p -= nzdata_len;
			k = *p - 1;
			subSVT = VECTOR_ELT(SVT, k);
			if (j == 1)
				break;
			SVT = subSVT;
		}
		/* 'subSVT' is an integer vector of counts or a list. */
		if (IS_INTEGER(subSVT)) {
			subSVT = PROTECT(
				alloc_list_of_appendable_leaf_vectors(
						INTEGER(subSVT),
						LENGTH(subSVT),
						TYPEOF(nzdata))
			);
			SET_VECTOR_ELT(SVT, k, subSVT);
			UNPROTECT(1);
		}
		SVT = subSVT;
	}

	p = nzindex + nzdata_offset + nzdata_len;
	k = *p - 1;
	subSVT = VECTOR_ELT(SVT, k);

	/* 'subSVT' is an "appendable leaf vector". */
	ret = append_pos_val_pair_to_leaf_vector(subSVT,
						 nzindex[nzdata_offset],
						 nzdata, nzdata_offset,
						 copy_Rvector_elt_FUN);
	if (ret < 0)
		return -1;
	if (ret == 1) {
		/* Appendable leaf vector 'subSVT' is now full.
		   Replace it with a regular (i.e. non-appendable) "leaf
		   vector". */
		subSVT = PROTECT(
			new_leaf_vector(VECTOR_ELT(subSVT, 0),
					VECTOR_ELT(subSVT, 1))
		);
		SET_VECTOR_ELT(SVT, k, subSVT);
		UNPROTECT(1);
	}
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_build_SVT_from_COO_SparseArray(SEXP x_dim, SEXP x_nzindex, SEXP x_nzdata,
		SEXP as_integer)
{
	int as_int, x_ndim, nzdata_len, ans_len, i, ret;
	CopyRVectorElt_FUNType copy_Rvector_elt_FUN;
	SEXP x_nzindex_dim, ans;

	as_int = LOGICAL(as_integer)[0];
	if (as_int)
		error("'as.integer=TRUE' is not supported yet");

	copy_Rvector_elt_FUN = select_copy_Rvector_elt_FUN(TYPEOF(x_nzdata));
	if (copy_Rvector_elt_FUN == NULL)
		error("'x@nzdata' has invalid type");

	x_ndim = LENGTH(x_dim);
	nzdata_len = LENGTH(x_nzdata);

	/* Check 'x_nzindex' dimensions. */
	x_nzindex_dim = GET_DIM(x_nzindex);
	if (LENGTH(x_nzindex_dim) != 2)
		error("'x@nzindex' must be a matrix");
	if (INTEGER(x_nzindex_dim)[0] != nzdata_len)
		error("nrow(x@nzindex) != length(x@nzdata)");
	if (INTEGER(x_nzindex_dim)[1] != x_ndim)
		error("ncol(x@nzindex) != length(x@dim)");

	if (nzdata_len == 0)
		return R_NilValue;

	if (x_ndim == 1)
		return make_leaf_vector(INTEGER(x_nzindex), x_nzdata,
					INTEGER(x_dim)[0]);

	ans_len = INTEGER(x_dim)[x_ndim - 1];

	/* 1st pass: Grow the branches of the tree but don't add any
	   leaf vectors to it, only compute their lengths. */
	if (x_ndim == 2) {
		ans = PROTECT(new_Rvector(INTSXP, (R_xlen_t) ans_len));
	} else {
		ans = PROTECT(NEW_LIST(ans_len));
	}
	for (i = 0; i < nzdata_len; i++) {
		ret = grow_SVT(ans,
			       INTEGER(x_dim), x_ndim,
			       INTEGER(x_nzindex), nzdata_len, i);
		if (ret < 0) {
			UNPROTECT(1);
			error("the supplied matrix contains "
			      "out-of-bound values");
		}
	}

	/* 2nd pass: Add the leaf vectors to the tree. */
	if (x_ndim == 2)
		ans = PROTECT(
			alloc_list_of_appendable_leaf_vectors(
					INTEGER(ans), ans_len,
					TYPEOF(x_nzdata))
		);
	for (i = 0; i < nzdata_len; i++) {
		ret = store_nzpos_and_nzval_in_SVT(
				INTEGER(x_nzindex), nzdata_len, x_ndim,
				x_nzdata, i,
				ans,
				copy_Rvector_elt_FUN);
		if (ret < 0) {
			UNPROTECT(1);
			error("S4Arrays internal error in "
			      "C_from_COO_SparseArray_to_SVT():\n"
			      "  store_nzpos_and_nzval_in_SVT() "
			      "returned an error");
		}
	}

	if (x_ndim == 2)
		UNPROTECT(1);
	UNPROTECT(1);
	return ans;
}

