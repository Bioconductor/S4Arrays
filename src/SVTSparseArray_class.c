/****************************************************************************
 *             Low-level manipulation of SVTSparseArray objects             *
 ****************************************************************************/
#include "SVTSparseArray_class.h"

#include <limits.h>  /* for INT_MAX */
#include <string.h>  /* for memcpy() */


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

/* General purpose copy function.
   We only support the 7 SEXP types listed in 'supported_SVT_Rtypes' above.
   Also we don't need to handle long vectors so 'in_offset', 'out_offset',
   and 'nelt' only need to be of type int. */
static inline int copy_vector_elts(SEXP in,  int in_offset,
				   SEXP out, int out_offset,
				   int nelt)
{
	SEXPTYPE Rtype;
	void *dest, *src;
	size_t n;
	int k;

	Rtype = TYPEOF(in);
	if (TYPEOF(out) != Rtype)
		return -1;
	if (in_offset  + nelt > LENGTH(in))
		return -1;
	if (out_offset + nelt > LENGTH(out))
		return -1;

	switch (Rtype) {
	    case VECSXP:
		for (k = 0; k < nelt; k++) {
			SET_VECTOR_ELT(out, out_offset + k,
			    VECTOR_ELT(in,  in_offset  + k));
		}
		return 0;

	    case STRSXP:
		for (k = 0; k < nelt; k++) {
			SET_STRING_ELT(out, out_offset + k,
			    STRING_ELT(in,  in_offset  + k));
		}
		return 0;

	    case LGLSXP: case INTSXP:
		dest = INTEGER(out) + out_offset;
		src  = INTEGER(in) + in_offset;
		n    = sizeof(int) * nelt;
		break;

	    case REALSXP:
		dest = REAL(out) + out_offset;
		src  = REAL(in) + in_offset;
		n    = sizeof(double) * nelt;
		break;

	    case CPLXSXP:
		dest = COMPLEX(out) + out_offset;
		src  = COMPLEX(in) + in_offset;
		n    = sizeof(Rcomplex) * nelt;
		break;

	    case RAWSXP:
		dest = RAW(out) + out_offset;
		src  = RAW(in) + in_offset;
		n    = sizeof(Rbyte) * nelt;
		break;

	    default: return -1;
	}
	memcpy(dest, src, n);
	return 0;
}

static SEXP make_constant_INTEGER(int len, int val)
{
	SEXP ans;
	int k;

	ans = PROTECT(NEW_INTEGER(len));
	for (k = 0; k < len; k++)
		INTEGER(ans)[k] = val;
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * Basic manipulation of a leaf vector
 *
 * A "leaf vector" is a sparse vector represented by a list of 2 parallel
 * vectors: an integer vector of positions and a vector (atomic or list)
 * of nonzero values.
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

static SEXP alloc_leaf_vector(int lv_len, SEXPTYPE Rtype)
{
	SEXP lv_pos, lv_vals, ans;

	lv_pos  = PROTECT(make_constant_INTEGER(lv_len, -1));
	lv_vals = PROTECT(allocVector(Rtype, lv_len));
	ans = new_leaf_vector(lv_pos, lv_vals);
	UNPROTECT(2);
	return ans;
}

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

static inline int append_pos_val_pair_to_leaf_vector(SEXP lv,
		int pos, SEXP nzdata, int nzdata_offset)
{
	int lv_len, k;
	SEXP lv_pos, lv_vals;

	lv_len = split_leaf_vector(lv, &lv_pos, &lv_vals);
	for (k = 0; k < lv_len; k++) {
		if (INTEGER(lv_pos)[k] == -1)
			goto ok;
	}
	return -1;
    ok:
	INTEGER(lv_pos)[k] = pos;
	/* Using copy_vector_elts() to copy a single element is not efficient.
	   TODO: Find something more efficient. */
	return copy_vector_elts(nzdata, nzdata_offset, lv_vals, k, 1);
}


/****************************************************************************
 * C_get_SVTSparseArray_nzdata_length()
 */

/* Recursive. */
static R_xlen_t sum_leaf_vector_lengths_REC(SEXP svtree, int ndim)
{
	R_xlen_t ans;
	int svtree_len, k;
	SEXP sub_svtree;

	if (isNull(svtree))
		return 0;

	if (ndim == 1) {
		/* 'svtree' is a "leaf vector". */
		return XLENGTH(VECTOR_ELT(svtree, 0));
	}

	/* 'svtree' is a regular node (list). */
	ans = 0;
	svtree_len = LENGTH(svtree);
	for (k = 0; k < svtree_len; k++) {
		sub_svtree = VECTOR_ELT(svtree, k);
		ans += sum_leaf_vector_lengths_REC(sub_svtree, ndim - 1);
	}
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_get_SVTSparseArray_nzdata_length(SEXP x_dim, SEXP x_svtree)
{
	R_xlen_t nzdata_len;

	nzdata_len = sum_leaf_vector_lengths_REC(x_svtree, LENGTH(x_dim));
	if (nzdata_len > INT_MAX)
		return ScalarReal((double) nzdata_len);
	return ScalarInteger((int) nzdata_len);
}


/****************************************************************************
 * Going from SVTSparseArray objects to COOSparseArray objects
 */

/* Also checks the supplied 'type'. */
static SEXPTYPE get_Rtype_from_SVTSparseArray_type(SEXP type)
{
	static const char *msg;
	SEXP type0;
	SEXPTYPE Rtype;
	int ntypes, i;

	msg = "S4Arrays internal error "
	      "in get_Rtype_from_SVTSparseArray_type():\n"
	      "  SVTSparseArray object has invalid type";
	if (!IS_CHARACTER(type) || LENGTH(type) != 1)
		error(msg);
	type0 = STRING_ELT(type, 0);
	if (type0 == NA_STRING)
		error(msg);
	Rtype = str2type(CHAR(type0));
	ntypes = sizeof(supported_SVT_Rtypes) / sizeof(SEXPTYPE);
	for (i = 0; i < ntypes; i++)
		if (Rtype == supported_SVT_Rtypes[i])
			return Rtype;
	error(msg);
	return 0;
}

static SEXP alloc_nzdata(R_xlen_t nzdata_len, SEXP type)
{
	SEXPTYPE Rtype;

	Rtype = get_Rtype_from_SVTSparseArray_type(type);
	return allocVector(Rtype, nzdata_len);
}

/* Recursive. */
static int extract_nzindex_and_nzdata_from_svtree_REC(SEXP svtree,
		SEXP nzdata, int *nzdata_offset,
		int *nzindex, int nzindex_nrow, int nzindex_ncol,
		int *rowbuf, int rowbuf_offset)
{
	int svtree_len, k, ret, lv_len, *p, j;
	SEXP sub_svtree, lv_pos, lv_vals;

	if (isNull(svtree))
		return 0;

	if (rowbuf_offset > 0) {
		if (!isVectorList(svtree))  // IS_LIST() is broken
			return -1;
		svtree_len = LENGTH(svtree);
		for (k = 0; k < svtree_len; k++) {
			sub_svtree = VECTOR_ELT(svtree, k);
			rowbuf[rowbuf_offset] = k + 1;
			ret = extract_nzindex_and_nzdata_from_svtree_REC(
					sub_svtree,
					nzdata, nzdata_offset,
					nzindex, nzindex_nrow, nzindex_ncol,
					rowbuf, rowbuf_offset - 1);
			if (ret < 0)
				return -1;
		}
		return 0;
	}

	/* 'svtree' is a "leaf vector". */
	lv_len = split_leaf_vector(svtree, &lv_pos, &lv_vals);
	if (lv_len < 0)
		return -1;

	ret = copy_vector_elts(lv_vals, 0,
			       nzdata, *nzdata_offset,
			       LENGTH(lv_vals));
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
SEXP C_from_SVTSparseArray_to_COOSparseArray(SEXP x_dim,
		SEXP x_type, SEXP x_svtree)
{
	R_xlen_t nzdata_len;
	int nzindex_nrow, nzindex_ncol, *rowbuf, nzdata_offset, ret;
	SEXP nzindex, nzdata, ans;

	nzdata_len = sum_leaf_vector_lengths_REC(x_svtree, LENGTH(x_dim));
	if (nzdata_len > INT_MAX)
		error("SVTSparseArray object contains too many nonzero "
		      "values to be turned into a COOSparseArray object");

	nzdata = PROTECT(alloc_nzdata(nzdata_len, x_type));

	nzindex_nrow = (int) nzdata_len;
	nzindex_ncol = LENGTH(x_dim);
	rowbuf = (int *) R_alloc(nzindex_ncol, sizeof(int));
	nzindex = PROTECT(allocMatrix(INTSXP, nzindex_nrow, nzindex_ncol));

	nzdata_offset = 0;
	ret = extract_nzindex_and_nzdata_from_svtree_REC(x_svtree,
			nzdata, &nzdata_offset,
			INTEGER(nzindex), nzindex_nrow, nzindex_ncol,
			rowbuf, nzindex_ncol - 1);
	if (ret < 0) {
		UNPROTECT(2);
		error("S4Arrays internal error "
		      "in C_from_SVTSparseArray_to_COOSparseArray():\n"
		      "  invalid SVTSparseArray object");
	}

	/* Sanity check (should never fail). */
	if (nzdata_offset != nzindex_nrow) {
		UNPROTECT(2);
		error("S4Arrays internal error "
		      "in C_from_SVTSparseArray_to_COOSparseArray():\n"
		      "  *out_offset != nzindex_nrow");
	}

	ans = PROTECT(NEW_LIST(2));
	SET_VECTOR_ELT(ans, 0, nzindex);
	SET_VECTOR_ELT(ans, 1, nzdata);
	UNPROTECT(3);
	return ans;
}


/****************************************************************************
 * Going from COOSparseArray objects to SVTSparseArray objects
 */

static int grow_svtree(SEXP svtree,
		const int *dim, int ndim,
		const int *nzindex, int nzdata_len, int nzdata_offset)
{
	const int *p;
	int j, k;
	SEXP sub_svtree;

	p = nzindex + nzdata_offset;
	if (*p < 1  || *p > dim[0])
		return -1;

	if (ndim >= 3) {
		p += (size_t) nzdata_len * ndim;
		for (j = ndim - 2; j >= 1; j--) {
			p -= nzdata_len;
			k = *p - 1;
			if (k < 0 || k >= LENGTH(svtree))
				return -1;
			sub_svtree = VECTOR_ELT(svtree, k);
			if (j == 1)
				break;
			/* 'sub_svtree' is NULL or a list. */
			if (isNull(sub_svtree)) {
				sub_svtree = PROTECT(NEW_LIST(dim[j]));
				SET_VECTOR_ELT(svtree, k, sub_svtree);
				UNPROTECT(1);
			}
			svtree = sub_svtree;
		}
		/* 'sub_svtree' is NULL or an integer vector of counts. */
		if (isNull(sub_svtree)) {
			sub_svtree = PROTECT(make_constant_INTEGER(dim[j], 0));
			SET_VECTOR_ELT(svtree, k, sub_svtree);
			UNPROTECT(1);
		}
		svtree = sub_svtree;
	}

	p = nzindex + nzdata_offset + nzdata_len;
	k = *p - 1;
	if (k < 0 || k >= LENGTH(svtree))
		return -1;
	INTEGER(svtree)[k]++;
	return 0;
}

static SEXP alloc_list_of_leaf_vectors(const int *lv_lens, int lv_lens_len,
		SEXPTYPE Rtype)
{
	SEXP ans, ans_elt;
	int k, lv_len;

	ans = PROTECT(NEW_LIST(lv_lens_len));
	for (k = 0; k < lv_lens_len; k++) {
		lv_len = lv_lens[k];
		if (lv_len != 0) {
			ans_elt = PROTECT(alloc_leaf_vector(lv_len, Rtype));
			SET_VECTOR_ELT(ans, k, ans_elt);
			UNPROTECT(1);
		}
	}
	UNPROTECT(1);
	return ans;
}

static int store_nzpos_and_nzval_in_svtree(
		const int *nzindex, int nzdata_len, int nzindex_ncol,
		SEXP nzdata, int nzdata_offset,
		SEXP svtree)
{
	const int *p;
	int j, k;
	SEXP sub_svtree;

	if (nzindex_ncol >= 3) {
		p = nzindex + nzdata_offset +
			      (size_t) nzdata_len * nzindex_ncol;
		for (j = nzindex_ncol - 2; j >= 1; j--) {
			p -= nzdata_len;
			k = *p - 1;
			sub_svtree = VECTOR_ELT(svtree, k);
			if (j == 1)
				break;
			svtree = sub_svtree;
		}
		/* 'sub_svtree' is an integer vector of counts or a list. */
		if (IS_INTEGER(sub_svtree)) {
			sub_svtree = PROTECT(
				alloc_list_of_leaf_vectors(INTEGER(sub_svtree),
							   LENGTH(sub_svtree),
							   TYPEOF(nzdata))
			);
			SET_VECTOR_ELT(svtree, k, sub_svtree);
			UNPROTECT(1);
		}
		svtree = sub_svtree;
	}

	p = nzindex + nzdata_offset + nzdata_len;
	k = *p - 1;
	sub_svtree = VECTOR_ELT(svtree, k);

	/* 'sub_svtree' is a "leaf vector". */
	return append_pos_val_pair_to_leaf_vector(sub_svtree,
						  nzindex[nzdata_offset],
						  nzdata, nzdata_offset);
}

/* --- .Call ENTRY POINT --- */
SEXP C_from_COOSparseArray_to_SVTSparseArray(SEXP x_dim,
		SEXP x_nzindex, SEXP x_nzdata)
{
	int ndim, nzdata_len, ans_len, i, ret;
	SEXP x_nzindex_dim, ans;

	ndim = LENGTH(x_dim);
	nzdata_len = LENGTH(x_nzdata);

	/* Check 'x_nzindex' dimensions. */
	x_nzindex_dim = GET_DIM(x_nzindex);
	if (LENGTH(x_nzindex_dim) != 2)
		error("'x@nzindex' must be a matrix");
	if (INTEGER(x_nzindex_dim)[0] != nzdata_len)
		error("nrow(x@nzindex) != length(x@nzdata)");
	if (INTEGER(x_nzindex_dim)[1] != ndim)
		error("ncol(x@nzindex) != length(x@dim)");

	if (nzdata_len == 0)
		return R_NilValue;

	if (ndim == 1)
		return make_leaf_vector(INTEGER(x_nzindex), x_nzdata,
					INTEGER(x_dim)[0]);

	ans_len = INTEGER(x_dim)[ndim - 1];

	/* 1st pass: Grow the branches of the tree but don't add any
	   leaf vectors to it, only compute their lengths. */
	if (ndim == 2) {
		ans = PROTECT(make_constant_INTEGER(ans_len, 0));
	} else {
		ans = PROTECT(NEW_LIST(ans_len));
	}
	for (i = 0; i < nzdata_len; i++) {
		ret = grow_svtree(ans,
				  INTEGER(x_dim), ndim,
				  INTEGER(x_nzindex), nzdata_len, i);
		if (ret < 0) {
			UNPROTECT(1);
			error("the supplied matrix contains "
			      "out-of-bound values");
		}
	}

	/* 2nd pass: Add the leaf vectors to the tree. */
	if (ndim == 2)
		ans = PROTECT(
			alloc_list_of_leaf_vectors(INTEGER(ans), ans_len,
						   TYPEOF(x_nzdata))
		);
	for (i = 0; i < nzdata_len; i++) {
		ret = store_nzpos_and_nzval_in_svtree(
					  INTEGER(x_nzindex), nzdata_len, ndim,
					  x_nzdata, i,
					  ans);
		if (ret < 0) {
			UNPROTECT(1);
			error("S4Arrays internal error "
			      "in C_from_matrix_to_SelectionTree():\n"
			      "  add_matrix_row_to_svtree() "
			      "returned an error");
		}
	}

	if (ndim == 2)
		UNPROTECT(1);
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * C_make_SVTSparseArray_from_dgCMatrix()
 */

static SEXP make_leaf_vector_from_dgCMatrix_col(SEXP x_i, SEXP x_x,
						int offset, int lv_len)
{
	SEXP lv_pos, lv_vals, ans;
	int k;

	lv_pos  = PROTECT(NEW_INTEGER(lv_len));
	lv_vals = PROTECT(NEW_NUMERIC(lv_len));
	for (k = 0; k < lv_len; k++) {
		INTEGER(lv_pos)[k]  = INTEGER(x_i)[offset] + 1;
		REAL(lv_vals)[k]    = REAL(x_x)[offset];
		offset++;
	}
	ans = new_leaf_vector(lv_pos, lv_vals);
	UNPROTECT(2);
	return ans;
}

SEXP C_make_SVTSparseArray_from_dgCMatrix(SEXP x, SEXP as_integer)
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
				make_leaf_vector_from_dgCMatrix_col(x_i, x_x,
							offset, lv_len)
			);
			SET_VECTOR_ELT(ans, j, lv);
			UNPROTECT(1);
		}
	}
	UNPROTECT(1);
	return ans;
}

