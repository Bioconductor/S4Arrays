/****************************************************************************
 *             Low-level manipulation of SVTSparseArray objects             *
 ****************************************************************************/
#include "SVTSparseArray_class.h"

#include <limits.h>  /* for INT_MAX */
#include <string.h>  /* for memcpy() */


/* A "leaf vector" is a sparse vector represented by a list of 2 parallel
   vectors: an integer vector of positions and a vector (atomic or list)
   of nonzero values. */
static int split_leaf_vector(SEXP svtree, SEXP *nzpos, SEXP *nzvals)
{
	R_xlen_t nzpos_len;

	/* Sanity checks (should never fail). */
	if (!isVectorList(svtree))  // IS_LIST() is broken
		return -1;
	if (LENGTH(svtree) != 2)
		return -1;
	*nzpos = VECTOR_ELT(svtree, 0);
	*nzvals = VECTOR_ELT(svtree, 1);
	if (!IS_INTEGER(*nzpos))
		return -1;
	nzpos_len = XLENGTH(*nzpos);
	if (nzpos_len > INT_MAX)
		return -1;
	if (XLENGTH(*nzvals) != nzpos_len)
		return -1;
	return (int) nzpos_len;
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

/* All the atomic types + "list". */
static const SEXPTYPE supported_Rtypes[] = {
	LGLSXP,   // "logical"
	INTSXP,   // "integer"
	REALSXP,  // "double"
	CPLXSXP,  // "complex"
	STRSXP,   // "character"
	RAWSXP,   // "raw"

	VECSXP    // "list"
};

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
	ntypes = sizeof(supported_Rtypes) / sizeof(SEXPTYPE);
	for (i = 0; i < ntypes; i++)
		if (Rtype == supported_Rtypes[i])
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

static inline int copy_nzvals_to_nzdata(SEXP nzvals, SEXP nzdata,
		int nzdata_offset)
{
	SEXPTYPE Rtype;
	int nzvals_len, k;

	Rtype = TYPEOF(nzvals);
	if (TYPEOF(nzdata) != Rtype)
		return -1;
	nzvals_len = LENGTH(nzvals);
	if (nzdata_offset > LENGTH(nzdata) - nzvals_len)
		return -1;

	switch (Rtype) {
	    case LGLSXP: case INTSXP:
		memcpy(INTEGER(nzdata) + nzdata_offset,
		       INTEGER(nzvals),
		       (size_t) nzvals_len * sizeof(int));
		break;

	    case REALSXP:
		memcpy(REAL(nzdata) + nzdata_offset,
		       REAL(nzvals),
		       (size_t) nzvals_len * sizeof(double));
		break;

	    case CPLXSXP:
		memcpy(COMPLEX(nzdata) + nzdata_offset,
		       COMPLEX(nzvals),
		       (size_t) nzvals_len * sizeof(Rcomplex));
		break;

	    case STRSXP:
		for (k = 0; k < nzvals_len; k++) {
			SET_STRING_ELT(nzdata, nzdata_offset + k,
				       STRING_ELT(nzvals, k));
		}
		break;

	    case RAWSXP:
		memcpy(RAW(nzdata) + nzdata_offset,
		       RAW(nzvals),
		       (size_t) nzvals_len * sizeof(Rbyte));
		break;

	    case VECSXP:
		for (k = 0; k < nzvals_len; k++) {
			SET_VECTOR_ELT(nzdata, nzdata_offset + k,
				       VECTOR_ELT(nzvals, k));
		}
		break;
	}
	return 0;
}

/* Recursive. */
static int extract_nzindex_and_nzdata_from_svtree_REC(SEXP svtree,
		SEXP nzdata, int *nzdata_offset,
		int *nzindex, int nzindex_nrow, int nzindex_ncol,
		int *rowbuf, int rowbuf_offset)
{
	int svtree_len, k, ret, leaf_len, *p, j;
	SEXP sub_svtree, nzpos, nzvals;

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
	leaf_len = split_leaf_vector(svtree, &nzpos, &nzvals);
	if (leaf_len < 0)
		return -1;

	ret = copy_nzvals_to_nzdata(nzvals, nzdata, *nzdata_offset);
	if (ret < 0)
		return -1;

	for (k = 0; k < leaf_len; k++) {
		rowbuf[0] = INTEGER(nzpos)[k];

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

#ifdef _TURLUTUTU_

static SEXP make_depth1_svtree(const int *in, int in_len, int refdim1)
{
	SEXP ans;
	int i, v;

	ans = PROTECT(NEW_INTEGER(in_len));
	for (i = 0; i < in_len; i++) {
		v = in[i];
		if (v < 1 || v > refdim1) {
			UNPROTECT(1);
			error("the supplied matrix contains "
			      "out-of-bound values");
		}
		INTEGER(ans)[i] = v;
	}
	UNPROTECT(1);
	return ans;
}

static int grow_svtree(const int *in,
		int in_nrow, int in_ncol, int in_offset,
		SEXP svtree, const int *refdim)
{
	const int *p;
	int j, k;
	SEXP sub_svtree;

	p = in + in_offset;
	if (*p < 1  || *p > refdim[0])
		return -1;

	if (in_ncol >= 3) {
		p += in_nrow * in_ncol;
		for (j = in_ncol - 2; j >= 1; j--) {
			p -= in_nrow;
			k = *p - 1;
			if (k < 0 || k >= LENGTH(svtree))
				return -1;
			sub_svtree = VECTOR_ELT(svtree, k);
			if (j == 1)
				break;
			/* 'sub_svtree' is NULL or a list. */
			if (isNull(sub_svtree)) {
				sub_svtree = PROTECT(NEW_LIST(refdim[j]));
				SET_VECTOR_ELT(svtree, k, sub_svtree);
				UNPROTECT(1);
			}
			svtree = sub_svtree;
		}
		/* 'sub_svtree' is NULL or an integer vector of counts. */
		if (isNull(sub_svtree)) {
			sub_svtree = PROTECT(NEW_INTEGER(refdim[j]));
			SET_VECTOR_ELT(svtree, k, sub_svtree);
			UNPROTECT(1);
		}
		svtree = sub_svtree;
	}

	p = in + in_offset + in_nrow;
	k = *p - 1;
	if (k < 0 || k >= LENGTH(svtree))
		return -1;
	INTEGER(svtree)[k]++;
	return 0;
}

static SEXP make_terminal_list(SEXP x)
{
	int x_len, k, ans_elt_len;
	SEXP ans, ans_elt;

	x_len = LENGTH(x);
	ans = PROTECT(NEW_LIST(x_len));
	for (k = 0; k < x_len; k++) {
		ans_elt_len = INTEGER(x)[k];
		if (ans_elt_len != 0) {
			ans_elt = PROTECT(NEW_INTEGER(ans_elt_len));
			SET_VECTOR_ELT(ans, k, ans_elt);
			UNPROTECT(1);
		}
	}
	UNPROTECT(1);
	return ans;
}

static int add_matrix_row_to_svtree(const int *in,
		int in_nrow, int in_ncol, int in_offset,
		SEXP svtree)
{
	const int *p;
	int j, k;
	SEXP sub_svtree;

	if (in_ncol >= 3) {
		p = in + in_offset + in_nrow * in_ncol;
		for (j = in_ncol - 2; j >= 1; j--) {
			p -= in_nrow;
			k = *p - 1;
			sub_svtree = VECTOR_ELT(svtree, k);
			if (j == 1)
				break;
			svtree = sub_svtree;
		}
		/* 'sub_svtree' is an integer vector of counts or a list. */
		if (IS_INTEGER(sub_svtree)) {
			sub_svtree = PROTECT(make_terminal_list(sub_svtree));
			SET_VECTOR_ELT(svtree, k, sub_svtree);
			UNPROTECT(1);
		}
		svtree = sub_svtree;
	}

	p = in + in_offset + in_nrow;
	k = *p - 1;
	sub_svtree = VECTOR_ELT(svtree, k);
	for (k = 0; k < LENGTH(sub_svtree); k++) {
		if (INTEGER(sub_svtree)[k] == 0)
			break;
	}
	p = in + in_offset;
	INTEGER(sub_svtree)[k] = *p;
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_from_matrix_to_SelectionTree(SEXP m, SEXP refdim)
{
	SEXP m_dim, ans;
	int m_nrow, m_ncol, ans_len, i, ret;

	m_dim = GET_DIM(m);
	if (LENGTH(m_dim) != 2)
		error("'m' must be a matrix");

	m_nrow = INTEGER(m_dim)[0];
	m_ncol = INTEGER(m_dim)[1];
	if (m_ncol != LENGTH(refdim))
		error("ncol(m) != length(refdim)");

	if (m_nrow == 0)
		return R_NilValue;

	if (m_ncol == 1)
		return make_depth1_svtree(INTEGER(m), m_nrow,
					  INTEGER(refdim)[0]);

	ans_len = INTEGER(refdim)[m_ncol - 1];

	/* 1st pass: Grow the branches of the tree but don't add any
	   leaves to it, only count them. */
	if (m_ncol == 2) {
		ans = PROTECT(NEW_INTEGER(ans_len));
	} else {
		ans = PROTECT(NEW_LIST(ans_len));
	}
	for (i = 0; i < m_nrow; i++) {
		ret = grow_svtree(INTEGER(m),
				  m_nrow, m_ncol, i,
				  ans, INTEGER(refdim));
		if (ret < 0) {
			UNPROTECT(1);
			error("the supplied matrix contains "
			      "out-of-bound values");
		}
	}

	/* 2nd pass: Add the leaves to the tree. */
	if (m_ncol == 2)
		ans = PROTECT(make_terminal_list(ans));
	for (i = 0; i < m_nrow; i++) {
		ret = add_matrix_row_to_svtree(INTEGER(m),
					       m_nrow, m_ncol, i,
					       ans);
		if (ret < 0) {
			UNPROTECT(1);
			error("S4Arrays internal error "
			      "in C_from_matrix_to_SelectionTree():\n"
			      "  add_matrix_row_to_svtree() "
			      "returned an error");
		}
	}

	if (m_ncol == 2)
		UNPROTECT(1);
	UNPROTECT(1);
	return ans;
}

#endif  /* _TURLUTUTU_ */

SEXP C_from_COOSparseArray_to_SVTSparseArray(SEXP x_dim,
		SEXP x_nzindex, SEXP x_nzdata)
{
	return R_NilValue;
}

