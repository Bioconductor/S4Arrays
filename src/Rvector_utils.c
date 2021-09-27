/****************************************************************************
 *       Some low-level convenience utilities to deal with R vectors        *
 *                      Maybe should go to S4Vectors?                       *
 ****************************************************************************/
#include "Rvector_utils.h"


size_t _get_Rtype_size(SEXPTYPE Rtype)
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
SEXP _new_Rvector(SEXPTYPE Rtype, R_xlen_t len)
{
	SEXP ans;
	size_t Rtype_size;

	ans = PROTECT(allocVector(Rtype, len));
	/* allocVector() does NOT initialize the vector elements, except
	   for a list or a character vector. */
	if (Rtype != STRSXP && Rtype != VECSXP) {
		Rtype_size = _get_Rtype_size(Rtype);
		if (Rtype_size == 0) {
			UNPROTECT(1);
			error("S4Arrays internal error in _new_Rvector():\n"
			      "    type \"%s\" is not supported",
			      type2char(Rtype));
		}
		memset(DATAPTR(ans), 0, Rtype_size * XLENGTH(ans));
	}
	UNPROTECT(1);
	return ans;
}

/* Like allocArray() but with initialization of the array elements and
   addition of the dimnames. */
SEXP _new_Rarray(SEXPTYPE Rtype, SEXP dim, SEXP dimnames)
{
	SEXP ans;
	size_t Rtype_size;

	ans = PROTECT(allocArray(Rtype, dim));
	/* allocArray() is just a thin wrapper for allocVector() and the
	   latter does NOT initialize the vector elements, except for a
	   list or a character vector. */
	if (Rtype != STRSXP && Rtype != VECSXP) {
		Rtype_size = _get_Rtype_size(Rtype);
		if (Rtype_size == 0) {
			UNPROTECT(1);
			error("S4Arrays internal error in _new_Rarray():\n"
			      "    type \"%s\" is not supported",
			      type2char(Rtype));
		}
		memset(DATAPTR(ans), 0, Rtype_size * XLENGTH(ans));
	}
	SET_DIMNAMES(ans, dimnames);
	UNPROTECT(1);
	return ans;
}

CopyRVectorElt_FUNType _select_copy_Rvector_elt_FUN(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case LGLSXP: case INTSXP: return _copy_INTEGER_elt;
	    case REALSXP:             return _copy_NUMERIC_elt;
	    case CPLXSXP:             return _copy_COMPLEX_elt;
	    case RAWSXP:              return _copy_RAW_elt;
	    case STRSXP:              return _copy_CHARACTER_elt;
	    case VECSXP:              return _copy_LIST_elt;
	}
	return NULL;
}

CopyRVectorElts_FUNType _select_copy_Rvector_elts_FUN(SEXPTYPE Rtype)
{
	switch (Rtype) {
	    case LGLSXP: case INTSXP: return _copy_INTEGER_elts;
	    case REALSXP:             return _copy_NUMERIC_elts;
	    case CPLXSXP:             return _copy_COMPLEX_elts;
	    case RAWSXP:              return _copy_RAW_elts;
	    case STRSXP:              return _copy_CHARACTER_elts;
	    case VECSXP:              return _copy_LIST_elts;
	}
	return NULL;
}


/****************************************************************************
 * _copy_selected_Rsubvec_elts()
 */

void _copy_selected_ints(const int *in,
		const int *selection, int n, int *out)
{
	for (int k = 0; k < n; k++, selection++, out++)
		*out = in[*selection];
	return;
}

void _copy_selected_doubles(const double *in,
		const int *selection, int n, double *out)
{
	for (int k = 0; k < n; k++, selection++, out++)
		*out = in[*selection];
	return;
}

void _copy_selected_Rcomplexes(const Rcomplex *in,
		const int *selection, int n, Rcomplex *out)
{
	for (int k = 0; k < n; k++, selection++, out++)
		*out = in[*selection];
	return;
}

void _copy_selected_Rbytes(const Rbyte *in,
		const int *selection, int n, Rbyte *out)
{
	for (int k = 0; k < n; k++, selection++, out++)
		*out = in[*selection];
	return;
}

/* The selection is assumed to have the same length as 'out_Rvector'.
   Only for a 'selection' and 'out_Rvector' of length <= INT_MAX.
   Do NOT use on a 'selection' or 'out_Rvector' of length > INT_MAX! */
void _copy_selected_Rsubvec_elts(
	SEXP in_Rvector, R_xlen_t in_offset,
	const int *selection, SEXP out_Rvector)
{
	SEXPTYPE Rtype;
	int out_len;
	CopyRVectorElt_FUNType copy_Rvector_elt_FUN;

	Rtype = TYPEOF(in_Rvector);
	out_len = LENGTH(out_Rvector);  /* also the length of the selection */

	/* Optimized for LGLSXP, INTSXP, REALSXP, CPLXSXP, and RAWSXP. */
	switch (TYPEOF(in_Rvector)) {
	    case LGLSXP: case INTSXP:
		_copy_selected_ints(INTEGER(in_Rvector) + in_offset,
				selection, out_len, INTEGER(out_Rvector));
		return;
	    case REALSXP:
		_copy_selected_doubles(REAL(in_Rvector) + in_offset,
				selection, out_len, REAL(out_Rvector));
		return;
	    case CPLXSXP:
		_copy_selected_Rcomplexes(COMPLEX(in_Rvector) + in_offset,
				selection, out_len, COMPLEX(out_Rvector));
		return;
	    case RAWSXP:
		_copy_selected_Rbytes(RAW(in_Rvector) + in_offset,
				selection, out_len, RAW(out_Rvector));
		return;
	}

	/* STRSXP and VECSXP cases. */
	copy_Rvector_elt_FUN = _select_copy_Rvector_elt_FUN(Rtype);
	if (copy_Rvector_elt_FUN == NULL)
		error("S4Arrays internal error in "
		      "copy_selected_Rsubvec_elts():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));

	for (R_xlen_t k = 0; k < out_len; k++, selection++) {
		R_xlen_t offset = in_offset + *selection;
		copy_Rvector_elt_FUN(in_Rvector, offset, out_Rvector, k);
	}
	return;
}


/****************************************************************************
 * _copy_Rvector_elts_from_selected_lloffsets()
 */

static void copy_ints_from_selected_lloffsets(const int *in,
		const long long *lloffsets, const int *lloffset_selection,
		int n, int *out)
{
	for (int k = 0; k < n; k++, lloffset_selection++, out++)
		*out = in[lloffsets[*lloffset_selection]];
	return;
}

static void copy_doubles_from_selected_lloffsets(const double *in,
		const long long *lloffsets, const int *lloffset_selection,
		int n, double *out)
{
	for (int k = 0; k < n; k++, lloffset_selection++, out++)
		*out = in[lloffsets[*lloffset_selection]];
	return;
}

static void copy_Rcomplexes_from_selected_lloffsets(const Rcomplex *in,
		const long long *lloffsets, const int *lloffset_selection,
		int n, Rcomplex *out)
{
	for (int k = 0; k < n; k++, lloffset_selection++, out++)
		*out = in[lloffsets[*lloffset_selection]];
	return;
}

static void copy_Rbytes_from_selected_lloffsets(const Rbyte *in,
		const long long *lloffsets, const int *lloffset_selection,
		int n, Rbyte *out)
{
	for (int k = 0; k < n; k++, lloffset_selection++, out++)
		*out = in[lloffsets[*lloffset_selection]];
	return;
}

/* The selection is assumed to have the same length as 'out_Rvector'.
   Only for an 'lloffset_selection' and 'out_Rvector' of length <= INT_MAX.
   Don't use on an 'lloffset_selection' or 'out_Rvector' of length > INT_MAX! */
void _copy_Rvector_elts_from_selected_lloffsets(SEXP in_Rvector,
		const long long *lloffsets, const int *lloffset_selection,
		SEXP out_Rvector)
{
	SEXPTYPE Rtype;
	int out_len;
	CopyRVectorElt_FUNType copy_Rvector_elt_FUN;

	Rtype = TYPEOF(in_Rvector);
	out_len = LENGTH(out_Rvector);  /* also the length of the selection */

	/* Optimized for LGLSXP, INTSXP, REALSXP, CPLXSXP, and RAWSXP. */
	switch (TYPEOF(in_Rvector)) {
	    case LGLSXP: case INTSXP:
		copy_ints_from_selected_lloffsets(INTEGER(in_Rvector),
				lloffsets, lloffset_selection, out_len,
				INTEGER(out_Rvector));
		return;
	    case REALSXP:
		copy_doubles_from_selected_lloffsets(REAL(in_Rvector),
				lloffsets, lloffset_selection, out_len,
				REAL(out_Rvector));
		return;
	    case CPLXSXP:
		copy_Rcomplexes_from_selected_lloffsets(COMPLEX(in_Rvector),
				lloffsets, lloffset_selection, out_len,
				COMPLEX(out_Rvector));
		return;
	    case RAWSXP:
		copy_Rbytes_from_selected_lloffsets(RAW(in_Rvector),
				lloffsets, lloffset_selection, out_len,
				RAW(out_Rvector));
		return;
	}

	/* STRSXP and VECSXP cases. */
	copy_Rvector_elt_FUN = _select_copy_Rvector_elt_FUN(Rtype);
	if (copy_Rvector_elt_FUN == NULL)
		error("S4Arrays internal error in "
		      "copy_Rvector_elts_from_selected_lloffsets():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));

	for (R_xlen_t k = 0; k < out_len; k++, lloffset_selection++)
		copy_Rvector_elt_FUN(in_Rvector, lloffsets[*lloffset_selection],
				     out_Rvector, k);
	return;
}

