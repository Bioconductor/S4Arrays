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
 * _collect_offsets_of_nonzero_Rsubvec_elts()
 */

static int collect_offsets_of_nonzero_int_elts(
		const int *in, int in_len, int *offs_buf)
{
	int *off_p = offs_buf;
	for (int offset = 0; offset < in_len; offset++, in++)
		if (*in != 0)
			*(off_p++) = offset;
	return (int) (off_p - offs_buf);
}

static int collect_offsets_of_nonzero_double_elts(
		const double *in, int in_len, int *offs_buf)
{
	int *off_p = offs_buf;
	for (int offset = 0; offset < in_len; offset++, in++)
		if (*in != 0.0)
			*(off_p++) = offset;
	return (int) (off_p - offs_buf);
}

#define	IS_NONZERO_RCOMPLEX(x) ((x)->r != 0.0 || (x)->i != 0.0)
static int collect_offsets_of_nonzero_Rcomplex_elts(
		const Rcomplex *in, int in_len, int *offs_buf)
{
	int *off_p = offs_buf;
	for (int offset = 0; offset < in_len; offset++, in++)
		if (IS_NONZERO_RCOMPLEX(in))
			*(off_p++) = offset;
	return (int) (off_p - offs_buf);
}

static int collect_offsets_of_nonzero_Rbyte_elts(
		const Rbyte *in, int in_len, int *offs_buf)
{
	int *off_p = offs_buf;
	for (int offset = 0; offset < in_len; offset++, in++)
		if (*in != 0)
			*(off_p++) = offset;
	return (int) (off_p - offs_buf);
}

#define	IS_NONEMPTY_CHARSXP(x) ((x) == NA_STRING || XLENGTH(x) != 0)
static int collect_offsets_of_nonempty_character_elts(
		SEXP Rvector, R_xlen_t subvec_offset, int subvec_len,
		int *offs_buf)
{
	int *off_p = offs_buf;
	for (int offset = 0; offset < subvec_len; offset++, subvec_offset++) {
		SEXP Rvector_elt = STRING_ELT(Rvector, subvec_offset);
		if (IS_NONEMPTY_CHARSXP(Rvector_elt))
			*(off_p++) = offset;
	}
	return (int) (off_p - offs_buf);
}

static int collect_offsets_of_nonnull_list_elts(
		SEXP Rvector, R_xlen_t subvec_offset, int subvec_len,
		int *offs_buf)
{
	int *off_p = offs_buf;
	for (int offset = 0; offset < subvec_len; offset++, subvec_offset++) {
		SEXP Rvector_elt = VECTOR_ELT(Rvector, subvec_offset);
		if (Rvector_elt != R_NilValue)
			*(off_p++) = offset;
	}
	return (int) (off_p - offs_buf);
}

/* Only looks at the subvector of 'Rvector' made of the range of elements
   defined by 'subvec_offset' and 'subvec_len'.
   Offsets of nonzero elements are collected with respect to this subvector.
   Caller must make sure to supply an 'offs_buf' buffer that is big enough
   to store all the collected offsets. Safe choice is to make the buffer of
   length 'subvec_len'.
   Note that even though 'Rvector' can be a long vector, the subvector
   defined by 'subvec_offset/subvec_len' cannot i.e. 'subvec_len' must be
   supplied as an int.
   Returns the number of collected offsets. */
int _collect_offsets_of_nonzero_Rsubvec_elts(
		SEXP Rvector, R_xlen_t subvec_offset, int subvec_len,
		int *offs_buf)
{
	SEXPTYPE Rtype;

	Rtype = TYPEOF(Rvector);
	switch (Rtype) {
	    case LGLSXP: case INTSXP:
		return collect_offsets_of_nonzero_int_elts(
				INTEGER(Rvector) + subvec_offset,
				subvec_len, offs_buf);
	    case REALSXP:
		return collect_offsets_of_nonzero_double_elts(
				REAL(Rvector) + subvec_offset,
				subvec_len, offs_buf);
	    case CPLXSXP:
		return collect_offsets_of_nonzero_Rcomplex_elts(
				COMPLEX(Rvector) + subvec_offset,
				subvec_len, offs_buf);
	    case RAWSXP:
		return collect_offsets_of_nonzero_Rbyte_elts(
				RAW(Rvector) + subvec_offset,
				subvec_len, offs_buf);
	    case STRSXP:
		return collect_offsets_of_nonempty_character_elts(
				Rvector, subvec_offset,
				subvec_len, offs_buf);
	    case VECSXP:
		return collect_offsets_of_nonnull_list_elts(
				Rvector, subvec_offset,
				subvec_len, offs_buf);
	}
	error("S4Arrays internal error in "
	      "_collect_offsets_of_nonzero_Rsubvec_elts():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return -1;  /* will never reach this */
}


/****************************************************************************
 * _reset_selected_Rvector_elts()
 */

static void reset_selected_int_elts(int *x,
		const int *selection, int n)
{
	for (int i = 0; i < n; i++, selection++)
		x[*selection] = 0;
	return;
}

static void reset_selected_double_elts(double *x,
		const int *selection, int n)
{
	for (int i = 0; i < n; i++, selection++)
		x[*selection] = 0.0;
	return;
}

static void reset_selected_Rcomplex_elts(Rcomplex *x,
		const int *selection, int n)
{
	Rcomplex *x_elt;

	for (int i = 0; i < n; i++, selection++) {
		x_elt = x + *selection;
		x_elt->r = x_elt->i = 0.0;
	}
	return;
}

static void reset_selected_Rbyte_elts(Rbyte *x,
		const int *selection, int n)
{
	for (int i = 0; i < n; i++, selection++)
		x[*selection] = 0.0;
	return;
}

static void reset_selected_character_elts(SEXP x,
		const int *selection, int n)
{
	SEXP x0;

	x0 = PROTECT(mkChar(""));
	for (int i = 0; i < n; i++, selection++)
		SET_STRING_ELT(x, *selection, x0);
	UNPROTECT(1);
	return;
}

static void reset_selected_list_elts(SEXP x,
		const int *selection, int n)
{
	for (int i = 0; i < n; i++, selection++)
		SET_VECTOR_ELT(x, *selection, R_NilValue);
	return;
}

void _reset_selected_Rvector_elts(SEXP Rvector, const int *selection, int n)
{
	SEXPTYPE Rtype;

	Rtype = TYPEOF(Rvector);
	switch (Rtype) {
	    case LGLSXP: case INTSXP:
		return reset_selected_int_elts(INTEGER(Rvector),
				selection, n);
	    case REALSXP:
		return reset_selected_double_elts(REAL(Rvector),
				selection, n);
	    case CPLXSXP:
		return reset_selected_Rcomplex_elts(COMPLEX(Rvector),
				selection, n);
	    case RAWSXP:
		return reset_selected_Rbyte_elts(RAW(Rvector),
				selection, n);
	    case STRSXP:
		return reset_selected_character_elts(Rvector,
				selection, n);
	    case VECSXP:
		return reset_selected_list_elts(Rvector,
				selection, n);
	}
	error("S4Arrays internal error in "
	      "_reset_selected_Rvector_elts():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return;
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
	switch (Rtype) {
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
		      "_copy_selected_Rsubvec_elts():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));

	for (R_xlen_t k = 0; k < out_len; k++, selection++) {
		R_xlen_t offset = in_offset + *selection;
		copy_Rvector_elt_FUN(in_Rvector, offset, out_Rvector, k);
	}
	return;
}


/****************************************************************************
 * _copy_Rvector_elts_to_offsets()
 */

void _copy_ints_to_offsets(const int *in,
		const int *selection, int n, int *out)
{
	for (int k = 0; k < n; k++, selection++, in++)
		out[*selection] = *in;
	return;
}

void _copy_doubles_to_offsets(const double *in,
		const int *selection, int n, double *out)
{
	for (int k = 0; k < n; k++, selection++, in++)
		out[*selection] = *in;
	return;
}

void _copy_Rcomplexes_to_offsets(const Rcomplex *in,
		const int *selection, int n, Rcomplex *out)
{
	for (int k = 0; k < n; k++, selection++, in++)
		out[*selection] = *in;
	return;
}

void _copy_Rbytes_to_offsets(const Rbyte *in,
		const int *selection, int n, Rbyte *out)
{
	for (int k = 0; k < n; k++, selection++, in++)
		out[*selection] = *in;
	return;
}

/* The selection is assumed to have the same length as 'in_Rvector'.
   Only for a 'selection' and 'in_Rvector' of length <= INT_MAX.
   Do NOT use on a 'selection' or 'in_Rvector' of length > INT_MAX! */
void _copy_Rvector_elts_to_offsets(
		SEXP in_Rvector, const int *selection,
		SEXP out_Rvector, R_xlen_t out_offset)
{
	SEXPTYPE Rtype;
	int in_len;
	CopyRVectorElt_FUNType copy_Rvector_elt_FUN;

	Rtype = TYPEOF(in_Rvector);
	in_len = LENGTH(in_Rvector);  /* also the length of the selection */

	/* Optimized for LGLSXP, INTSXP, REALSXP, CPLXSXP, and RAWSXP. */
	switch (Rtype) {
	    case LGLSXP: case INTSXP:
		_copy_ints_to_offsets(INTEGER(in_Rvector),
				selection, in_len,
				INTEGER(out_Rvector) + out_offset);
		return;
	    case REALSXP:
		_copy_doubles_to_offsets(REAL(in_Rvector),
				selection, in_len,
				REAL(out_Rvector) + out_offset);
		return;
	    case CPLXSXP:
		_copy_Rcomplexes_to_offsets(COMPLEX(in_Rvector),
				selection, in_len,
				COMPLEX(out_Rvector) + out_offset);
		return;
	    case RAWSXP:
		_copy_Rbytes_to_offsets(RAW(in_Rvector),
				selection, in_len,
				RAW(out_Rvector) + out_offset);
		return;
	}

	/* STRSXP and VECSXP cases. */
	copy_Rvector_elt_FUN = _select_copy_Rvector_elt_FUN(Rtype);
	if (copy_Rvector_elt_FUN == NULL)
		error("S4Arrays internal error in "
		      "_copy_Rvector_elts_to_offsets():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));

	for (R_xlen_t k = 0; k < in_len; k++, selection++) {
		R_xlen_t offset = out_offset + *selection;
		copy_Rvector_elt_FUN(in_Rvector, k, out_Rvector, offset);
	}
	return;
}


/****************************************************************************
 * _copy_Rvector_elts_from_selected_offsets()
 * _copy_Rvector_elts_from_selected_lloffsets()
 */

static void copy_ints_from_selected_offsets(const int *in,
		const int *offsets, const int *offset_selection,
		int n, int *out)
{
	for (int k = 0; k < n; k++, offset_selection++, out++)
		*out = in[offsets[*offset_selection]];
	return;
}
static void copy_ints_from_selected_lloffsets(const int *in,
		const long long *lloffsets, const int *lloffset_selection,
		int n, int *out)
{
	for (int k = 0; k < n; k++, lloffset_selection++, out++)
		*out = in[lloffsets[*lloffset_selection]];
	return;
}

static void copy_doubles_from_selected_offsets(const double *in,
		const int *offsets, const int *offset_selection,
		int n, double *out)
{
	for (int k = 0; k < n; k++, offset_selection++, out++)
		*out = in[offsets[*offset_selection]];
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

static void copy_Rcomplexes_from_selected_offsets(const Rcomplex *in,
		const int *offsets, const int *offset_selection,
		int n, Rcomplex *out)
{
	for (int k = 0; k < n; k++, offset_selection++, out++)
		*out = in[offsets[*offset_selection]];
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

static void copy_Rbytes_from_selected_offsets(const Rbyte *in,
		const int *offsets, const int *offset_selection,
		int n, Rbyte *out)
{
	for (int k = 0; k < n; k++, offset_selection++, out++)
		*out = in[offsets[*offset_selection]];
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
   Only for an 'offset_selection' and 'out_Rvector' of length <= INT_MAX.
   Don't use on an 'offset_selection' or 'out_Rvector' of length > INT_MAX! */
void _copy_Rvector_elts_from_selected_offsets(SEXP in_Rvector,
		const int *offsets, const int *offset_selection,
		SEXP out_Rvector)
{
	SEXPTYPE Rtype;
	int out_len;
	CopyRVectorElt_FUNType copy_Rvector_elt_FUN;

	Rtype = TYPEOF(in_Rvector);
	out_len = LENGTH(out_Rvector);  /* also the length of the selection */

	/* Optimized for LGLSXP, INTSXP, REALSXP, CPLXSXP, and RAWSXP. */
	switch (Rtype) {
	    case LGLSXP: case INTSXP:
		copy_ints_from_selected_offsets(INTEGER(in_Rvector),
				offsets, offset_selection, out_len,
				INTEGER(out_Rvector));
		return;
	    case REALSXP:
		copy_doubles_from_selected_offsets(REAL(in_Rvector),
				offsets, offset_selection, out_len,
				REAL(out_Rvector));
		return;
	    case CPLXSXP:
		copy_Rcomplexes_from_selected_offsets(COMPLEX(in_Rvector),
				offsets, offset_selection, out_len,
				COMPLEX(out_Rvector));
		return;
	    case RAWSXP:
		copy_Rbytes_from_selected_offsets(RAW(in_Rvector),
				offsets, offset_selection, out_len,
				RAW(out_Rvector));
		return;
	}

	/* STRSXP and VECSXP cases. */
	copy_Rvector_elt_FUN = _select_copy_Rvector_elt_FUN(Rtype);
	if (copy_Rvector_elt_FUN == NULL)
		error("S4Arrays internal error in "
		      "_copy_Rvector_elts_from_selected_offsets():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));

	for (R_xlen_t k = 0; k < out_len; k++, offset_selection++)
		copy_Rvector_elt_FUN(
			in_Rvector, (R_xlen_t) offsets[*offset_selection],
			out_Rvector, k);
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
	switch (Rtype) {
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
		      "_copy_Rvector_elts_from_selected_lloffsets():\n"
		      "    type \"%s\" is not supported", type2char(Rtype));

	for (R_xlen_t k = 0; k < out_len; k++, lloffset_selection++)
		copy_Rvector_elt_FUN(
			in_Rvector, (R_xlen_t) lloffsets[*lloffset_selection],
			out_Rvector, k);
	return;
}


/****************************************************************************
 * _get_summarize_opcode()
 */

int _get_summarize_opcode(SEXP op, SEXPTYPE Rtype)
{
	const char *s;

	if (!IS_CHARACTER(op) || LENGTH(op) != 1)
		error("'op' must be a single string");
	op = STRING_ELT(op, 0);
	if (op == NA_STRING)
		error("'op' cannot be NA");
	s = CHAR(op);
	if (Rtype != LGLSXP && Rtype != INTSXP && Rtype != REALSXP)
		error("%s() does not support SparseArray objects "
		      "of type \"%s\"", s, type2char(Rtype));
	if (strcmp(s, "min") == 0)
		return MIN_OPCODE;
	if (strcmp(s, "max") == 0)
		return MAX_OPCODE;
	if (strcmp(s, "range") == 0)
		return RANGE_OPCODE;
	if (strcmp(s, "sum") == 0)
		return SUM_OPCODE;
	if (strcmp(s, "prod") == 0)
		return PROD_OPCODE;
	if (strcmp(s, "sum_squares") == 0)
		return SUM_SQUARES_OPCODE;
	if (Rtype == REALSXP)
		error("%s() does not support SparseArray objects "
		      "of type \"%s\"", s, type2char(Rtype));
	if (strcmp(s, "any") == 0)
		return ANY_OPCODE;
	if (strcmp(s, "all") == 0)
		return ALL_OPCODE;
	error("'op' must be one of: \"min\", \"max\", \"range\", "
	      "\"sum\", \"prod\", \"any\", \"all\",\n"
	      "                       \"sum_squares\"");
	return 0;
}


/****************************************************************************
 * Callback functions used for summarize operations
 *
 * All these callback functions return the "new status":
 *    0 = 'init' has not been set yet
 *    1 = 'init' has been set
 *    2 = 'init' has been set and we don't need to continue (break condition)
 *
 * IMPORTANT NOTE: Most of them ignore the supplied 'status', only
 * min/max/range_ints() don't. This is because 'init' should have been set
 * by select_summarize_FUN() before the callback function gets called (see
 * below in this file). So it doesn't matter whether the supplied 'status'
 * is 0 or 1, and it doesn't matter if the callback function returns a new
 * status of 0 or 1.
 * Note that the _init2SEXP() function below in this file will ignore the
 * final status anyways, except when 'opcode' is MIN/MAX/RANGE_OPCODE
 * and 'Rtype' is INTSXP.
 * So for the callback functions that ignore the supplied 'status', the only
 * thing that matters is whether the returned 'status' is 2 or not, so the
 * caller knows whether to bail out or not (break condition).
 */

#define DOUBLE_IS_NA(x) (R_IsNA(x) || R_IsNaN(x))

/* Does NOT ignore 'status'. */
static inline int min_ints(void *init, const int *x, int n,
		int na_rm, int status)
{
	int *int_init;

	int_init = (int *) init;
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm) {
				int_init[0] = *x;
				return 2;
			}
			continue;
		}
		if (status == 0 || *x < int_init[0]) {
			int_init[0] = *x;
			status = 1;
		}
	}
	return status;
}

/* Ignores 'status'. */
static inline int min_doubles(void *init, const double *x, int n,
		int na_rm, int status)
{
	double *double_init;

	double_init = (double *) init;
	for (int i = 0; i < n; i++, x++) {
		if (DOUBLE_IS_NA(*x)) {
			if (!na_rm) {
				double_init[0] = *x;
				if (R_IsNA(*x))
					return 2;
			}
			continue;
		}
		if (!R_IsNaN(double_init[0]) && *x < double_init[0])
			double_init[0] = *x;
	}
	return status;
}

/* Does NOT ignore 'status'. */
static inline int max_ints(void *init, const int *x, int n,
		int na_rm, int status)
{
	int *int_init;

	int_init = (int *) init;
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm) {
				int_init[0] = *x;
				return 2;
			}
			continue;
		}
		if (status == 0 || *x > int_init[0]) {
			int_init[0] = *x;
			status = 1;
		}
	}
	return status;
}

/* Ignores 'status'. */
static inline int max_doubles(void *init, const double *x, int n,
		int na_rm, int status)
{
	double *double_init;

	double_init = (double *) init;
	for (int i = 0; i < n; i++, x++) {
		if (DOUBLE_IS_NA(*x)) {
			if (!na_rm) {
				double_init[0] = *x;
				if (R_IsNA(*x))
					return 2;
			}
			continue;
		}
		if (!R_IsNaN(double_init[0]) && *x > double_init[0])
			double_init[0] = *x;
	}
	return status;
}

/* Does NOT ignore 'status'. */
static inline int range_ints(void *init, const int *x, int n,
		int na_rm, int status)
{
	int *int_init;

	int_init = (int *) init;
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm) {
				int_init[0] = int_init[1] = *x;
				return 2;
			}
			continue;
		}
		if (status == 0) {
			int_init[0] = int_init[1] = *x;
			status = 1;
		} else {
			if (*x < int_init[0])
				int_init[0] = *x;
			if (*x > int_init[1])
				int_init[1] = *x;
		}
	}
	return status;
}

/* Ignores 'status'. */
static inline int range_doubles(void *init, const double *x, int n,
		int na_rm, int status)
{
	double *double_init;

	double_init = (double *) init;
	for (int i = 0; i < n; i++, x++) {
		if (DOUBLE_IS_NA(*x)) {
			if (!na_rm) {
				double_init[0] = double_init[1] = *x;
				if (R_IsNA(*x))
					return 2;
			}
			continue;
		}
		if (!R_IsNaN(double_init[0])) {
			if (*x < double_init[0])
				double_init[0] = *x;
			if (*x > double_init[1])
				double_init[1] = *x;
		}
	}
	return status;
}

/* Ignores 'status'. */
static inline int sum_ints(void *init, const int *x, int n,
		int na_rm, int status)
{
	double *double_init;

	double_init = (double *) init;
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm) {
				double_init[0] = NA_REAL;
				return 2;
			}
			continue;
		}
		double_init[0] += (double) *x;
	}
	return status;
}

/* Ignores 'status'. */
static inline int sum_doubles(void *init, const double *x, int n,
		int na_rm, int status)
{
	double *double_init;

	double_init = (double *) init;
	for (int i = 0; i < n; i++, x++) {
		if (DOUBLE_IS_NA(*x)) {
			if (!na_rm) {
				double_init[0] = *x;
				if (R_IsNA(*x))
					return 2;
			}
			continue;
		}
		if (!R_IsNaN(double_init[0]))
			double_init[0] += *x;
	}
	return status;
}

/* Ignores 'status'. */
static inline int prod_ints(void *init, const int *x, int n,
		int na_rm, int status)
{
	double *double_init;

	double_init = (double *) init;
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm) {
				double_init[0] = NA_REAL;
				return 2;
			}
			continue;
		}
		double_init[0] *= (double) *x;
	}
	return status;
}

/* Ignores 'status'. */
static inline int prod_doubles(void *init, const double *x, int n,
		int na_rm, int status)
{
	double *double_init;

	double_init = (double *) init;
	for (int i = 0; i < n; i++, x++) {
		if (DOUBLE_IS_NA(*x)) {
			if (!na_rm) {
				double_init[0] = *x;
				if (R_IsNA(*x))
					return 2;
			}
			continue;
		}
		if (!R_IsNaN(double_init[0]))
			double_init[0] *= *x;
	}
	return status;
}

/* Ignores 'status'. */
static inline int any_ints(void *init, const int *x, int n,
		int na_rm, int status)
{
	int *int_init;

	int_init = (int *) init;
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm)
				int_init[0] = *x;
			continue;
		}
		if (*x != 0) {
			int_init[0] = *x;
			return 2;
		}
	}
	return status;
}

/* Ignores 'status'. */
static inline int all_ints(void *init, const int *x, int n,
		int na_rm, int status)
{
	int *int_init;

	int_init = (int *) init;
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm)
				int_init[0] = *x;
			continue;
		}
		if (*x == 0) {
			int_init[0] = *x;
			return 2;
		}
	}
	return status;
}

/* Ignores 'status'. */
static inline int sum_int_squares(void *init, const int *x, int n,
		int na_rm, int status)
{
	double *double_init, y;

	double_init = (double *) init;
	for (int i = 0; i < n; i++, x++) {
		if (*x == NA_INTEGER) {
			if (!na_rm) {
				double_init[0] = NA_REAL;
				return 2;
			}
			continue;
		}
		y = (double) *x - double_init[1];
		double_init[0] += y * y;
	}
	return status;
}

/* Ignores 'status'. */
static inline int sum_double_squares(void *init, const double *x, int n,
		int na_rm, int status)
{
	double *double_init, y;

	double_init = (double *) init;
	for (int i = 0; i < n; i++, x++) {
		if (DOUBLE_IS_NA(*x)) {
			if (!na_rm) {
				double_init[0] = *x;
				if (R_IsNA(*x))
					return 2;
			}
			continue;
		}
		if (!R_IsNaN(double_init[0])) {
			y = (double) *x - double_init[1];
			double_init[0] += y * y;
		}
	}
	return status;
}

/* One of '*summarize_ints_FUN' or '*summarize_doubles_FUN' will be set
   to NULL and the other one to a non-NULL value. */
static void select_summarize_FUN(int opcode, SEXPTYPE Rtype, double center,
		SummarizeInts_FUNType *summarize_ints_FUN,
		SummarizeDoubles_FUNType *summarize_doubles_FUN,
		void *init)
{
	int *int_init = (int *) init;
	double *double_init = (double *) init;

	*summarize_ints_FUN = NULL;
	*summarize_doubles_FUN = NULL;
	if (opcode == ANY_OPCODE) {
		*summarize_ints_FUN = any_ints;
		int_init[0] = 0;
		return;
	}
	if (opcode == ALL_OPCODE) {
		*summarize_ints_FUN = all_ints;
		int_init[0] = 1;
		return;
	}
	if (opcode == SUM_OPCODE) {
		if (Rtype == REALSXP) {
			*summarize_doubles_FUN = sum_doubles;
		} else {
			*summarize_ints_FUN = sum_ints;
		}
		double_init[0] = 0.0;
		return;
	}
	if (opcode == PROD_OPCODE) {
		if (Rtype == REALSXP) {
			*summarize_doubles_FUN = prod_doubles;
		} else {
			*summarize_ints_FUN = prod_ints;
		}
		double_init[0] = 1.0;
		return;
	}
	if (opcode == SUM_SQUARES_OPCODE) {
		if (Rtype == REALSXP) {
			*summarize_doubles_FUN = sum_double_squares;
		} else {
			*summarize_ints_FUN = sum_int_squares;
		}
		double_init[0] = 0.0;
		double_init[1] = center;
		return;
	}
	if (Rtype == REALSXP) {
		switch (opcode) {
		    case MIN_OPCODE:
			*summarize_doubles_FUN = min_doubles;
			double_init[0] = R_PosInf;
			break;
		    case MAX_OPCODE:
			*summarize_doubles_FUN = max_doubles;
			double_init[0] = R_NegInf;
			break;
		    case RANGE_OPCODE:
			*summarize_doubles_FUN = range_doubles;
			double_init[0] = R_PosInf;
			double_init[1] = R_NegInf;
			break;
		}
		return;
	}
	/* NO initial value! */
	switch (opcode) {
	    case MIN_OPCODE:   *summarize_ints_FUN = min_ints;   break;
	    case MAX_OPCODE:   *summarize_ints_FUN = max_ints;   break;
	    case RANGE_OPCODE: *summarize_ints_FUN = range_ints; break;
	}
	return;
}

SummarizeOp _init_SummarizeOp(int opcode, SEXPTYPE Rtype,
		int na_rm, double center, void *init)
{
	SummarizeOp summarize_op;
	SummarizeInts_FUNType summarize_ints_FUN;
	SummarizeDoubles_FUNType summarize_doubles_FUN;

	select_summarize_FUN(opcode, Rtype, center,
		&summarize_ints_FUN, &summarize_doubles_FUN, init);

	summarize_op.opcode = opcode;
	summarize_op.Rtype = Rtype;
	summarize_op.na_rm = na_rm;
	summarize_op.center = center;
	summarize_op.summarize_ints_FUN = summarize_ints_FUN;
	summarize_op.summarize_doubles_FUN = summarize_doubles_FUN;
	return summarize_op;
}

int _apply_summarize_op(const SummarizeOp *summarize_op,
		void *init, const void *x, int n, int status)
{
	if (summarize_op->Rtype == INTSXP) {
		status = summarize_op->summarize_ints_FUN(
					init, (const int *) x, n,
					summarize_op->na_rm, status);
	} else {
		status = summarize_op->summarize_doubles_FUN(
					init, (const double *) x, n,
					summarize_op->na_rm, status);
	}
	return status;
}

SEXP _init2SEXP(int opcode, SEXPTYPE Rtype, void *init, int status)
{
	int *int_init = (int *) init;
	double *double_init = (double *) init;
	SEXP ans;

	if (opcode == ANY_OPCODE || opcode == ALL_OPCODE)
		return ScalarLogical(int_init[0]);
	if (opcode == MIN_OPCODE || opcode == MAX_OPCODE) {
		if (Rtype == REALSXP)
			return ScalarReal(double_init[0]);
		if (status == 0) {
			return ScalarReal(opcode == MIN_OPCODE ? R_PosInf
							       : R_NegInf);
		} else {
			return ScalarInteger(int_init[0]);
		}
	}
	if (opcode == RANGE_OPCODE) {
		if (Rtype == REALSXP) {
			ans = PROTECT(NEW_NUMERIC(2));
			REAL(ans)[0] = double_init[0];
			REAL(ans)[1] = double_init[1];
		} else {
			if (status == 0) {
				ans = PROTECT(NEW_NUMERIC(2));
				REAL(ans)[0] = R_PosInf;
				REAL(ans)[1] = R_NegInf;
			} else {
				ans = PROTECT(NEW_INTEGER(2));
				INTEGER(ans)[0] = int_init[0];
				INTEGER(ans)[1] = int_init[1];
			}
		}
		UNPROTECT(1);
		return ans;
	}
	/* 'opcode' is either SUM_OPCODE, PROD_OPCODE, or SUM_SQUARES_OPCODE. */
	if (Rtype == REALSXP)
		return ScalarReal(double_init[0]);
	/* Direct comparison with NA_REAL is safe. No need to use R_IsNA(). */
	if (double_init[0] == NA_REAL)
		return ScalarInteger(NA_INTEGER);
	if (double_init[0] <= INT_MAX && double_init[0] >= -INT_MAX)
		return ScalarInteger((int) double_init[0]);
	return ScalarReal(double_init[0]);
}


/****************************************************************************
 * _count_Rvector_NAs() and _Rvector_has_any_NA()
 */

static int count_NA_int_elts(const int *x, int n)
{
	int count;

	for (int i = count = 0; i < n; i++, x++)
		if (*x == NA_INTEGER)
			count++;
	return count;
}
static int any_NA_int_elt(const int *x, int n)
{
	for (int i = 0; i < n; i++, x++)
		if (*x == NA_INTEGER)
			return 1;
	return 0;
}

static int count_NA_double_elts(const double *x, int n)
{
	int count;

	for (int i = count = 0; i < n; i++, x++)
		if (DOUBLE_IS_NA(*x))
			count++;
	return count;
}
static int any_NA_double_elt(const double *x, int n)
{
	for (int i = 0; i < n; i++, x++)
		if (DOUBLE_IS_NA(*x))
			return 1;
	return 0;
}

#define RCOMPLEX_IS_NA(x) (DOUBLE_IS_NA((x)->r) || DOUBLE_IS_NA((x)->i))
static int count_NA_Rcomplex_elts(const Rcomplex *x, int n)
{
	int count;

	for (int i = count = 0; i < n; i++, x++)
		if (RCOMPLEX_IS_NA(x))
			count++;
	return count;
}
static int any_NA_Rcomplex_elt(const Rcomplex *x, int n)
{
	for (int i = 0; i < n; i++, x++)
		if (RCOMPLEX_IS_NA(x))
			return 1;
	return 0;
}

static int count_NA_character_elts(SEXP x)
{
	int n, count;

	n = LENGTH(x);
	for (int i = count = 0; i < n; i++)
		if (STRING_ELT(x, i) == NA_STRING)
			count++;
	return count;
}
static int any_NA_character_elt(SEXP x)
{
	int n;

	n = LENGTH(x);
	for (int i = 0; i < n; i++)
		if (STRING_ELT(x, i) == NA_STRING)
			return 1;
	return 0;
}

static inline int is_single_NA(SEXP x)
{
	SEXPTYPE Rtype;

	Rtype = TYPEOF(x);
	switch (Rtype) {
	    case LGLSXP: case INTSXP:
		if (LENGTH(x) == 1 && INTEGER(x)[0] == NA_INTEGER)
			return 1;
	    break;
	    case REALSXP:
		if (LENGTH(x) == 1 && DOUBLE_IS_NA(REAL(x)[0]))
			return 1;
	    break;
	    case CPLXSXP:
		if (LENGTH(x) == 1 && RCOMPLEX_IS_NA(COMPLEX(x)))
			return 1;
	    break;
	    case STRSXP:
		if (LENGTH(x) == 1 && STRING_ELT(x, 0) == NA_STRING)
			return 1;
	    break;
	}
	return 0;
}

static int count_NA_list_elts(SEXP x)
{
	int n, count;

	n = LENGTH(x);
	for (int i = count = 0; i < n; i++)
		if (is_single_NA(VECTOR_ELT(x, i)))
			count++;
	return count;
}
static int any_NA_list_elt(SEXP x)
{
	int n;

	n = LENGTH(x);
	for (int i = 0; i < n; i++) {
		if (is_single_NA(VECTOR_ELT(x, i)))
			return 1;
	}
	return 0;
}

int _count_Rvector_NAs(SEXP Rvector)
{
	SEXPTYPE Rtype;
	int n;

	Rtype = TYPEOF(Rvector);
	n = LENGTH(Rvector);
	switch (Rtype) {
	    case LGLSXP: case INTSXP:
			  return count_NA_int_elts(INTEGER(Rvector), n);
	    case REALSXP: return count_NA_double_elts(REAL(Rvector), n);
	    case CPLXSXP: return count_NA_Rcomplex_elts(COMPLEX(Rvector), n);
	    case RAWSXP:  return 0;
	    case STRSXP:  return count_NA_character_elts(Rvector);
	    case VECSXP:  return count_NA_list_elts(Rvector);
	}
	error("S4Arrays internal error in "
	      "_count_Rvector_NAs():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return -1;  /* will never reach this */
}

int _Rvector_has_any_NA(SEXP Rvector)
{
	SEXPTYPE Rtype;
	int n;

	Rtype = TYPEOF(Rvector);
	n = LENGTH(Rvector);
	switch (Rtype) {
	    case LGLSXP: case INTSXP:
			  return any_NA_int_elt(INTEGER(Rvector), n);
	    case REALSXP: return any_NA_double_elt(REAL(Rvector), n);
	    case CPLXSXP: return any_NA_Rcomplex_elt(COMPLEX(Rvector), n);
	    case RAWSXP:  return 0;
	    case STRSXP:  return any_NA_character_elt(Rvector);
	    case VECSXP:  return any_NA_list_elt(Rvector);
	}
	error("S4Arrays internal error in "
	      "_Rvector_has_any_NA():\n"
	      "    type \"%s\" is not supported", type2char(Rtype));
	return -1;  /* will never reach this */
}

