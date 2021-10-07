/****************************************************************************
 *                "Summary" methods for SparseArray objects                 *
 ****************************************************************************/
#include "SparseArray_Summary.h"

#include "Rvector_utils.h"


/****************************************************************************
 * C_SVT_SparseArray_Summary()
 */

/* Recursive. */
static int REC_summarize_SVT(SEXP SVT, const int *dim, int ndim,
		SummarizeInts_FUNType summarize_ints_FUN,
		SummarizeDoubles_FUNType summarize_doubles_FUN,
		void *init, int na_rm, int status)
{
	return status;
}

SEXP C_SVT_SparseArray_Summary(SEXP x_dim, SEXP x_type, SEXP x_SVT,
		SEXP op, SEXP na_rm)
{
	SEXPTYPE Rtype;
	int opcode, narm0, status;
	SummarizeInts_FUNType summarize_ints_FUN;
	SummarizeDoubles_FUNType summarize_doubles_FUN;
	double init[2];  /* 'init' will store 1 or 2 ints or doubles */

	Rtype = _get_Rtype_from_Rstring(x_type);
        if (Rtype == 0)
		error("S4Arrays internal error in "
		      "C_SVT_SparseArray_Summary():\n"
		      "    SVT_SparseArray object has invalid type");

	opcode = _get_opcode(op, Rtype);

	if (!(IS_LOGICAL(na_rm) && LENGTH(na_rm) == 1))
		error("'na.rm' must be TRUE or FALSE");
	narm0 = LOGICAL(na_rm)[0];

	_select_Summary_FUN(opcode, Rtype,
		&summarize_ints_FUN, &summarize_doubles_FUN, init);

	status = REC_summarize_SVT(x_SVT, INTEGER(x_dim), LENGTH(x_dim),
			summarize_ints_FUN,
			summarize_doubles_FUN,
			init, narm0, 0);

	return _init2SEXP(opcode, Rtype, init, status);
}


/****************************************************************************
 * C_count_SVT_SparseArray_NAs()
 */

/* Recursive. */
static R_xlen_t REC_count_SVT_NAs(SEXP SVT, int ndim)
{
	R_xlen_t na_count;
	int SVT_len, i;
	SEXP subSVT;

	if (SVT == R_NilValue)
		return 0;

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		return _count_Rvector_NAs(VECTOR_ELT(SVT, 1));
	}

	/* 'SVT' is a regular node (list). */
	na_count = 0;
	SVT_len = LENGTH(SVT);
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		na_count += REC_count_SVT_NAs(subSVT, ndim - 1);
	}
	return na_count;
}

/* --- .Call ENTRY POINT --- */
SEXP C_count_SVT_SparseArray_NAs(SEXP x_dim, SEXP x_type, SEXP x_SVT)
{
	R_xlen_t na_count;

	na_count = REC_count_SVT_NAs(x_SVT, LENGTH(x_dim));
	if (na_count > INT_MAX)
		return ScalarReal((double) na_count);
	return ScalarInteger((int) na_count);
}


/****************************************************************************
 * C_anyNA_SVT_SparseArray()
 */

/* Recursive. */
static int REC_SVT_has_any_NA(SEXP SVT, int ndim)
{
	int SVT_len, i;
	SEXP subSVT;

	if (SVT == R_NilValue)
		return 0;

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		return _Rvector_has_any_NA(VECTOR_ELT(SVT, 1));
	}

	/* 'SVT' is a regular node (list). */
	SVT_len = LENGTH(SVT);
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		if (REC_SVT_has_any_NA(subSVT, ndim - 1))
			return 1;
	}
	return 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_anyNA_SVT_SparseArray(SEXP x_dim, SEXP x_type, SEXP x_SVT)
{
	return ScalarLogical(REC_SVT_has_any_NA(x_SVT, LENGTH(x_dim)));
}

