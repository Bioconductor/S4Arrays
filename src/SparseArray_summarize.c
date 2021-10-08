/****************************************************************************
 *              Summarization methods for SparseArray objects               *
 ****************************************************************************/
#include "SparseArray_summarize.h"

#include "Rvector_utils.h"
#include "leaf_vector_utils.h"


/****************************************************************************
 * C_summarize_SVT_SparseArray()
 */


/* Recursive. */
static int REC_summarize_SVT(SEXP SVT, const int *dim, int ndim,
		SummarizeInts_FUNType summarize_ints_FUN,
		SummarizeDoubles_FUNType summarize_doubles_FUN,
		void *init, int na_rm, int status, int *has_null_leaves)
{
	int SVT_len, i;
	SEXP subSVT;

	if (SVT == R_NilValue) {
		*has_null_leaves = 1;
		return status;
	}

	if (ndim == 1) {
		/* 'SVT' is a "leaf vector". */
		return _summarize_leaf_vector(SVT, dim[0],
					summarize_ints_FUN,
					summarize_doubles_FUN,
					init, na_rm, status);
	}

	/* 'SVT' is a regular node (list). */
	SVT_len = LENGTH(SVT);
	for (i = 0; i < SVT_len; i++) {
		subSVT = VECTOR_ELT(SVT, i);
		status = REC_summarize_SVT(subSVT, dim, ndim - 1,
					summarize_ints_FUN,
					summarize_doubles_FUN,
					init, na_rm, status,
					has_null_leaves);
		if (status == 2)
			break;
	}
	return status;
}

static int summarize_SVT(SEXP SVT, const int *dim, int ndim, SEXPTYPE Rtype,
		SummarizeInts_FUNType summarize_ints_FUN,
		SummarizeDoubles_FUNType summarize_doubles_FUN,
		void *init, int na_rm)
{
	int has_null_leaves, status;

	if (SVT == R_NilValue)
		return 0;
	
	has_null_leaves = 0;
	status = REC_summarize_SVT(SVT, dim, ndim,
				   summarize_ints_FUN,
				   summarize_doubles_FUN,
				   init, na_rm, 0,
				   &has_null_leaves);
	if (status == 2 || !has_null_leaves)
		return status;

	if (Rtype == INTSXP) {
		int zero = 0;
		status = summarize_ints_FUN(init, &zero, 1, na_rm, status);
	} else {
		double zero = 0.0;
		status = summarize_doubles_FUN(init, &zero, 1, na_rm, status);
	}
	return status;
}

SEXP C_summarize_SVT_SparseArray(SEXP x_dim, SEXP x_type, SEXP x_SVT,
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
		      "C_summarize_SVT_SparseArray():\n"
		      "    SVT_SparseArray object has invalid type");

	opcode = _get_opcode(op, Rtype);

	if (!(IS_LOGICAL(na_rm) && LENGTH(na_rm) == 1))
		error("'na.rm' must be TRUE or FALSE");
	narm0 = LOGICAL(na_rm)[0];

	_select_summarize_FUN(opcode, Rtype,
			&summarize_ints_FUN,
			&summarize_doubles_FUN,
			init);

	status = summarize_SVT(x_SVT, INTEGER(x_dim), LENGTH(x_dim), Rtype,
			summarize_ints_FUN,
			summarize_doubles_FUN,
			init, narm0);

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

