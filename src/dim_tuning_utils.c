/****************************************************************************
 *                           Dim tuning utilities                           *
 ****************************************************************************/
#include "dim_tuning_utils.h"


/****************************************************************************
 * Dim tuning and the 'dim_tuner' argument
 *
 * Dim tuning
 * ----------
 * Dim tuning is the act of adding and/or dropping ineffective dimensions,
 * (i.e. dimensions that have an extent of 1) to/from an array-like object.
 * Note that dim tuning doesn't change the length (which is prod(dim(.)))
 * or alter the content of the object, and is always reversible (except when
 * it drops ineffective dimensions with dimnames on them).
 *
 * The 'dim_tuner' argument
 * ------------------------
 * The exact action to perform on the dimensions of the object is encoded
 * in the 'dim_tuner' argument. This is an integer vector where each value
 * represents one of three possible operations:
 *   o  0: Keep the dimension.
 *   o -1: Drop the dimension. This operation is allowed only if the
 *         dimension to drop is ineffective (i.e. has an extent of 1).
 *   o  1: Add ineffective dimension.
 * Note that the 'dim_tuner' vector can contain any number of 1's, but the
 * number of non-positive values (i.e. number of 0 and -1 values together)
 * must match the number of dimensions of the array-like object to tune.
 *
 * Additionally, 'dim_tuner' must contain at least one 0. In other words,
 * the tuning must retain at least one of the original dimensions of the
 * object.
 *
 * Reverse tuning
 * --------------
 * To revert a dim tuning, simply tune again with '-dim_tuner' (i.e. minus
 * 'dim_tuner'). More precisely, for S4Arrays:::tune_dims() this means that
 * 'dim2' will always be identical to 'dim' here:
 *
 *      tuned_dim <- S4Arrays:::tune_dims(dim, dim_tuner)
 *      dim2 <- S4Arrays:::tune_dims(tuned_dim, -dim_tuner)
 *
 * This should be TRUE for any vector of dimensions 'dim' and any 'dim_tuner'
 * vector compatible with 'dim'.
 */

#define	KEEP_DIM        0
#define	DROP_DIM       -1
#define	ADD_DIM         1

/* Return the "new" number of dimensions i.e. the number of dims that we
   will get after tuning the current vector of dimensions. Note that this
   is simply the number of non-negative values in 'ops' (i.e. number of
   0 and 1 values together). */
static int validate_dim_tuner(const int *ops, int nops,
		const int *dims, int ndim)
{
	int along1, along2, nkept, r, op;

	along1 = along2 = nkept = 0;
	for (r = 0; r < nops; r++) {
		op = ops[r];  /* ADD_DIM, KEEP_DIM, or DROP_DIM */
		if (op == ADD_DIM) {
			along2++;
			continue;
		}
		if (along1 >= ndim)
			error("S4Arrays internal error in "
			      "validate_dim_tuner():\n"
			      "    number of 0 (KEEP) or -1 (DROP) values "
			      "in 'dim_tuner' is > 'length(dim(x))'");
		if (op == KEEP_DIM) {
			along2++;
			nkept++;
			along1++;
			continue;
		}
		if (op != DROP_DIM)
			error("S4Arrays internal error in "
			      "validate_dim_tuner():\n"
			      "    'dim_tuner' can only contain 0 (KEEP), "
			      "-1 (DROP), or 1 (ADD) values");
		if (dims[along1] != 1)
			error("S4Arrays internal error in "
			      "validate_dim_tuner():\n"
			      "    'dim_tuner[%d]' (= -1) is "
			      "mapped to 'dim(x)[%d]' (= %d)\n"
			      "    which cannot be dropped",
			      r + 1, along1 + 1, dims[along1]);
		along1++;
	}
	if (along1 < ndim)
		error("S4Arrays internal error in "
		      "validate_dim_tuner():\n"
		      "    number of 0 (KEEP) or -1 (DROP) values "
		      "in 'dim_tuner' is < 'length(dim(x))'");
	if (nkept == 0)
		error("S4Arrays internal error in "
		      "validate_dim_tuner():\n"
		      "    'dim_tuner' must contain at least one 0");
	return along2;
}


/****************************************************************************
 * C_tune_dims() and C_tune_dimnames()
 */

/* 'dim_names' must be a character vector (or R_NilValue) that contains the
   names on the vector of dimensions. This is NOT the same as the dimnames! */
static SEXP tune_dims(const int *dims, SEXP dim_names,
		const int *ops, int nops, int ans_len)
{
	SEXP ans, ans_names = R_NilValue;
	int along1, along2, r, op;

	ans = PROTECT(NEW_INTEGER(ans_len));
	if (dim_names != R_NilValue)
		ans_names = PROTECT(NEW_CHARACTER(ans_len));
	along1 = along2 = 0;
	for (r = 0; r < nops; r++) {
		op = ops[r];  /* ADD_DIM, KEEP_DIM, or DROP_DIM */
		if (op == ADD_DIM) {
			INTEGER(ans)[along2] = 1;
			along2++;
			continue;
		}
		if (op == KEEP_DIM) {
			INTEGER(ans)[along2] = dims[along1];
			if (dim_names != R_NilValue)
				SET_STRING_ELT(ans_names, along2,
					       STRING_ELT(dim_names, along1));
			along2++;
		}
		along1++;
	}
	if (dim_names != R_NilValue) {
		SET_NAMES(ans, ans_names);
		UNPROTECT(1);
	}
	UNPROTECT(1);
	return ans;
}

static SEXP tune_dimnames(SEXP dimnames,
		const int *ops, int nops, int ans_len)
{
	int along1, along2, r, op;
	SEXP ans;

	ans = PROTECT(NEW_LIST(ans_len));
	along1 = along2 = 0;
	for (r = 0; r < nops; r++) {
		op = ops[r];  /* ADD_DIM, KEEP_DIM, or DROP_DIM */
		if (op == ADD_DIM) {
			along2++;
			continue;
		}
		if (op == KEEP_DIM) {
			SET_VECTOR_ELT(ans, along2,
				       VECTOR_ELT(dimnames, along1));
			along2++;
		}
		along1++;
	}
	UNPROTECT(1);
	return ans;
}

/* Return the length of the tuned 'dimnames' or 0 if the tuning retains no
   dimnames. */
static int compute_tuned_dimnames_length(SEXP dimnames,
		const int *ops, int nops)
{
	int ndim, along1, along2, any_retained, r, op;

	if (dimnames == R_NilValue)
		return 0;
	ndim = LENGTH(dimnames);
	along1 = along2 = any_retained = 0;
	for (r = 0; r < nops; r++) {
		op = ops[r];  /* ADD_DIM, KEEP_DIM, or DROP_DIM */
		if (op == ADD_DIM) {
			along2++;
			continue;
		}
		if (along1 >= ndim)
			error("S4Arrays internal error in "
			      "compute_tuned_dimnames_length():\n"
			      "    number of 0 (KEEP) or -1 (DROP) values "
			      "in 'dim_tuner' is > 'length(dim(x))'");
		if (op == KEEP_DIM) {
			if (VECTOR_ELT(dimnames, along1) != R_NilValue)
				any_retained = 1;
			along2++;
		}
		along1++;
	}
	return any_retained ? along2 : 0;
}

/* --- .Call ENTRY POINT --- */
SEXP C_tune_dims(SEXP dim, SEXP dim_tuner)
{
	int ndim, nops, ans_len;
	const int *dims, *ops;

	ndim = LENGTH(dim);
	dims = INTEGER(dim);
	nops = LENGTH(dim_tuner);
	ops = INTEGER(dim_tuner);
	ans_len = validate_dim_tuner(ops, nops, dims, ndim);
	return tune_dims(dims, GET_NAMES(dim), ops, nops, ans_len);
}

/* --- .Call ENTRY POINT --- */
SEXP C_tune_dimnames(SEXP dimnames, SEXP dim_tuner)
{
	int nops, ans_len;
	const int *ops;

	nops = LENGTH(dim_tuner);
	ops = INTEGER(dim_tuner);
	ans_len = compute_tuned_dimnames_length(dimnames, ops, nops);
	if (ans_len == 0)
		return R_NilValue;
	return tune_dimnames(dimnames, ops, nops, ans_len);
}

