/****************************************************************************
 *                        C_rowsum() and C_colsum()                         *
 ****************************************************************************/
#include "rowsum.h"

#include "S4Vectors_interface.h"

#include <limits.h>  /* for INT_MAX */


static void check_group(SEXP group, int x_nrow, int ngroup)
{
	if (!IS_INTEGER(group))
		error("the grouping vector must be "
		      "an integer vector or factor");
	if (LENGTH(group) != x_nrow)
		error("the grouping vector must have one element "
		      "per row in 'x' for rowsum()\n  and one element "
		      "per column in 'x' for colsum()");
	for (int i = 0; i < x_nrow; i++) {
		int g = INTEGER(group)[i];
		if (g == NA_INTEGER) {
			if (ngroup < 1)
				error("'ngroup' must be >= 1 when 'group' "
				      "contains missing values");
		} else {
			if (g < 1 || g > ngroup)
				error("all non-NA values in 'group' must "
				      "be >= 1 and <= 'ngroup'");
		}
	}
	return;
}

static void compute_rowsum_double(const double *x, int x_nrow, int x_ncol,
		const int *groups, int narm, double *out, int out_nrow)
{
	memset(out, 0, sizeof(double) * out_nrow * x_ncol);
	error("compute_rowsum_double() is not ready");
	return;
}

static void compute_rowsum_int(const int *x, int x_nrow, int x_ncol,
		const int *groups, int narm, int *out, int out_nrow)
{
	memset(out, 0, sizeof(int) * out_nrow * x_ncol);
	error("compute_rowsum_int() is not ready");
	return;
}

static void compute_colsum_double(const double *x, int x_nrow, int x_ncol,
		const int *groups, int narm, double *out, int out_ncol)
{
	memset(out, 0, sizeof(double) * x_nrow * out_ncol);
	for (int j = 0; j < x_ncol; j++) {
		int g = groups[j];
		if (g == NA_INTEGER)
			g = out_ncol;
		g--;  // from 1-base to 0-base
		double *out_p = out + g * x_nrow;
		for (int i = 0; i < x_nrow; i++, x++, out_p++) {
			/* ISNAN(): True for *both* NA and NaN.
			   See <R_ext/Arith.h> */
			if (narm && ISNAN(*x))
				continue;
			*out_p += *x;
		}
	}
	return;
}

static void compute_colsum_int(const int *x, int x_nrow, int x_ncol,
		const int *groups, int narm, int *out, int out_ncol)
{
	memset(out, 0, sizeof(int) * x_nrow * out_ncol);
	int overflow = 0;
	for (int j = 0; j < x_ncol; j++) {
		int g = groups[j];
		if (g == NA_INTEGER)
			g = out_ncol;
		g--;  // from 1-base to 0-base
		int *out_p = out + g * x_nrow;
		for (int i = 0; i < x_nrow; i++, x++, out_p++) {
			if (*out_p == NA_INTEGER)
				continue;
			if (*x == NA_INTEGER) {
				if (!narm)
					*out_p = NA_INTEGER;
				continue;
			}
			double y = (double) *out_p + *x;
			if (-INT_MAX <= y && y <= INT_MAX) {
				*out_p = (int) y;
			} else {
				overflow = 1;
				*out_p = NA_INTEGER;
			}
		}
	}
	if (overflow)
		warning("NAs produced by integer overflow");
	return;
}


/****************************************************************************
 * C_rowsum() and C_colsum()
 */

/* --- .Call ENTRY POINT --- */
SEXP C_rowsum(SEXP x, SEXP group, SEXP ngroup, SEXP na_rm)
{
	SEXP x_dim = GET_DIM(x);
	if (x_dim == R_NilValue || LENGTH(x_dim) != 2)
		error("input object must have 2 dimensions");
	int x_nrow = INTEGER(x_dim)[0];
	int x_ncol = INTEGER(x_dim)[1];
	int narm = LOGICAL(na_rm)[0];

	int ans_nrow = INTEGER(ngroup)[0];
	check_group(group, x_nrow, ans_nrow);

	reset_ovflow_flag();
	/* Only to detect a potential integer overflow. The returned value
	   is actually not needed so we ignore it. */
	safe_int_mult(ans_nrow, x_ncol);
	if (get_ovflow_flag())
		error("too many groups (matrix of sums will be too big)");

	SEXP ans;
	/* Note that base::rowsum() only supports numeric matrices i.e.
	   matrices of type() "double" or "integer", so we do the same. */
	SEXPTYPE x_Rtype = TYPEOF(x);
	if (x_Rtype == REALSXP) {
		ans = PROTECT(allocMatrix(REALSXP, ans_nrow, x_ncol));
		compute_rowsum_double(REAL(x), x_nrow, x_ncol,
				      INTEGER(group), narm,
				      REAL(ans), ans_nrow);
	} else if (x_Rtype == INTSXP) {
		ans = PROTECT(allocMatrix(INTSXP, ans_nrow, x_ncol));
		compute_rowsum_int(INTEGER(x), x_nrow, x_ncol,
				   INTEGER(group), narm,
				   INTEGER(ans), ans_nrow);
	} else {
		error("rowsum() and colsum() do not support "
		      "matrices of type \"%s\" at the moment",
		      type2char(x_Rtype));
	}

	UNPROTECT(1);
	return ans;
}

/* --- .Call ENTRY POINT --- */
SEXP C_colsum(SEXP x, SEXP group, SEXP ngroup, SEXP na_rm)
{
	SEXP x_dim = GET_DIM(x);
	if (x_dim == R_NilValue || LENGTH(x_dim) != 2)
		error("input object must have 2 dimensions");
	int x_nrow = INTEGER(x_dim)[0];
	int x_ncol = INTEGER(x_dim)[1];
	int narm = LOGICAL(na_rm)[0];

	int ans_ncol = INTEGER(ngroup)[0];
	check_group(group, x_ncol, ans_ncol);

	reset_ovflow_flag();
	/* Only to detect a potential integer overflow. The returned value
	   is actually not needed so we ignore it. */
	safe_int_mult(x_nrow, ans_ncol);
	if (get_ovflow_flag())
		error("too many groups (matrix of sums will be too big)");

	SEXP ans;
	/* Note that base::rowsum() only supports numeric matrices i.e.
	   matrices of type() "double" or "integer", so we do the same. */
	SEXPTYPE x_Rtype = TYPEOF(x);
	if (x_Rtype == REALSXP) {
		ans = PROTECT(allocMatrix(REALSXP, x_nrow, ans_ncol));
		compute_colsum_double(REAL(x), x_nrow, x_ncol,
				      INTEGER(group), narm,
				      REAL(ans), ans_ncol);
	} else if (x_Rtype == INTSXP) {
		ans = PROTECT(allocMatrix(INTSXP, x_nrow, ans_ncol));
		compute_colsum_int(INTEGER(x), x_nrow, x_ncol,
				   INTEGER(group), narm,
				   INTEGER(ans), ans_ncol);
	} else {
		error("rowsum() and colsum() do not support "
		      "matrices of type \"%s\" at the moment",
		      type2char(x_Rtype));
	}

	UNPROTECT(1);
	return ans;
}

