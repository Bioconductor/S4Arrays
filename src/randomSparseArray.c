/****************************************************************************
 *                   Generate a random SparseArray object                   *
 ****************************************************************************/
#include "randomSparseArray.h"

#include <R_ext/Random.h>
#include <math.h>  /* for exp() */


/****************************************************************************
 * C_simple_rpois()
 */

#define	CUMSUM_DPOIS_MAX_LENGTH 32  /* enough to support 0 <= lambda <= 4 */

static int compute_cumsum_dpois(double *cumsum_dpois, double lambda)
{
	long double csdp, dp;
	int n;

	csdp = dp = exp(-lambda);
	if (csdp >= 1.0)  /* means 'lambda' is zero or very small */
		return 0;
	cumsum_dpois[0] = (double) csdp;
	//printf("cumsum_dpois[%d]=%2.12f\n", 0, cumsum_dpois[0]);
	for (n = 1; n < CUMSUM_DPOIS_MAX_LENGTH; n++) {
		dp *= lambda / n;
		csdp += dp;
		if ((double) csdp == cumsum_dpois[n - 1])
			return n;
		cumsum_dpois[n] = (double) csdp;
		//printf("cumsum_dpois[%d]=%2.12f\n", n, cumsum_dpois[n]);
	}
	return -1;
}

/* Returns a value >= 0 and <= 'cumsum_dpois_len'. */
static inline int bsearch_cumsum_dpois(double u,
		const double *cumsum_dpois, int cumsum_dpois_len)
{
	int k1, k2, k;
	double csdp;

	if (cumsum_dpois_len == 0)
		return 0;

	/* Compare with 'cumsum_dpois[0]'. */
	if (u < cumsum_dpois[0])
		return 0;

	/* Compare with 'cumsum_dpois[cumsum_dpois_len-1]'. */
	k2 = cumsum_dpois_len - 1;
	csdp = cumsum_dpois[k2];
	if (u >= csdp)
		return cumsum_dpois_len;

	/* Binary search.
	   Seems that using >> 1 instead of / 2 is faster, even when compiling
	   with 'gcc -O2' (one would hope that the optimizer is able to do that
	   kind of optimization). */
	k1 = 0;
	while ((k = (k1 + k2) >> 1) != k1) {
		csdp = cumsum_dpois[k];
		if (u < csdp)
			k2 = k;
		else
			k1 = k;
	}
	return k2;
}

static int simple_rpois(double lambda)
{
	static double last_lambda = -1;
	static double cumsum_dpois[CUMSUM_DPOIS_MAX_LENGTH];
	static int cumsum_dpois_len;

	double u;

	if (lambda != last_lambda) {
		cumsum_dpois_len = compute_cumsum_dpois(cumsum_dpois, lambda);
		//printf("cumsum_dpois_len = %d\n", cumsum_dpois_len);
		if (cumsum_dpois_len < 0)
			error("'lambda' too big?");
		last_lambda = lambda;
	}
/*
	u = 0.0;
	int k = bsearch_cumsum_dpois(u, cumsum_dpois, cumsum_dpois_len);
	printf("u = %2.12f --> k = %d\n", u, k);
	for (int n = 0; n < cumsum_dpois_len; n++) {
		printf("cumsum_dpois[%d] = %2.12f\n", n, cumsum_dpois[n]);
		u = cumsum_dpois[n] - 0.001;
		k = bsearch_cumsum_dpois(u, cumsum_dpois, cumsum_dpois_len);
		printf("u = %2.12f --> k = %d\n", u, k);
		u = cumsum_dpois[n];
		k = bsearch_cumsum_dpois(u, cumsum_dpois, cumsum_dpois_len);
		printf("u = %2.12f --> k = %d\n", u, k);
		u = cumsum_dpois[n] + 0.001;
		k = bsearch_cumsum_dpois(u, cumsum_dpois, cumsum_dpois_len);
		printf("u = %2.12f --> k = %d\n", u, k);
	}
*/
	u = unif_rand();
	return bsearch_cumsum_dpois(u, cumsum_dpois, cumsum_dpois_len);
}

/* --- .Call ENTRY POINT --- */
SEXP C_simple_rpois(SEXP n, SEXP lambda)
{
	int n0, i;
	double lambda0;
	SEXP ans;

	if (!IS_INTEGER(n) || LENGTH(n) != 1)
		error("'n' must be a single integer");
	n0 = INTEGER(n)[0];
	if (n0 < 0)
		error("'n' cannot be negative");

	if (!IS_NUMERIC(lambda) || LENGTH(lambda) != 1)
		error("'lambda' must be a single numeric value");
	lambda0 = REAL(lambda)[0];
	if (lambda0 < 0)
		error("'lambda' cannot be negative");

	ans = PROTECT(NEW_INTEGER(n0));
	GetRNGstate();
	for (i = 0; i < n0; i++)
		INTEGER(ans)[i] = simple_rpois(lambda0);
	PutRNGstate();
	UNPROTECT(1);
	return ans;
}


/****************************************************************************
 * C_poissonSparseArray()
 */
/* --- .Call ENTRY POINT --- */
SEXP C_poissonSparseArray(SEXP dim, SEXP lambda)
{
	return R_NilValue;
}

