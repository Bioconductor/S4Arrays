/****************************************************************************
 *                   Generate a random SparseArray object                   *
 ****************************************************************************/
#include "randomSparseArray.h"

#include <R_ext/Random.h>


/****************************************************************************
 * C_poissonSparseArray()
 */

/* --- .Call ENTRY POINT --- */
SEXP C_poissonSparseArray(SEXP dim, SEXP lambda)
{
	GetRNGstate();
	for (int i = 0; i < 20; i++)
		printf("i=%d unif_rand()=%2.8f\n", i + 1, unif_rand());
	PutRNGstate();
	return R_NilValue;
}

