/****************************************************************************
 *                       Extract a SparseArray subset                       *
 ****************************************************************************/
#include "extract_sparse_array.h"


/* --- .Call ENTRY POINT --- */
SEXP C_extract_SVT_SparseArray_subset(SEXP x_dim, SEXP x_SVT, SEXP index)
{
	SEXP ans_dim, ans_SVT, ans;

	ans_dim = x_dim;
	ans_SVT = x_SVT;

	ans = PROTECT(NEW_LIST(2));
	SET_VECTOR_ELT(ans, 0, ans_dim);
	if (ans_SVT != R_NilValue) {
		SET_VECTOR_ELT(ans, 1, ans_SVT);
		UNPROTECT(1);
	}
	UNPROTECT(2);
	return ans;
}

