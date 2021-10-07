#ifndef _SPARSEARRAY_SUMMARY_H_
#define _SPARSEARRAY_SUMMARY_H_

#include <Rdefines.h>

SEXP C_SVT_SparseArray_Summary(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP op,
	SEXP na_rm
);

SEXP C_count_SVT_SparseArray_NAs(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT
);

SEXP C_anyNA_SVT_SparseArray(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT
);

#endif  /* _SPARSEARRAY_SUMMARY_H_ */

