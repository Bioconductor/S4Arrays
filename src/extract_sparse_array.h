#ifndef _EXTRACT_SPARSE_ARRAY_H_
#define _EXTRACT_SPARSE_ARRAY_H_

#include <Rdefines.h>

SEXP C_subset_SVT_SparseArray(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP index
);

#endif  /* _EXTRACT_SPARSE_ARRAY_H_ */

