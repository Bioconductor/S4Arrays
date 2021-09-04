#ifndef _SPARSEVECTORTREE_CLASS_H_
#define _SPARSEVECTORTREE_CLASS_H_

#include <Rdefines.h>

SEXP C_get_SparseVectorTree_nzdata_length(
	SEXP x_dim,
	SEXP x_tree
);

SEXP C_from_SparseVectorTree_to_COOSparseArray(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_tree
);

SEXP C_from_COOSparseArray_to_SparseVectorTree(
	SEXP x_dim,
	SEXP x_nzindex,
        SEXP x_nzdata
);

#endif  /* _SPARSEVECTORTREE_CLASS_H_ */

