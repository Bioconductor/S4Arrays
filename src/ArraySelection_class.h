#ifndef _ARRAY_SELECTION_CLASS_H_
#define _ARRAY_SELECTION_CLASS_H_

#include <Rdefines.h>

SEXP C_get_SelectionTree_length(
	SEXP x_refdim,
	SEXP x_selection
);

SEXP C_from_SelectionTree_to_matrix(
	SEXP x_refdim,
	SEXP x_selection
);

SEXP C_from_matrix_to_SelectionTree(
	SEXP m,
	SEXP refdim
);

#endif  /* _ARRAY_SELECTION_CLASS_H_ */

