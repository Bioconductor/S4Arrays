#ifndef _SPARSEARRAY_SUBASSIGNMENT_H_
#define _SPARSEARRAY_SUBASSIGNMENT_H_

#include <Rdefines.h>

SEXP C_subassign_SVT_by_Mindex(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP Mindex,
	SEXP vals
);

SEXP C_subassign_SVT_by_Lindex(
	SEXP x_dim,
	SEXP x_type,
	SEXP x_SVT,
	SEXP Lindex,
	SEXP vals
);

#endif  /* _SPARSEARRAY_SUBASSIGNMENT_H_ */

