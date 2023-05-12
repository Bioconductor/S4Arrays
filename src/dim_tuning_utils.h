#ifndef _DIM_TUNING_UTILS_H_
#define	_DIM_TUNING_UTILS_H_

#include <Rdefines.h>

SEXP C_tune_dims(
	SEXP dim,
	SEXP dim_selector
);

SEXP C_tune_dimnames(
	SEXP dimnames,
	SEXP dim_selector
);

#endif  /* _DIM_TUNING_UTILS_H_ */

