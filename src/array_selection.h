#ifndef _ARRAY_SELECTION_H_
#define _ARRAY_SELECTION_H_

#include <Rdefines.h>

SEXP C_Lindex2Mindex(SEXP Lindex, SEXP dim, SEXP use_names);
SEXP C_Mindex2Lindex(SEXP Mindex, SEXP dim, SEXP use_names, SEXP as_integer);

#endif  /* _ARRAY_SELECTION_H_ */

