#ifndef _ROWSUM_H_
#define _ROWSUM_H_

#include <Rdefines.h>

SEXP C_rowsum(SEXP x, SEXP group, SEXP ngroup, SEXP na_rm);
SEXP C_colsum(SEXP x, SEXP group, SEXP ngroup, SEXP na_rm);

#endif  /* _ROWSUM_H_ */

