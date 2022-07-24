#ifndef _READ_SPARSE_CSV_H_
#define _READ_SPARSE_CSV_H_

#include <Rdefines.h>

SEXP C_readSparseCSV_as_SVT_SparseMatrix(
	SEXP filexp,
	SEXP sep,
	SEXP transpose,
	SEXP csv_ncol,
	SEXP tmpenv
);

SEXP C_readSparseCSV_as_COO_SparseMatrix(
	SEXP filexp,
	SEXP sep
);

#endif  /* _READ_SPARSE_CSV_H_ */

