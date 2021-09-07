#include <R_ext/Rdynload.h>

#include "ArraySelection_class.h"
#include "SVT_SparseArray_class.h"
#include "readSparseCSV.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* ArraySelection_class.c */
	CALLMETHOD_DEF(C_get_SelectionTree_length, 2),
	CALLMETHOD_DEF(C_from_SelectionTree_to_matrix, 2),
	CALLMETHOD_DEF(C_from_matrix_to_SelectionTree, 2),

/* SVT_SparseArray_class.c */
	CALLMETHOD_DEF(C_get_SVT_SparseArray_nzdata_length, 2),
	CALLMETHOD_DEF(C_from_SVT_SparseArray_to_COO_SparseArray, 3),
	CALLMETHOD_DEF(C_from_COO_SparseArray_to_SVT_SparseArray, 3),
	CALLMETHOD_DEF(C_make_SVT_SparseArray_from_dgCMatrix, 2),
	CALLMETHOD_DEF(C_from_SVT_SparseArray_to_CsparseMatrix, 3),
	CALLMETHOD_DEF(C_from_SVT_SparseArray_to_array, 4),
	CALLMETHOD_DEF(C_from_array_to_SVT_SparseArray, 1),

/* readSparseCSV.c */
	CALLMETHOD_DEF(C_readSparseCSV, 2),

	{NULL, NULL, 0}
};

void R_init_SparseArray(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	return;
}

