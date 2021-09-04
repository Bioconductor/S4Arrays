#include <R_ext/Rdynload.h>

#include "ArraySelection_class.h"
#include "SparseVectorTree_class.h"
#include "readSparseCSV.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* ArraySelection_class.c */
	CALLMETHOD_DEF(C_get_SelectionTree_length, 2),
	CALLMETHOD_DEF(C_from_SelectionTree_to_matrix, 2),
	CALLMETHOD_DEF(C_from_matrix_to_SelectionTree, 2),

/* SparseVectorTree_class.c */
	CALLMETHOD_DEF(C_get_SparseVectorTree_nzdata_length, 2),
	CALLMETHOD_DEF(C_from_SparseVectorTree_to_COOSparseArray, 3),
	CALLMETHOD_DEF(C_from_COOSparseArray_to_SparseVectorTree, 3),

/* readSparseCSV.c */
	CALLMETHOD_DEF(C_readSparseCSV, 2),

	{NULL, NULL, 0}
};

void R_init_SparseArray(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	return;
}

