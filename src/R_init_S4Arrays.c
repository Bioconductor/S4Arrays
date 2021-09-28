#include <R_ext/Rdynload.h>

#include "abind.h"
#include "sparseMatrix_utils.h"
#include "ArraySelection_class.h"
#include "array_selection.h"
#include "SVT_SparseArray_class.h"
#include "SparseArray_subassignment.h"
#include "SparseArray_subsetting.h"
#include "SparseArray_combine.h"
#include "SparseMatrix_mult.h"
#include "readSparseCSV.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* abind.c */
	CALLMETHOD_DEF(C_abind, 3),

/* sparseMatrix_utils.c */
	CALLMETHOD_DEF(C_rowsum_dgCMatrix, 4),
	CALLMETHOD_DEF(C_colMins_dgCMatrix, 2),
	CALLMETHOD_DEF(C_colMaxs_dgCMatrix, 2),
	CALLMETHOD_DEF(C_colRanges_dgCMatrix, 2),
	CALLMETHOD_DEF(C_colVars_dgCMatrix, 2),

/* ArraySelection_class.c */
	CALLMETHOD_DEF(C_get_SelectionTree_length, 2),
	CALLMETHOD_DEF(C_from_SelectionTree_to_matrix, 2),
	CALLMETHOD_DEF(C_from_matrix_to_SelectionTree, 2),

/* array_selection.c */
	CALLMETHOD_DEF(C_Lindex2Mindex, 3),
	CALLMETHOD_DEF(C_Mindex2Lindex, 4),

/* SVT_SparseArray_class.c */
	CALLMETHOD_DEF(C_get_SVT_SparseArray_nzdata_length, 2),
	CALLMETHOD_DEF(C_set_SVT_SparseArray_type, 4),
	CALLMETHOD_DEF(C_from_SVT_SparseArray_to_Rarray, 4),
	CALLMETHOD_DEF(C_build_SVT_from_Rarray, 2),
	CALLMETHOD_DEF(C_from_SVT_SparseMatrix_to_CsparseMatrix, 3),
	CALLMETHOD_DEF(C_build_SVT_from_CsparseMatrix, 2),
	CALLMETHOD_DEF(C_from_SVT_SparseArray_to_COO_SparseArray, 3),
	CALLMETHOD_DEF(C_build_SVT_from_COO_SparseArray, 4),
	CALLMETHOD_DEF(C_transpose_SVT_SparseMatrix, 3),

/* SparseArray_subassignment.c */
	CALLMETHOD_DEF(C_subassign_SVT_by_Mindex, 5),
	CALLMETHOD_DEF(C_subassign_SVT_by_Lindex, 5),

/* SparseArray_subsetting.c */
	CALLMETHOD_DEF(C_drop_SVT_SparseArray_ineffective_dims, 4),
	CALLMETHOD_DEF(C_subset_SVT_SparseArray, 4),

/* SparseArray_combine.c */
	CALLMETHOD_DEF(C_abind_SVT_SparseArray_objects, 3),

/* SparseMatrix_mult.c */
	CALLMETHOD_DEF(C_SVT_SparseMatrix_crossprod, 5),

/* readSparseCSV.c */
	CALLMETHOD_DEF(C_readSparseCSV, 2),

	{NULL, NULL, 0}
};

void R_init_SparseArray(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	return;
}

