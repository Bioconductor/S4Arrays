#include <R_ext/Rdynload.h>

#include "abind.h"
#include "ArraySelection_class.h"
#include "array_selection.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* abind.c */
	CALLMETHOD_DEF(C_abind, 3),

/* ArraySelection_class.c */
	CALLMETHOD_DEF(C_get_SelectionTree_length, 2),
	CALLMETHOD_DEF(C_from_SelectionTree_to_matrix, 2),
	CALLMETHOD_DEF(C_from_matrix_to_SelectionTree, 2),

/* array_selection.c */
	CALLMETHOD_DEF(C_Lindex2Mindex, 3),
	CALLMETHOD_DEF(C_Mindex2Lindex, 4),

	{NULL, NULL, 0}
};

void R_init_S4Arrays(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	return;
}

