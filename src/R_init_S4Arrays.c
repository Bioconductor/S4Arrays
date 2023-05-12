#include <R_ext/Rdynload.h>

#include "abind.h"
#include "array_selection.h"
#include "dim_tuning_utils.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* abind.c */
	CALLMETHOD_DEF(C_abind, 3),

/* array_selection.c */
	CALLMETHOD_DEF(C_Lindex2Mindex, 3),
	CALLMETHOD_DEF(C_Mindex2Lindex, 4),

/* dim_tuning_utils.c */
	CALLMETHOD_DEF(C_tune_dims, 2),
	CALLMETHOD_DEF(C_tune_dimnames, 2),

	{NULL, NULL, 0}
};

void R_init_S4Arrays(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, 0);
	return;
}

