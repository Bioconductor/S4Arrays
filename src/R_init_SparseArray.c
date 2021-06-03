#include <R_ext/Rdynload.h>

#include "readSparseCSV.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {

/* readSparseCSV.c */
	CALLMETHOD_DEF(C_readSparseCSV, 2),

	{NULL, NULL, 0}
};

void R_init_SparseArray(DllInfo *info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	return;
}

