
// See Sections 5.4 and 6.16 of Writing R Extensions

#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "sped.h"

static R_CallMethodDef callMethods[]  = {
    {"descent", (DL_FUNC) &descent, 4},
    {NULL, NULL, 0}
};

void attribute_visible R_init_sped(DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}

