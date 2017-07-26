#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _BayesianTools_vsemC(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_BayesianTools_vsemC", (DL_FUNC) &_BayesianTools_vsemC, 2},
    {NULL, NULL, 0}
};

void R_init_BayesianTools(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
