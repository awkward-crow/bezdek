#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP fcm_standard_c(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP fcm_infinity_c(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef call_methods[] = {
    {"fcm_standard_c", (DL_FUNC) &fcm_standard_c, 5},
    {"fcm_infinity_c", (DL_FUNC) &fcm_infinity_c, 6},
    {NULL, NULL, 0}
};

void R_init_bezdek(DllInfo *dll) {
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
