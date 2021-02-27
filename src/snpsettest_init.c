#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void qfc(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* .Call calls */
extern SEXP _snpsettest_match_cpp(SEXP xSEXP, SEXP tableSEXP);
extern SEXP _snpsettest_get_ev_from_evd(SEXP matSEXP);
extern SEXP _snpsettest_get_ev_from_svd(SEXP matSEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_snpsettest_match_cpp", (DL_FUNC) &_snpsettest_match_cpp, 2},
  {"_snpsettest_get_ev_from_evd", (DL_FUNC) &_snpsettest_get_ev_from_evd, 1},
  {"_snpsettest_get_ev_from_svd", (DL_FUNC) &_snpsettest_get_ev_from_svd, 1},
  {NULL, NULL, 0}
};

static const R_CMethodDef CEntries[] = {
  {"qfc", (DL_FUNC) &qfc, 11},
  {NULL, NULL, 0}
};

void R_init_snpsettest(DllInfo *info) {
  R_registerRoutines(info, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

