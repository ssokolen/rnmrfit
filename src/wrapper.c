#define R_NO_REMAP
#define STRICT_R_HEADERS
#include <Rinternals.h>

// Import C headers for rust API
#include "rnmrfit/api.h"

// Actual Wrappers
SEXP fit_1d_wrapper(SEXP x, SEXP y, SEXP p, SEXP lb, SEXP ub,  
		    SEXP n, SEXP nl, SEXP nb, SEXP np, SEXP basis,
		    SEXP eq, SEXP iq, SEXP neq, SEXP niq){

  fit_1d(REAL(x), REAL(y), REAL(p), REAL(lb), REAL(ub), 
 	 Rf_asInteger(n), Rf_asInteger(nl), Rf_asInteger(nb), Rf_asInteger(np),
	 REAL(basis),
	 REAL(eq), REAL(iq),
	 Rf_asInteger(neq), Rf_asInteger(niq));
  
  return p;
}

SEXP eval_1d_wrapper(SEXP x, SEXP y, SEXP p,   
	 	     SEXP n, SEXP nl, SEXP nb, SEXP np, SEXP basis){

  eval_1d(REAL(x), REAL(y), REAL(p),  
     	  Rf_asInteger(n), Rf_asInteger(nl), Rf_asInteger(nb), Rf_asInteger(np),
	  REAL(basis));
  
  return y;
}

// Standard R package stuff
static const R_CallMethodDef CallEntries[] = {
  {"fit_1d_wrapper", (DL_FUNC) &fit_1d_wrapper, 14},
  {"eval_1d_wrapper", (DL_FUNC) &eval_1d_wrapper, 8},
  {NULL, NULL, 0}
};

void R_init_rnmrfit(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
