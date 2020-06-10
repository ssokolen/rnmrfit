#define R_NO_REMAP
#define STRICT_R_HEADERS
#include <Rinternals.h>

// Import C headers for rust API
#include "wrapper/api.h"

// Actual Wrappers
SEXP fit_1d_wrapper(SEXP x, SEXP y, SEXP knots, SEXP p, SEXP lb, SEXP ub,  
		    SEXP n, SEXP nl, SEXP nb, SEXP np, SEXP nk,
		    SEXP eq, SEXP iq, SEXP neq, SEXP niq){

  fit_1d(REAL(x), REAL(y), REAL(knots), REAL(p), REAL(lb), REAL(ub), 
 	 Rf_asInteger(n), Rf_asInteger(nl), Rf_asInteger(nb), Rf_asInteger(np), Rf_asInteger(nk),
	 REAL(eq), REAL(iq),
	 Rf_asInteger(neq), Rf_asInteger(niq));
  
  return p;
}

SEXP eval_1d_wrapper(SEXP x, SEXP y, SEXP knots, SEXP p,   
	 	     SEXP n, SEXP nl, SEXP nb, SEXP np, SEXP nk, SEXP basis){

  eval_1d(REAL(x), REAL(y), REAL(knots), REAL(p),  
     	  Rf_asInteger(n), Rf_asInteger(nl), Rf_asInteger(nb), Rf_asInteger(np), Rf_asInteger(nk));
  
  return y;
}

// Actual Wrappers
SEXP fit_2d_wrapper(SEXP x_direct, SEXP x_indirect, SEXP y, 
		    SEXP resonances, SEXP dimensions, SEXP knots, 
		    SEXP p, SEXP lb, SEXP ub,  
		    SEXP n, SEXP nl, SEXP nb, SEXP np, SEXP nk,
		    SEXP eq, SEXP iq, SEXP neq, SEXP niq){

  fit_2d(REAL(x_direct), REAL(x_indirect), REAL(y), 
	 INTEGER(resonances), INTEGER(dimensions), REAL(knots), 
	 REAL(p), REAL(lb), REAL(ub), 
 	 Rf_asInteger(n), Rf_asInteger(nl), Rf_asInteger(nb), Rf_asInteger(np), Rf_asInteger(nk),
	 REAL(eq), REAL(iq),
	 Rf_asInteger(neq), Rf_asInteger(niq));
  
  return p;
}

// Actual Wrappers
SEXP eval_2d_wrapper(SEXP x_direct, SEXP x_indirect, SEXP y, 
	 	     SEXP resonances, SEXP dimensions, SEXP knots, 
		     SEXP p, SEXP n, SEXP nl, SEXP nb, SEXP np, SEXP nk){

  eval_2d(REAL(x_direct), REAL(x_indirect), REAL(y), 
	  INTEGER(resonances), INTEGER(dimensions), REAL(knots), 
	  REAL(p),  Rf_asInteger(n), Rf_asInteger(nl), Rf_asInteger(nb), Rf_asInteger(np), Rf_asInteger(nk));
  
  return y;
}

// Standard R package stuff
static const R_CallMethodDef CallEntries[] = {
  {"fit_1d_wrapper", (DL_FUNC) &fit_1d_wrapper, 15},
  {"eval_1d_wrapper", (DL_FUNC) &eval_1d_wrapper, 9},
  {"fit_2d_wrapper", (DL_FUNC) &fit_2d_wrapper, 18},
  {"eval_2d_wrapper", (DL_FUNC) &eval_2d_wrapper, 12},
  {NULL, NULL, 0}
};

void R_init_rnmrfit(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
