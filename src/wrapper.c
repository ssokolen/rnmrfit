#define R_NO_REMAP
#define STRICT_R_HEADERS
#include <R.h>
#include <Rinternals.h>

// Import C headers for rust API
#include "wrapper/api.h"

// 1D
SEXP fit_1d_wrapper(SEXP x, SEXP y, SEXP knots, SEXP p, SEXP lb, SEXP ub,  
		    SEXP n, SEXP nl, SEXP nb, SEXP np, SEXP nk,
		    SEXP eq, SEXP iq, SEXP neq, SEXP niq,
		    SEXP alg, SEXP xtr, SEXP mxt,
		    SEXP out){

  fit_1d(REAL(x), REAL(y), REAL(knots), REAL(p), REAL(lb), REAL(ub), 
 	 Rf_asInteger(n), Rf_asInteger(nl), Rf_asInteger(nb), Rf_asInteger(np), Rf_asInteger(nk),
	 REAL(eq), REAL(iq),
	 Rf_asInteger(neq), Rf_asInteger(niq),
	 Rf_asInteger(alg), Rf_asReal(xtr), Rf_asReal(mxt),
	 INTEGER(out));

  SEXP vec = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(vec, 0, p);
  SET_VECTOR_ELT(vec, 1, out);

  UNPROTECT(1);
  return vec;
}

SEXP eval_1d_wrapper(SEXP x, SEXP y, SEXP knots, SEXP p,   
	 	     SEXP n, SEXP nl, SEXP nb, SEXP np, SEXP nk){

  eval_1d(REAL(x), REAL(y), REAL(knots), REAL(p),  
     	  Rf_asInteger(n), Rf_asInteger(nl), Rf_asInteger(nb), Rf_asInteger(np), Rf_asInteger(nk));
  
  return y;
}

SEXP baseline_1d_wrapper(SEXP x, SEXP y, SEXP knots, SEXP p,   
	 	     	 SEXP n, SEXP nb, SEXP nk){

  baseline_1d(REAL(x), REAL(y), REAL(knots), REAL(p),  
     	      Rf_asInteger(n), Rf_asInteger(nb), Rf_asInteger(nk));
  
  return y;
}

SEXP phase_1d_wrapper(SEXP x, SEXP y, SEXP p, SEXP n, SEXP np){

  phase_1d(REAL(x), REAL(y), REAL(p), Rf_asInteger(n), Rf_asInteger(np));
  
  return y;
}

// 2D
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

SEXP eval_2d_wrapper(SEXP x_direct, SEXP x_indirect, SEXP y, 
	 	     SEXP resonances, SEXP dimensions, SEXP knots, 
		     SEXP p, SEXP n, SEXP nl, SEXP nb, SEXP np, SEXP nk){

  eval_2d(REAL(x_direct), REAL(x_indirect), REAL(y), 
	  INTEGER(resonances), INTEGER(dimensions), REAL(knots), 
	  REAL(p),  Rf_asInteger(n), Rf_asInteger(nl), Rf_asInteger(nb), Rf_asInteger(np), Rf_asInteger(nk));
  
  return y;
}

SEXP baseline_2d_wrapper(SEXP x_direct, SEXP x_indirect, SEXP y, SEXP knots, 
		     	 SEXP p, SEXP n, SEXP nb, SEXP nk){

  baseline_2d(REAL(x_direct), REAL(x_indirect), REAL(y), REAL(knots), 
	      REAL(p),  Rf_asInteger(n), Rf_asInteger(nb), Rf_asInteger(nk));
  
  return y;
}

SEXP phase_2d_wrapper(SEXP x_direct, SEXP x_indirect, SEXP y, 
		      SEXP p, SEXP n, SEXP np){

  phase_2d(REAL(x_direct), REAL(x_indirect), REAL(y),
  	   REAL(p),  Rf_asInteger(n), Rf_asInteger(np));
  
  return y;
}

// Standard R package stuff
static const R_CallMethodDef CallEntries[] = {
  {"fit_1d_wrapper", (DL_FUNC) &fit_1d_wrapper, 19},
  {"eval_1d_wrapper", (DL_FUNC) &eval_1d_wrapper, 9},
  {"baseline_1d_wrapper", (DL_FUNC) &baseline_1d_wrapper, 7},
  {"phase_1d_wrapper", (DL_FUNC) &phase_1d_wrapper, 5},
  {"fit_2d_wrapper", (DL_FUNC) &fit_2d_wrapper, 18},
  {"eval_2d_wrapper", (DL_FUNC) &eval_2d_wrapper, 12},
  {"baseline_2d_wrapper", (DL_FUNC) &baseline_2d_wrapper, 8},
  {"phase_2d_wrapper", (DL_FUNC) &phase_2d_wrapper, 6},
  {NULL, NULL, 0}
};

void R_init_rnmrfit(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
