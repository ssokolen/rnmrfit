// Common header for all rnmrfit C++ code

#ifndef _RNMRFIT_
#define _RNMRFIT_

#include <vector>
#include <math.h>
#include <stdlib.h>
#include <Rcpp.h>

//==============================================================================
// Prototypes

void lorentz(double p, double wl, double h, 
             int i_res, int i_par, double *grad, void *data);

void gauss(double p, double wg, double h, 
           int i_res, int i_par, double *grad, void *data);

void voigt(double p, double wl, double h,  double f,
           int i_res, int i_par, double *grad, void *data);

//==============================================================================
// Data structures

typedef struct {
  std::vector< double > x;	
  std::vector< std::vector < std::complex<double> > > peak_fit;
  std::vector< std::vector < std::complex<double> > > peak_partial;
} data_lineshape;

typedef struct {
  data_lineshape lineshape;
  std::vector< double > x;	
  std::vector< std::vector<double> > y;	
  std::vector< double > y_mod;	
  std::vector< double > y_fit;	
  std::vector< double > y_dif;	
  std::vector< std::vector<double> > basis;	
  int n_par, n_peaks, n_baseline, n_phase, count;
} data_1d;

typedef struct {
  std::vector< bool > sign;	
  std::vector< int > peak_number_1;	
  std::vector< int > peak_number_2;	
  double offset;
} data_constraint;

#endif

