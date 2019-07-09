#include <vector>
#include <math.h>
#include <stdlib.h>
#include <Rcpp.h>
#include "rnmrfit.h"



//==============================================================================
// 1D lineshape function
//------------------------------------------------------------------------------



// [[Rcpp::export]]
std::vector < std::complex<double> > 
  lineshape_1d(Rcpp::NumericVector x_direct, Rcpp::NumericVector par) {

  using namespace std;

  //--------------------
  // Converting and storing data (to make use of generic fit functions)
  
  int n_points = x_direct.size();
  int n_par = par.size();
  int n_peaks = n_par/4;

  vector< double > x(n_points, 0);

  // Converting Rcpp objects to std::vector to keep things standardized
  for (int i = 0; i < n_points; i++) {
    x.at(i) = x_direct.at(i);
  }

  data_lineshape data_direct[1];
  
  data_direct[0].x = x;

  data_direct[0].peak_fit = 
    vector< vector < complex<double> > >
      ( 1, vector< complex<double> > 
        ( n_points, complex<double> (0,0) ) );

  //-------------------- 
  // Loop over peaks
  
  for (int i = 0; i < n_peaks; i++) {

    double p = par[i*4];
    double w = par[i*4+1];
    double h = par[i*4+2];
    double f = par[i*4+3];

    // If the fraction gauss value is 0, proceed as Lorentz (w becomes wl)
    if ( f < 1e-6 ) {
		  lorentz(p, w, h, 0, i*4, NULL, &data_direct[0]);	
    } else if ( f < (1 - 1e-6) ) {
		  voigt(p, w, h, f, 0, i*4, NULL, &data_direct[0]);	
    } else if ( f < (1 + 1e-6) ) {
		  gauss(p, w, h, 0, i*4, NULL, &data_direct[0]);	
    } else {
      Rcpp::stop("Fraction gauss outside valid [0,1] range. Aborting");
    }
  }
  
  //--------------------
  // For 1D peaks, the data_direct.peak_fit can be read directly into y


  vector< complex<double> > y_val(n_points, complex<double>(0, 0));
  vector< complex<double> > peak_fit = data_direct[0].peak_fit.at(0);

  for (int i = 0; i < n_points; i++) {
    y_val.at(i) = peak_fit.at(i);
  }

  return y_val;
}



//==============================================================================
// 2D lineshape function
//------------------------------------------------------------------------------



// [[Rcpp::export]]
std::vector < std::complex<double> > 
  lineshape_2d(Rcpp::NumericVector x_direct, Rcpp::NumericVector x_indirect,
               Rcpp::IntegerVector xi_direct, Rcpp::IntegerVector xi_indirect,
               Rcpp::NumericVector par, Rcpp::IntegerVector i_res,
               Rcpp::IntegerVector i_dim) {

  using namespace std;

  //--------------------
  // Double checking x vector lengths
  
  if ( xi_direct.size() != xi_indirect.size() ) {
    Rcpp::stop("Direct and indirect dimension vectors must be same length.");
  }

  //--------------------
  // Converting and storing data (to make use of generic fit functions)
  
  int n_res = max(i_dim) + 1;
  int n_points = xi_direct.size();
  int n_par = par.size();
  int n_peaks = n_par/4;

  // Converting Rcpp objects to std::vector to keep things standardized
  
  vector< double > x_1(x_direct.size(), 0);
  for (int i = 0; i < x_direct.size(); i++) {
    x_1.at(i) = x_direct.at(i);
  }

  vector< double > x_2(x_indirect.size(), 0);
  for (int i = 0; i < x_indirect.size(); i++) {
    x_2.at(i) = x_indirect.at(i);
  }

  data_lineshape data_direct[2];
  
  data_direct[0].x = x_1;
  data_direct[1].x = x_2;

  data_direct[0].peak_fit = 
    vector< vector < complex<double> > >
      ( n_res, vector< complex<double> > 
        ( n_points, complex<double> (0,0) ) );
  data_direct[1].peak_fit = 
    vector< vector < complex<double> > >
      ( n_res, vector< complex<double> > 
        ( n_points, complex<double> (0,0) ) );

  //-------------------- 
  // Loop over peaks
  
  for (int i = 0; i < n_peaks; i++) {

    double p = par[i*4];
    double w = par[i*4+1];
    double h = par[i*4+2];
    double f = par[i*4+3];

    // If the fraction gauss value is 0, proceed as Lorentz (w becomes wl)
    if ( f < 1e-6 ) {
		  lorentz(p, w, h, i_res.at(i), i*4, NULL, &data_direct[i_dim.at(i)]);	
    } else if ( f < (1 - 1e-6) ) {
		  voigt(p, w, h, f, i_res.at(i), i*4, NULL, &data_direct[i_dim.at(i)]);	
    } else if ( f < (1 + 1e-6) ) {
		  gauss(p, w, h, i_res.at(i), i*4, NULL, &data_direct[i_dim.at(i)]);	
    } else {
      Rcpp::stop("Fraction gauss outside valid [0,1] range. Aborting");
    }
  }
  
  //--------------------
  // For 2D peaks, the y value is the sum of products between direct and
  // indirect dimensions for matching resonances

  vector< complex<double> > y_val(n_points, complex<double>(0, 0));
  vector< complex<double> > peak_fit_1 = data_direct[0].peak_fit.at(0);

  for (int i = 0; i < n_points; i++) {
    y_val.at(i) = peak_fit.at(i);
  }

  return y_val;
}
