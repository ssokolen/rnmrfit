#include <vector>
#include <math.h>
#include <stdlib.h>
#include <Rcpp.h>
#include "rnmrfit.h"

using namespace Rcpp;




//==============================================================================
// 1D lineshape function
//------------------------------------------------------------------------------



// [[Rcpp::export]]
void lineshape_1d(NumericVector x_direct, NumericMatrix y, NumericVector par) {

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
      stop("Fraction gauss outside valid [0,1] range. Aborting");
    }
  }
  
  //--------------------
  // For 1D peaks, the data_direct.peak_fit can be read directly into y

  vector< complex<double> > peak_fit = data_direct[0].peak_fit.at(0);

  for (int i = 0; i < n_points; i++) {
    y(i, 0) = real(peak_fit.at(i));
    y(i, 1) = imag(peak_fit.at(i));
  }

  return;
}



//==============================================================================
// 2D lineshape function
//------------------------------------------------------------------------------



// [[Rcpp::export]]
void lineshape_2d(NumericVector x_direct, NumericVector x_indirect,
                  IntegerVector xi_direct, IntegerVector xi_indirect,
                  NumericMatrix y, NumericVector par, IntegerVector i_res,
                  IntegerVector i_dim) {

  using namespace std;

  //--------------------
  // Double checking x vector lengths

  if ( xi_direct.size() != xi_indirect.size() ) {
    stop("Direct and indirect dimension vectors must be same length.");
  }

  //--------------------
  // Converting and storing data (to make use of generic fit functions)

  int n_res = max(i_res) + 1;
  int n_points = xi_direct.size();
  int n_direct = x_direct.size();
  int n_indirect = x_indirect.size();
  int n_par = par.size();
  int n_peaks = n_par/4;

  // Converting Rcpp objects to std::vector to keep things standardized

  data_lineshape data_direct[2];


  data_direct[0].x = vector< double > (n_direct, 0);
  for (int i = 0; i < n_direct; i++) {
    data_direct[0].x.at(i) = x_direct.at(i);
  }
  data_direct[0].peak_fit = 
    vector< vector < complex<double> > >
      ( n_res, vector< complex<double> > 
        ( n_points, complex<double> (0,0) ) );

  data_direct[1].x = vector< double > (n_indirect, 0);
  for (int i = 0; i < n_indirect; i++) {
    data_direct[1].x.at(i) = x_indirect.at(i);
  }
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
      stop("Fraction gauss outside valid [0,1] range. Aborting");
    }
  }

  //--------------------
  // For 2D peaks, the y value is the sum of products between direct and
  // indirect dimensions for matching resonances
  
  // Referencing 2D data to 1D slices
  int i_direct, i_indirect;
  double dr, di, ir, ii;

  // Looping over points
  for (int i = 0; i < n_points; i++) {
    y(i, 0) = 0;
    y(i, 1) = 0;
    y(i, 2) = 0;
    y(i, 3) = 0;

    i_direct = xi_direct.at(i);
    i_indirect = xi_indirect.at(i);

    // Looping over resonances
    for (int j = 0; j < n_res; j++) {
      dr = real(data_direct[0].peak_fit.at(j).at(i_direct));
      di = imag(data_direct[0].peak_fit.at(j).at(i_direct));
      ir = real(data_direct[1].peak_fit.at(j).at(i_indirect));
      ii = imag(data_direct[1].peak_fit.at(j).at(i_indirect));

      y(i, 0) += dr * ir;
      y(i, 1) += dr * ii;
      y(i, 2) += ir * di;
      y(i, 3) += di * ii;
    }
  }

  return;
}
