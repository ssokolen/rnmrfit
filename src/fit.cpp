#include <vector>
#include <math.h>
#include <stdlib.h>
#include <Rcpp.h>
#include <nloptrAPI.h>
#include "Faddeeva.h"
#include "rnmrfit.h"


//' @importFrom Rcpp sourceCpp
//' @import nloptr
//' @useDynLib rnmrfit



//==============================================================================
// Constraint functions
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
double constrain_position(unsigned n, const double *x, 
                          double *grad, void *data) {

	using namespace std;

  // Unpacking data structure
  data_constraint *d = (data_constraint *) data;
  
  vector< bool > sign = d->sign;	
  vector< int > peak_number = d->peak_number_1;	
  double offset = d->offset;

  // Ensuring that the gradient is zero for all unaffected terms
  if ( grad ) {
    for (int i = 0; i < n; i++) {
      grad[i] = 0;
    }
  }

  // Evaluating
  double eval = 0;

  for (int i = 0; i < peak_number.size(); i++) {

    // Converting peak number to index
    int j = (peak_number.at(i) - 1)*4;

    // Extracting out sign
    bool neg = sign.at(i);

    if ( neg ) { eval -= x[j]; }
    else { eval += x[j]; }

    if (grad) {
      if ( neg ) { grad[j] = -1; }
      else { grad[j] = 1; }
    }
  }
  
  return (eval - offset);
}



//------------------------------------------------------------------------------
double constrain_width(unsigned n, const double *x, 
                      double *grad, void *data) {

	using namespace std;

  // Unpacking data structure
  data_constraint *d = (data_constraint *) data;
  
  vector< bool > sign = d->sign;	
  vector< int > peak_number = d->peak_number_1;	
  double offset = d->offset;

  // Ensuring that the gradient is zero for all unaffected terms
  if ( grad ) {
    for (int i = 0; i < n; i++) {
      grad[i] = 0;
    }
  }

  // Unlike position the width constraints divides the sum
  // of all "positive" indexes by the sum of all "negative" indexes so
  // it's necessary to do 2 passes: one to calculate the evaluation and another
  // to calculate the gradients.

  // First pass to calculate evaluation
  double numerator = 0;
  double denominator = 0;

  for (int i = 0; i < peak_number.size(); i++) {

    // Converting peak number to index of width
    int j = (peak_number.at(i) - 1)*4 + 1;

    // Extracting out sign
    bool neg = sign.at(i);

    if ( neg ) { denominator += x[j]; }
    else { numerator += x[j]; }
  }

  double eval = numerator/denominator;

  // Second pass to calculate gradients
  if ( grad ) {

    // All gradient terms are the same
    double grad_pos = 1/denominator;
    double grad_neg = -eval/denominator;

    for (int i = 0; i < peak_number.size(); i++) {

      // Converting peak number to index of width
      int j = (peak_number.at(i) - 1)*4 + 1;

      // Extracting out sign
      bool neg = sign.at(i);

      if ( neg ) { grad[j] = grad_neg; }
      else {  grad[j] = grad_pos; }
    }
  }
  
  return (eval - offset);
}

//------------------------------------------------------------------------------
double constrain_height(unsigned n, const double *x, 
                        double *grad, void *data) {

	using namespace std;

  // Unpacking data structure
  data_constraint *d = (data_constraint *) data;
  
  vector< bool > sign = d->sign;	
  vector< int > peak_number = d->peak_number_1;	
  double offset = d->offset;

  // Ensuring that the gradient is zero for all unaffected terms
  if ( grad ) {
    for (int i = 0; i < n; i++) {
      grad[i] = 0;
    }
  }

  // Unlike position the height constraints divides the sum
  // of all "positive" indexes by the sum of all "negative" indexes so
  // it's necessary to do 2 passes: one to calculate the evaluation and another
  // to calculate the gradients.

  // First pass to calculate evaluation
  double numerator = 0;
  double denominator = 0;

  for (int i = 0; i < peak_number.size(); i++) {

    // Converting peak number to index of height
    int j = (peak_number.at(i) - 1)*4 + 2;

    // Extracting out sign
    bool neg = sign.at(i);

    if ( neg ) { denominator += x[j]; }
    else { numerator += x[j]; }
  }

  double eval = numerator/denominator;

  // Second pass to calculate gradients
  if ( grad ) {

    // All gradient terms are the same
    double grad_pos = 1/denominator;
    double grad_neg = -eval/denominator;

    for (int i = 0; i < peak_number.size(); i++) {

      // Converting peak number to index of height
      int j = (peak_number.at(i) - 1)*4 + 2;

      // Extracting out sign
      bool neg = sign.at(i);

      if ( neg ) { grad[j] = grad_neg; }
      else {  grad[j] = grad_pos; }
    }
  }
  
  return (eval - offset);
}

//------------------------------------------------------------------------------
double constrain_fraction(unsigned n, const double *x, 
                          double *grad, void *data) {

	using namespace std;

  // Unpacking data structure
  data_constraint *d = (data_constraint *) data;
  
  vector< bool > sign = d->sign;	
  vector< int > peak_number = d->peak_number_1;	
  double offset = d->offset;

  // Ensuring that the gradient is zero for all unaffected terms
  if ( grad ) {
    for (int i = 0; i < n; i++) {
      grad[i] = 0;
    }
  }

  // Evaluating
  double eval = 0;

  for (int i = 0; i < peak_number.size(); i++) {

    // Converting peak number to index offset to fraction
    int j = (peak_number.at(i) - 1)*4 + 3;

    // Extracting out sign
    bool neg = sign.at(i);

    if ( neg ) { eval -= x[j]; }
    else { eval += x[j]; }

    if (grad) {
      if ( neg ) { grad[j] = -1; }
      else { grad[j] = 1; }
    }
  }
  
  return (eval - offset);
}


//------------------------------------------------------------------------------
double constrain_area(unsigned n, const double *x, 
                      double *grad, void *data) {

	using namespace std;

  // Unpacking data structure
  data_constraint *d = (data_constraint *) data;
  
  vector< bool > sign = d->sign;	
  vector< int > peak_number = d->peak_number_1;	
  double offset = d->offset;

  // Ensuring that the gradient is zero for all unaffected terms
  if ( grad ) {
    for (int i = 0; i < n; i++) {
      grad[i] = 0;
    }
  }

  // Unlike position the area constraints divides the sum
  // of all "positive" indexes by the sum of all "negative" indexes so
  // it's necessary to do 2 passes: one to calculate the evaluation and another
  // to calculate the gradients.

  // First pass to calculate evaluation
  double numerator = 0;
  double denominator = 0;

  for (int i = 0; i < peak_number.size(); i++) {

    // Converting peak number to index
    int j = (peak_number.at(i) - 1)*4;

    double w = x[j+1];
    double h = x[j+2];
    double f = x[j+3];
    double area = 0;
    double wg;
    complex<double> z, F;

    // Area calculation depends on fraction
    if ( f < 1e-6 ) { 
      area = M_PI * w * h; 
    } else if ( f < (1-1e-6) ) { 
      wg = w*f/(1 - f);
      z = complex<double> (0, w/(sqrt(2) * wg));
      F = Faddeeva::w( z );
      area = sqrt(2 * M_PI) * wg * h / real(F);
    } else {
      area = sqrt(2 * M_PI) * w * h;
    }

    // Extracting out sign
    bool neg = sign.at(i);

    if ( neg ) { denominator += area; }
    else { numerator += area; }
  }

  double eval = numerator/denominator;

  // Second pass to calculate gradients
  if ( grad ) {

    // Calculating some common terms
    double dfda_pos = 1/denominator;
    double dfda_neg = -eval/denominator;
    double dfda = 0;

    for (int i = 0; i < peak_number.size(); i++) {

      // Converting peak number to index
      int j = (peak_number.at(i) - 1)*4;

      double w = x[j+1];
      double h = x[j+2];
      double f = x[j+3];
      double wg;
      complex<double> z, F;

      // Extracting out sign
      bool neg = sign.at(i);

      if ( neg ) { dfda = dfda_neg; }
      else { dfda = dfda_pos; }

      // Gradient calculations depend on fraction
      if ( f < 1e-6 ) {
        // Width
        grad[j+1] = dfda * M_PI * h;
        // Height
        grad[j+2] = dfda * M_PI * w;
      } else if ( f < (1-1e-6) ) { 
        wg = w*f/(1-f);
        z = complex<double> (0, w/(sqrt(2) * wg));
        F = Faddeeva::w( z );
        // Width
        grad[j+1] = dfda * 2 * w * f * sqrt(2 * M_PI) / ( real(F) * (1-f) );
        // Height
        grad[j+2] = dfda * sqrt(2 * M_PI) * wg / real(F);
        // Fraction gauss
        grad[j+3] = dfda * w*w * sqrt(2 * M_PI) * 
                    ( 1/real(F*(1-f)*(1-f)) + 1/real(F*f*f) - 
                      2/real(F*F*(1-f)*f*sqrt( 2 * M_PI )) );
      } else if ( f < (1+1e-6) ) { 
        // Width
        grad[j+1] = dfda * sqrt(2 * M_PI) * h;
        // Height
        grad[j+2] = dfda * sqrt(2 * M_PI) * w;
      } else {
        Rcpp::stop("Fraction gauss fell out of [0,1] range. Aborting");
      }
    }
  }
  
  return (eval - offset);
}




//==============================================================================
// Lineshape functions
//------------------------------------------------------------------------------
// All lineshape functions take a series of peak parameters and use those
// parameters to increment fit data in a data_lineshape structure. To reuse the
// same functions for 1D and 2D analysis, the fits are divided by resonance
// index, which should be kept at 0 for 1D usage.



//------------------------------------------------------------------------------
void lorentz(double p, double wl, double h, 
             int i_res, int i_par, double *grad, void *data) {

	using namespace std;

  // Unpacking data structure
  data_lineshape *d = (data_lineshape *) data;
  
  vector< double > x = d->x;
  
  vector< complex<double> > *peak_fit = &(d->peak_fit.at(i_res));
  vector< vector < complex<double> > > *peak_partial = &(d->peak_partial);

  // Main terms
  complex<double> fo, dfdz, dfdp, dfdwl, dzdwl, dzdp;

  // Commonly used intermediate values
  double z, z2, z2p1;

  // Looping through each x value
  for(int i = 0; i < x.size(); i++ ) {

    // Pre-calculating common values
    z = (p - x.at(i))/wl;
    z2 = z*z;
    z2p1 = z2 + 1;

    // Calculating an intermediate f with no h term
    fo = ( complex<double> (1, z) ) / 
         ( complex<double> (z2p1, 0) );
       
    // Tacking on the h
    (*peak_fit).at(i) += ( complex<double> (h, 0) ) * fo;

    // Only compute derivatives if required
    if ( grad ) {
      dfdz = ( complex<double> (-2*z, 1-z2) ) * 
             ( complex<double> (h/(z2p1*z2p1), 0) );

      // Derivatives of z for chain rule
      dzdp = complex<double> (1/wl, 0);
      dzdwl = complex<double> (-z/wl, 0);

      // Position gradient
      dfdp = dfdz*dzdp;
      (*peak_partial).at(i).at(i_par) = dfdp;

      // Lorentz width gradient
      dfdwl = dfdz*dzdwl;
      (*peak_partial).at(i).at(i_par+1) = dfdwl;

      // Height gradient
      (*peak_partial).at(i).at(i_par+2) = fo;

      // Fraction gradient
      (*peak_partial).at(i).at(i_par+3) = 0;
    }
  }

  return; 
}

      

//------------------------------------------------------------------------------
void gauss(double p, double wg, double h, 
             int i_res, int i_par, double *grad, void *data) {

	using namespace std;

  // Unpacking data structure
  data_lineshape *d = (data_lineshape *) data;
  
  vector< double > x = d->x;
  
  vector< complex<double> > *peak_fit = &(d->peak_fit.at(i_res));
  vector< vector < complex<double> > > *peak_partial = &(d->peak_partial);

  // Main terms
  complex<double> fo, dfdz, dfdp, dfdwg, dzdwg, dzdp;

  // Commonly used intermediate values
  double z;

  // Looping through each x value
  for(int i = 0; i < x.size(); i++ ) {

    // Main terms
    z = (p - x.at(i))/(sqrt(2.0)*wg);
    
    // Calculating an intermediate f with no h term
    fo = Faddeeva::w( complex<double> (z, 0) );
       
    // Tacking on the h
    (*peak_fit).at(i) += ( complex<double> (h, 0) ) * fo;
    
    // Only compute derivatives if required
    if ( grad ) {
      dfdz = ( complex<double> (h, 0) ) *
             ( complex<double> ( -2*z, 0 ) * fo + 
               complex<double> ( 0, 2/sqrt(M_PI) ) );

      // Derivatives of z for chain rule
      dzdp = complex<double> (1/(sqrt(2.0)*wg), 0);
      dzdwg = complex<double> (-z/wg, 0);

      // Position gradient
      dfdp = dfdz * dzdp;
      (*peak_partial).at(i).at(i_par) = dfdp;

      // Gauss width gradient
      dfdwg = dfdz * dzdwg;
      (*peak_partial).at(i).at(i_par+1) = dfdwg;

      // Height gradient
      (*peak_partial).at(i).at(i_par+2) = fo;

      // Fraction gradient
      (*peak_partial).at(i).at(i_par+3) = 0;
    }
  }

  return; 
}



void voigt(double p, double wl, double h, double f, 
           int i_res, int i_par, double *grad, void *data) {

	using namespace std;

  // Unpacking data structure
  data_lineshape *d = (data_lineshape *) data;
  
  vector< double > x = d->x;
  
  vector< complex<double> > *peak_fit = &(d->peak_fit.at(i_res));
  vector< vector < complex<double> > > *peak_partial = &(d->peak_partial);

  // Main terms
  complex<double> z, f0, dfdz, dfdp, dfdwl, dfdf, dzdp, dzdwl, dzdf;

  // Scaling terms
  complex<double> zn, fn, fn2, dfndzn, dzndf;

  // Pre-calculations
  double wg = wl * f / (1 - f);

  //Normalising factors
  zn = complex<double> ( 0, wl / (sqrt(2.0) * wg) );
  fn = Faddeeva::w( zn );
  fn2 = fn*fn;

  dfndzn = ( complex<double> ( -2, 0 ) * zn * fn + 
             complex<double> ( 0, 2/sqrt(M_PI) ) );  

  dzndf = -zn * wg / (wl * f*f);

  // Looping through each x value
  for(int i = 0; i < x.size(); i++ ) {

    // Pre-calculating common values
    z = ( complex<double> ( p - x.at(i) , wl ) ) / 
        ( complex<double> ( (sqrt(2.0) * wg), 0 ) );

    // Calculating an intermediate f with no h term
    f0 = Faddeeva::w( z );
       
    // Tacking on the h
    (*peak_fit).at(i) += ( complex<double> (h, 0) ) * f0 / fn;

    // Only compute derivatives if required
    if ( grad ) {
      
      dfdz = ( complex<double> ( -2, 0 ) * z * f0 + 
               complex<double> ( 0, 2/sqrt(M_PI) ) );  

      // Derivatives of z for chain rule
      dzdp = complex<double> ( 1/( sqrt(2.0)*wg ), 0 );
      dzdwl = complex<double> ( 0, 1/( sqrt(2.0)*wg ) );
      dzdf = -z * wg / (wl * f*f);

      // Position gradient
      dfdp = dfdz * dzdp * h / fn;
      (*peak_partial).at(i).at(i_par) = dfdp;

      // Lorentz width gradient
      dfdwl = h * (dfdz * fn - dfndzn * f0) / fn2 * dzdwl;
      (*peak_partial).at(i).at(i_par+1) = dfdwl;
      
	    // Height gradient
      (*peak_partial).at(i).at(i_par+2) = f0 / fn;

      // Fraction gradient
      dfdf = h * (dfdz * dzdf * fn - dfndzn * dzndf * f0) / fn2;
      (*peak_partial).at(i).at(i_par+3) = dfdf;
      
    }
  }

  return; 
}



//==============================================================================
// 1D fitting

//------------------------------------------------------------------------------
// Objective function

double f_obj_1d(unsigned n, const double *par, double *grad, void *data) {
    
	using namespace std;

  // Unpacking lineshape_data structure
  data_1d *d = (data_1d *) data;

  vector< double > x = d->x;
  vector< vector<double> > y = d->y;
  vector< double > *y_mod = &(d->y_mod);
  vector< double > *y_fit = &(d->y_fit);
  vector< double > *y_dif = &(d->y_dif);

  vector< vector<double> > *basis = &(d->basis);
  vector< double > *basis_row;

  int n_par = d->n_par;
  int n_peaks = d->n_peaks;
  int n_baseline = d->n_baseline;
  int n_phase = d->n_phase;
  int *count = &(d->count);

  // Whereas n_peaks refers to the number of peak, the number of peak parameters
  // are also frequentely needed
  int n_peaks_par = n_peaks*4;

  data_lineshape *data_direct = &(d->lineshape);
  data_lineshape *d1 = (data_lineshape *) data_direct;

  vector< complex<double> > *peak_fit = &(d1->peak_fit.at(0));
  vector< vector < complex<double> > > *peak_partial = &(d1->peak_partial);

  *count += 1;

  //--------------------
  // Calculate lineshape peak_fit, and gradient partial derivatives

  // Reseting peak fits to zero to ensure proper accumulation
  for(int i = 0; i < (*peak_fit).size(); i++ ) {
    (*peak_fit).at(i) = complex<double> (0, 0);
  }
  
  // Loop over peaks
  for (int i = 0; i < n_peaks; i++) {


    double p = par[i*4];
    double w = par[i*4+1];
    double h = par[i*4+2];
    double f = par[i*4+3];

    Rcpp::Rcout << i << ", " << n_peaks << ", " << f << std::endl;

    // If the fraction gauss value is 0, proceed as Lorentz (w becomes wl)
    if ( f < 1e-6 ) {
		  lorentz(p, w, h, 0, i*4, grad, data_direct);	
    } else if ( f < (1 - 1e-6) ) {
		  voigt(p, w, h, f, 0, i*4, grad, data_direct);	
    } else if ( f < (1 + 1e-6) ) {
		  gauss(p, w, h, 0, i*4, grad, data_direct);	
    } else {
      Rcpp::stop("Fraction gauss fell out of [0,1] range. Aborting");
    }
  }

  //--------------------
	// Sum of squares
	double eval = 0;

  // For phase correction
  double theta = 0;

  // Ensuring that all the gradients start at 0 for accumulation
  if ( grad ) {
    for (int j = 0; j < n_par; j++) {
      grad[j] = 0;
    }
  }
  
  // For 1D fit, assume that indexes of y and lineshape data match up
  for (int i = 0; i < y.size(); i++) {

    // Applying phase correction if necessary
    if ( n_phase > 0 ) { 

      if ( n_phase == 1 ) {
        theta = par[n_par-1];
      } else {
        theta = par[n_par-2] + x.at(i)*par[n_par-1];
      }

      (*y_mod).at(0) = y.at(i).at(0)*cos(theta) + y.at(i).at(1)*sin(theta); 
      (*y_mod).at(1) = -y.at(i).at(0)*sin(theta) + y.at(i).at(1)*cos(theta); 

    } else { 
      (*y_mod).at(0) = y.at(i).at(0);
      (*y_mod).at(1) = y.at(i).at(1);
    }

    // Applying peak fit
    (*y_fit).at(0) = (*peak_fit).at(i).real();
    (*y_fit).at(1) = (*peak_fit).at(i).imag();

    // Applying baseline fit
    basis_row = &(basis->at(i));
    for (int j = 0; j < n_baseline; j++) {
      (*y_fit).at(0) += par[n_peaks_par+j]*(*basis_row).at(j);
      (*y_fit).at(1) += par[n_peaks_par+j+n_baseline]*(*basis_row).at(j);
    }

    // Calculating difference
    (*y_dif).at(0) = (*y_mod).at(0) - (*y_fit).at(0);
    (*y_dif).at(1) = (*y_mod).at(1) - (*y_fit).at(1);;

    // Incrementing sum of squares
	  eval += (*y_dif).at(0) * (*y_dif).at(0);
	  eval += (*y_dif).at(1) * (*y_dif).at(1);

    // Filling in gradients if necessary
    // Note that the 2 scaling term is applied outside the loop
    if ( grad ) {

      // First the peaks
      for (int j = 0; j < n_peaks_par; j++) {
        grad[j] -= (*y_dif).at(0) * (*peak_partial).at(i).at(j).real();
        grad[j] -= (*y_dif).at(1) * (*peak_partial).at(i).at(j).imag();
      }

      // Then the baseline
      for (int j = 0; j < n_baseline; j++) {
        grad[n_peaks_par+j] -= (*y_dif).at(0) * (*basis_row).at(j);
        grad[n_peaks_par+j+n_baseline] -= (*y_dif).at(1) * (*basis_row).at(j);
      }

      // Finally, the phase
      if ( n_phase > 0 ) {

        double dtheta = (*y_dif).at(0) * 
                        (-y.at(i).at(0)*sin(theta) + y.at(i).at(1)*cos(theta)) +
                        (*y_dif).at(1) * 
                        (-y.at(i).at(0)*cos(theta) - y.at(i).at(1)*sin(theta));

        if ( n_phase == 1 ) {
          grad[n_par-1] += dtheta; 
        } else {
          grad[n_par-2] += dtheta;
          grad[n_par-1] += dtheta * x.at(i);
        }
      }
    }
  }

  // Multiplying gradients by 2
  if ( grad ) {
    for (int j = 0; j < n_par; j++) {
      grad[j] *= 2;
    }
  }

  return eval;	
}



//------------------------------------------------------------------------------
// Main function called from R

//' @export
// [[Rcpp::export]]
double fit_lineshape_1d(
  const Rcpp::NumericVector x_val, const Rcpp::ComplexVector y_val,
  Rcpp::NumericVector par, Rcpp::NumericVector lb, Rcpp::NumericVector ub,
  Rcpp::NumericMatrix basis_val, Rcpp::List eq, Rcpp::List ineq, 
  int n_peaks, int n_baseline, int n_phase) {

  using namespace std;

  //--------------------
  // Converting and storing data
  
  int n_points = y_val.size();
  int n_par = par.size();

  // x and y data  
  vector< vector<double> > y(n_points, vector<double> (2,0) );
  vector< double > x(n_points, 0);

  vector< double > y_mod(2,0);
  vector< double > y_fit(2,0);
  vector< double > y_dif(2,0);

  std::vector< std::vector<double> > basis(n_points, 
      vector<double> (n_baseline,0) );	

  Rcpp::Function R_re("Re");
  Rcpp::Function R_im("Im");

  // Converting Rcpp objects to std::vector to keep things standardized
  // For 1D fit, x and y can be treated as paired
  for (int i = 0; i < n_points; i++) {
    y.at(i).at(0) = Rcpp::as<double>(R_re(y_val.at(i)));
    y.at(i).at(1) = Rcpp::as<double>(R_im(y_val.at(i)));

    x.at(i) = x_val.at(i);

    for (int j = 0; j < n_baseline; j++) {
      basis.at(i).at(j) = basis_val(i, j);  
    }

  }

  // First, packing lineshape data
  data_lineshape data_direct;
  
  data_direct.x = x;

  data_direct.peak_fit = 
    vector< vector < complex<double> > >
      ( 1, vector< complex<double> > 
        ( n_points, complex<double> (0,0) ) );
  
  data_direct.peak_partial = 
    vector< vector < complex<double> > > 
      ( n_points, vector< complex<double> > 
        ( n_par, complex<double> (0,0) ) );

  // And then packing general data
  data_1d data[1];

  data[0].lineshape = data_direct;
  data[0].x = x;
  data[0].y = y;
  data[0].y_mod = y_mod;
  data[0].y_fit = y_fit;
  data[0].y_dif = y_dif;
  data[0].basis = basis;
  data[0].n_par = n_par;
  data[0].n_peaks = n_peaks;
  data[0].n_baseline = n_baseline;
  data[0].n_phase = n_phase;
  data[0].count = 0;

  //--------------------
  // Initializing the optimizer

  nlopt_opt opt;

  opt = nlopt_create(NLOPT_LD_SLSQP, n_par);    
  // opt = nlopt_create(NLOPT_LN_COBYLA, n_par); 
  nlopt_set_min_objective(opt, f_obj_1d, &data[0]);
  nlopt_set_xtol_rel(opt, 1e-4);

  // Bounds
  nlopt_set_lower_bounds(opt, &(lb[0]));
  nlopt_set_upper_bounds(opt, &(ub[0]));

  //--------------------
  // Adding constraints

  // First the equality constraints
  int n_eq = eq.size();
  data_constraint data_eq[n_eq];

  for (int i = 0; i < n_eq; i++) {

    vector< double > constraints = eq.at(i); 
    int n_terms = constraints.size() - 2;

    data_eq[i].offset = constraints.at(1);

    data_eq[i].sign = vector< bool > ( n_terms, 0);
    data_eq[i].peak_number_1 = vector< int > ( n_terms, 0);

    for (int j = 0; j < n_terms; j++) {
      bool neg = signbit(constraints.at(j+2));
      data_eq[i].sign[j] = neg;
      if ( neg ) { data_eq[i].peak_number_1[j] = -round(constraints.at(j+2)); }
      else { data_eq[i].peak_number_1[j] = round(constraints.at(j+2)); }
    }

    // Adding constraint based on flag
    int flag = round(constraints.at(0));
    if ( flag == 0 ) {
      nlopt_add_equality_constraint(opt, constrain_position, &data_eq[i], 1e-8);
    } else if ( flag == 1) {
      nlopt_add_equality_constraint(opt, constrain_width, &data_eq[i], 1e-8);
    } else if ( flag == 2) {
      nlopt_add_equality_constraint(opt, constrain_height, &data_eq[i], 1e-8);
    } else if ( flag == 3) {
      nlopt_add_equality_constraint(opt, constrain_fraction, &data_eq[i], 1e-8);
    } else if ( flag == 4) {
      nlopt_add_equality_constraint(opt, constrain_area, &data_eq[i], 1e-8);
    } else {
      Rcpp::stop("Encountered undefined equality constraint. Aborting.");
    }
  }

  // The inequality constraint proceed exactly the same as the equality ones
  int n_ineq = ineq.size();
  data_constraint data_ineq[n_ineq];

  for (int i = 0; i < n_ineq; i++) {

    vector< double > constraints = ineq.at(i); 
    int n_terms = constraints.size() - 2;

    data_ineq[i].offset = constraints.at(1);

    data_ineq[i].sign = vector< bool > ( n_terms, 0);
    data_ineq[i].peak_number_1 = vector< int > ( n_terms, 0);

    for (int j = 0; j < n_terms; j++) {
      bool neg = signbit(constraints.at(j+2));
      data_ineq[i].sign[j] = neg;
      if ( neg ) { data_ineq[i].peak_number_1[j] = -round(constraints.at(j+2));}
      else { data_ineq[i].peak_number_1[j] = round(constraints.at(j+2)); }
    }

    // Adding constraint based on flag
    int flag = round(constraints.at(0));
    if ( flag == 0 ) {
      nlopt_add_inequality_constraint(opt, constrain_position, 
                                      &data_ineq[i], 1e-8);
    } else if ( flag == 1) {
      nlopt_add_inequality_constraint(opt, constrain_width, 
                                      &data_ineq[i], 1e-8);
    } else if ( flag == 2) {
      nlopt_add_inequality_constraint(opt, constrain_height, 
                                      &data_ineq[i], 1e-8);
    } else if ( flag == 3) {
      nlopt_add_inequality_constraint(opt, constrain_fraction, 
                                      &data_ineq[i], 1e-8);
    } else if ( flag == 4) {
      nlopt_add_inequality_constraint(opt, constrain_area, 
                                      &data_ineq[i], 1e-8);
    } else {
      Rcpp::stop("Encountered undefined inequality constraint. Aborting.");
    }
  }
    
  // Minimum objective value, upon return
  double minf; 						

  //--------------------
  // Running the fit

  data_1d *d = (data_1d *) data;
  int *count = &(d->count);

  if ( nlopt_optimize(opt, &(par[0]), &minf ) < 0) {
    Rcpp::Rcout << "nlopt failed!" << std::endl;
  } else {
    Rcpp::Rcout << "nlopt succeeded!" << std::endl;
    	//Checking minf
	Rcpp::Rcout << "Found minimum with" << std::endl;  
    Rcpp::Rcout << "par(" << par[0] << ", " << par[1] << ", " << par[2] << ", " << par[3] << ")" << std::endl;
  }
	Rcpp::Rcout << "Iterations: " << *count << std::endl;  
  nlopt_destroy(opt);

	return minf;
}








       
