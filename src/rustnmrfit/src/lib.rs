use nlopt::Algorithm;
use nlopt::{Nlopt, SuccessState, FailState};
use nlopt::Target;

use ndarray::prelude::*;

mod lineshape;
mod fit;

use fit::Fit1D;

/// Optimizes the fit of NMR peaks, baseline, and phase to y data using NLopt.
///
/// The y data is a 2D array with the first row encoding real components and the second
/// row imaginary ones.
///
/// WARNING: This function is meant to be used from the R package rnmrfit, YMMV.
pub fn fit_1d(x: Array<f64, Ix1>, y: Array<f64, Ix2>,
              p: Array<f64, Ix1>, lb: Array<f64, Ix1>, ub: Array<f64, Ix1>,
              nl: usize, nb: usize, np: usize, 
              basis: Option<Array<f64, Ix2>>) 
    -> (Array<f64, Ix1>, Result<(SuccessState, f64), (FailState, f64)>) {

    // Generate Fit object
    let fit = Fit1D::new(x, y, nl, nb, np, basis);

    // Define objective function
    fn eval(p: &[f64], grad: Option<&mut [f64]>, obj: &mut Fit1D) -> f64 {
        obj.eval(p, grad)
    }

    // Initializing nlopt object
    let n_par: usize = nl + nb*2 + np;
    let mut opt = Nlopt::new(Algorithm::Slsqp, n_par, eval, Target::Minimize, fit);

    // Fit requirements
    opt.set_maxeval(1000).unwrap();
    opt.set_xtol_rel(1e-4).unwrap();

    // Set simple bounds
    let lb = lb.to_vec();
    let ub = ub.to_vec();

    opt.set_lower_bounds(&lb[..]).unwrap();
    opt.set_upper_bounds(&ub[..]).unwrap();

    // Run the optimization
    let mut p = p.to_vec();
    let out = opt.optimize(&mut p[..]);

    // Reforming array from vector
    let p = Array::from_shape_vec((p.len(), ), p).unwrap();
    
    (p, out)
}

/// Evaluates the fit and gradient of NMR peaks, baseline, and phase.
fn eval_grad_1d(x: Array<f64, Ix1>, p: Array<f64, Ix1>, nl: usize, nb: usize, np: usize, 
                basis: Option<Array<f64, Ix2>>) 
    -> (Array<f64, Ix2>, Array<f64, Ix2>) {

    // Initialize output
    let n = x.len();
    let n_par = p.len();

    // Generate Fit object
    let y = Array::zeros((2, n));
    let mut fit = Fit1D::new(x, y, nl, nb, np, basis);

    // Evaluate fit using given parameters
    let mut grad = vec![0.0; n_par];
    fit.eval(&p.to_vec(), Some(&mut grad));

    // Extracting fit
    let mut y = Array::zeros((2, n));
    y.assign(&fit.y_fit);
    let mut grad = Array::zeros((2, n_par));
    grad.assign(&fit.grad);

    (y, grad)
}

/// Evaluates the fit of NMR peaks, baseline, and phase.
///
/// WARNING: This function is meant to be used from the R package rnmrfit, YMMV.
pub fn eval_1d(x: Array<f64, Ix1>, p: Array<f64, Ix1>, nl: usize, nb: usize, np: usize, 
           basis: Option<Array<f64, Ix2>>) 
    -> Array<f64, Ix2> {

    let (y, _) = eval_grad_1d(x, p, nl, nb, np, basis);
    y
}

/// Evaluates the gradient of NMR peaks, baseline, and phase.
///
/// WARNING: This function is meant to be used from the R package rnmrfit, YMMV.
pub fn grad_1d(x: Array<f64, Ix1>, p: Array<f64, Ix1>, nl: usize, nb: usize, np: usize, 
           basis: Option<Array<f64, Ix2>>) 
    -> Array<f64, Ix2> {

    let (_, grad) = eval_grad_1d(x, p, nl, nb, np, basis);
    grad
}
