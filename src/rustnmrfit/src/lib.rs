use nlopt::Algorithm;
use nlopt::{Nlopt, SuccessState, FailState};
use nlopt::Target;

use ndarray::prelude::*;

mod peak;
mod lineshape;
mod phase;
mod fit;
mod constraint;

use fit::{Fit1D};
use constraint::Constraint;
pub use peak::{Peak, Lorentz, Voigt}; 

/// Optimizes the fit of NMR peaks, baseline, and phase to y data using NLopt.
///
/// The y data is a 2D array with the first row encoding real components and the second
/// row imaginary ones.
///
/// WARNING: This function is meant to be used from the R package rnmrfit, YMMV.
pub fn fit_1d(x: Array<f64, Ix1>, y: Array<f64, Ix2>,
              p: Array<f64, Ix1>, lb: Array<f64, Ix1>, ub: Array<f64, Ix1>,
              nl: usize, nb: usize, np: usize,
              basis: Option<Array<f64, Ix2>>,
              eq: Option<Vec<(usize, f64, Vec<usize>, Vec<usize>)>>,
              iq: Option<Vec<(usize, f64, Vec<usize>, Vec<usize>)>>)

    -> (Array<f64, Ix1>, Result<(SuccessState, f64), (FailState, f64)>) {

    // Generate Fit object
    let fit = Fit1D::new(x, y, nl, nb, np, basis);

    // Initializing nlopt object
    let n_par: usize = nl + nb*2 + np;
    let mut opt = Nlopt::new(Algorithm::Slsqp, n_par, Fit1D::obj, Target::Minimize, fit);

    // Fit requirements
    opt.set_maxeval(1000).unwrap();
    opt.set_xtol_rel(1e-4).unwrap();

    // Set simple bounds
    let lb = lb.to_vec();
    let ub = ub.to_vec();

    opt.set_lower_bounds(&lb[..]).unwrap();
    opt.set_upper_bounds(&ub[..]).unwrap();

    // Add equality constraints
    if eq.is_some() {
        let eq = eq.unwrap();

        for i in 0 .. eq.len() {
        
            let flag = eq[i].0;
            let offset = eq[i].1;
            let lhs = eq[i].2.clone();
            let rhs = eq[i].3.clone();

            let con = Constraint::new(lhs, rhs, offset);
        
            if flag == 0 {
                opt.add_equality_constraint(Constraint::position, con, 1e-8).unwrap();
            } else if flag == 1 {
                opt.add_equality_constraint(Constraint::width, con, 1e-8).unwrap();
            } else if flag == 2 {
                opt.add_equality_constraint(Constraint::height, con, 1e-8).unwrap();
            } else if flag == 3 {
                opt.add_equality_constraint(Constraint::fraction, con, 1e-8).unwrap();
            } else if flag == 4 {
                opt.add_equality_constraint(Constraint::area, con, 1e-8).unwrap();
            } else {
                panic!("Constraint flag must be 0-4, inclusive, aborting.")
            }
        }
    }

    // Add inequality constraints
    if iq.is_some() {
        let iq = iq.unwrap();

        for i in 0 .. iq.len() {
        
            let flag = iq[i].0;
            let offset = iq[i].1;
            let lhs = iq[i].2.clone();
            let rhs = iq[i].3.clone();

            let con = Constraint::new(lhs, rhs, offset);
        
            if flag == 0 {
                opt.add_inequality_constraint(Constraint::position, con, 1e-8).unwrap();
            } else if flag == 1 {
                opt.add_inequality_constraint(Constraint::width, con, 1e-8).unwrap();
            } else if flag == 2 {
                opt.add_inequality_constraint(Constraint::height, con, 1e-8).unwrap();
            } else if flag == 3 {
                opt.add_inequality_constraint(Constraint::fraction, con, 1e-8).unwrap();
            } else if flag == 4 {
                opt.add_inequality_constraint(Constraint::area, con, 1e-8).unwrap();
            } else {
                panic!("Constraint flag must be 0-4, inclusive, aborting.")
            }
        }
    }

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

/*
/// Optimizes the fit of NMR peaks, baseline, and phase to y data using NLopt.
///
/// The y data is a 2D array with the first row encoding rr, ri, ir, and ii component in 4 rows
///
/// 
/// WARNING: This function is meant to be used from the R package rnmrfit, YMMV.
pub fn fit_2d(x_direct: Array<f64, Ix1>, x_indirect: Array<f64, Ix1>, y: Array<f64, Ix2>,
              xi_direct: Array<usize, Ix1>, xi_indirect: Array<usize, Ix1>,
              peak_resonances: Array<usize, Ix1>, peak_dimensions: Array<usize, Ix1>,
              p: Array<f64, Ix1>, lb: Array<f64, Ix1>, ub: Array<f64, Ix1>,
              nl: usize, nb: usize, np: usize,
              _: Option<Array<f64, Ix2>>,
              eq: Option<Vec<(usize, f64, Vec<usize>, Vec<usize>)>>,
              iq: Option<Vec<(usize, f64, Vec<usize>, Vec<usize>)>>)

    -> (Array<f64, Ix1>, Result<(SuccessState, f64), (FailState, f64)>) {

    // Generate Fit object
    let fit = Fit2D::new(x_direct, x_indirect, y, xi_direct, xi_indirect, 
                         peak_resonances, peak_dimensions, nl, nb, np, _);

    // Initializing nlopt object
    let n_par: usize = nl + nb*2 + np;
    let mut opt = Nlopt::new(Algorithm::Slsqp, n_par, Fit2D::obj, Target::Minimize, fit);

    // Fit requirements
    opt.set_maxeval(1000).unwrap();
    opt.set_xtol_rel(1e-4).unwrap();

    // Set simple bounds
    let lb = lb.to_vec();
    let ub = ub.to_vec();

    opt.set_lower_bounds(&lb[..]).unwrap();
    opt.set_upper_bounds(&ub[..]).unwrap();

    // Add equality constraints
    if eq.is_some() {
        let eq = eq.unwrap();

        for i in 0 .. eq.len() {
        
            let flag = eq[i].0;
            let offset = eq[i].1;
            let lhs = eq[i].2.clone();
            let rhs = eq[i].3.clone();

            let con = Constraint::new(lhs, rhs, offset);
        
            if flag == 0 {
                opt.add_equality_constraint(Constraint::position, con, 1e-8).unwrap();
            } else if flag == 1 {
                opt.add_equality_constraint(Constraint::width, con, 1e-8).unwrap();
            } else if flag == 2 {
                opt.add_equality_constraint(Constraint::height, con, 1e-8).unwrap();
            } else if flag == 3 {
                opt.add_equality_constraint(Constraint::fraction, con, 1e-8).unwrap();
            } else if flag == 4 {
                opt.add_equality_constraint(Constraint::area, con, 1e-8).unwrap();
            } else {
                panic!("Constraint flag must be 0-4, inclusive, aborting.")
            }
        }
    }

    // Add inequality constraints
    if iq.is_some() {
        let iq = iq.unwrap();

        for i in 0 .. iq.len() {
        
            let flag = iq[i].0;
            let offset = iq[i].1;
            let lhs = iq[i].2.clone();
            let rhs = iq[i].3.clone();

            let con = Constraint::new(lhs, rhs, offset);
        
            if flag == 0 {
                opt.add_inequality_constraint(Constraint::position, con, 1e-8).unwrap();
            } else if flag == 1 {
                opt.add_inequality_constraint(Constraint::width, con, 1e-8).unwrap();
            } else if flag == 2 {
                opt.add_inequality_constraint(Constraint::height, con, 1e-8).unwrap();
            } else if flag == 3 {
                opt.add_inequality_constraint(Constraint::fraction, con, 1e-8).unwrap();
            } else if flag == 4 {
                opt.add_inequality_constraint(Constraint::area, con, 1e-8).unwrap();
            } else {
                panic!("Constraint flag must be 0-4, inclusive, aborting.")
            }
        }
    }

    // Run the optimization
    let mut p = p.to_vec();
    let out = opt.optimize(&mut p[..]);

    // Reforming array from vector
    let p = Array::from_shape_vec((p.len(), ), p).unwrap();
    
    (p, out)
}

/// Evaluates the fit and gradient of NMR peaks, baseline, and phase.
fn eval_grad_2d(x_direct: Array<f64, Ix1>, x_indirect: Array<f64, Ix1>, 
                xi_direct: Array<usize, Ix1>, xi_indirect: Array<usize, Ix1>,
                peak_resonances: Array<usize, Ix1>, peak_dimensions: Array<usize, Ix1>,
                p: Array<f64, Ix1>, nl: usize, nb: usize, np: usize,
                _: Option<Array<f64, Ix2>>)
    -> (Array<f64, Ix2>, Array<f64, Ix2>) {

    // Initialize output
    let n = x.len();
    let n_par = p.len();

    // Generate Fit object
    let y = Array::zeros((4, n));
    let mut fit = Fit2D::new(x_direct, x_indirect, y, xi_direct, xi_indirect, 
                             peak_resonances, peak_dimensions, nl, nb, np, _);

    // Evaluate fit using given parameters
    let mut grad = vec![0.0; n_par];
    fit.eval(&p.to_vec(), Some(&mut grad));

    // Extracting fit
    let mut y = Array::zeros((4, n));
    y.assign(&fit.y_fit);
    let mut grad = Array::zeros((4, n_par));
    grad.assign(&fit.grad);

    (y, grad)
}

/// Evaluates the fit of NMR peaks, baseline, and phase.
///
/// WARNING: This function is meant to be used from the R package rnmrfit, YMMV.
pub fn eval_2d(x_direct: Array<f64, Ix1>, x_indirect: Array<f64, Ix1>, 
               xi_direct: Array<usize, Ix1>, xi_indirect: Array<usize, Ix1>,
               peak_resonances: Array<usize, Ix1>, peak_dimensions: Array<usize, Ix1>,
               p: Array<f64, Ix1>, nl: usize, nb: usize, np: usize,
               _: Option<Array<f64, Ix2>> 
    -> Array<f64, Ix2> {

    let (y, _) = eval_grad_2d(x_direct, x_indirect, xi_direct, xi_indirect,
                              peak_resonances, peak_dimensions, p, nl, nb, np, _);
    y
}

/// Evaluates the gradient of NMR peaks, baseline, and phase.
///
/// WARNING: This function is meant to be used from the R package rnmrfit, YMMV.
pub fn grad_2d(x_direct: Array<f64, Ix1>, x_indirect: Array<f64, Ix1>, 
               xi_direct: Array<usize, Ix1>, xi_indirect: Array<usize, Ix1>,
               peak_resonances: Array<usize, Ix1>, peak_dimensions: Array<usize, Ix1>,
               p: Array<f64, Ix1>, nl: usize, nb: usize, np: usize,
               _: Option<Array<f64, Ix2>> 
    -> Array<f64, Ix2> {

    let (_, grad) = eval_grad_2d(x_direct, x_indirect, xi_direct, xi_indirect,
                              peak_resonances, peak_dimensions, p, nl, nb, np, _);
    grad
}
*/
