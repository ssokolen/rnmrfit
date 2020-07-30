use nlopt::{Algorithm, Target, Nlopt, ObjFn, SuccessState, FailState};

use ndarray::prelude::*;

mod common;
mod constraint;
mod peak;
mod lineshape;
mod baseline;
mod phase;
mod fit;

pub use peak::{Peak, Lorentz, Voigt}; 

use common::{NMRFitComponent, NMRFit};
use constraint::Constraint;
use lineshape::{Lineshape1D, Lineshape2D};
use baseline::{Baseline1D, Baseline2D};
use phase::{Phase1D, Phase2D};
use fit::{Fit1D, Fit2D};

//==============================================================================
pub fn add_constraints<F: ObjFn<T>, T>(opt: &mut Nlopt<F, T>, lb: Array1<f64>, ub: Array1<f64>,
                                       eq: Option<Vec<(usize, f64, Vec<usize>, Vec<usize>)>>,
                                       iq: Option<Vec<(usize, f64, Vec<usize>, Vec<usize>)>>) {

    // Basic constraints
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

}

//==============================================================================
pub fn fit_1d(x: Array1<f64>, y: Array2<f64>, knots: Array1<f64>,
              p: Array1<f64>, lb: Array1<f64>, ub: Array1<f64>,
              nl: usize, nb: usize, np: usize,
              eq: Option<Vec<(usize, f64, Vec<usize>, Vec<usize>)>>,
              iq: Option<Vec<(usize, f64, Vec<usize>, Vec<usize>)>>)

    -> (Array<f64, Ix1>, Result<(SuccessState, f64), (FailState, f64)>) {

    // Generate Fit object
    let fit = Fit1D::new(x, y, knots, nl, nb, np);

    // Initializing nlopt object
    let n_par: usize = nl + nb*2 + np;
    let mut opt = Nlopt::new(Algorithm::Slsqp, n_par, Fit1D::obj, Target::Minimize, fit);
    //let mut opt = Nlopt::new(Algorithm::Cobyla, n_par, Fit1D::obj, Target::Minimize, fit);
    //let mut opt = Nlopt::new(Algorithm::Isres, n_par, Fit1D::obj, Target::Minimize, fit);

    // Common nlopt setup
    add_constraints(&mut opt, lb, ub, eq, iq);

    // Run the optimization
    let mut p = p.to_vec();
    let out = opt.optimize(&mut p[..]);

    // Reforming array from vector
    let p = Array::from_shape_vec((p.len(), ), p).unwrap();
    
    (p, out)
}

//==============================================================================
pub fn eval_1d(x: Array1<f64>, knots: Array1<f64>, p: Array1<f64>, 
               nl: usize, nb: usize, np: usize) 
    -> Array2<f64> {

    let x = x.into_shared();

    // First lineshape
    let mut lineshape = Lineshape1D::new(x.clone(), nl);
    let p_slice = p.slice(s![0 .. nl]).to_vec(); 
    lineshape.eval(&p_slice);

    // Then baseline
    let mut baseline = Baseline1D::new(x.clone(), knots, nb);
    let p_slice = p.slice(s![nl .. (nl + nb*2)]).to_vec(); 
    baseline.eval(&p_slice);

    // Output is the sum of lineshape and baseline components
    let y = &lineshape.y + &baseline.y;

    // Adding phase (only useful for simulating phase errors)
    let mut phase = Phase1D::new(x.clone(), y, np);
    let p_slice = p.slice(s![(nl + nb*2) .. (nl + nb*2 + np)]).to_vec(); 
    phase.eval(&p_slice);

    phase.y.clone() 
}

//==============================================================================
pub fn baseline_1d(x: Array1<f64>, knots: Array1<f64>, p: Array1<f64>, nb: usize) 
    -> Array2<f64> {

    let x = x.into_shared();

    let mut baseline = Baseline1D::new(x.clone(), knots, nb);
    let p_slice = p.slice(s![0 .. (nb*2)]).to_vec(); 
    baseline.eval(&p_slice);

    baseline.y.clone() 
}

//==============================================================================
pub fn phase_1d(x: Array1<f64>, y: Array2<f64>, p: Array1<f64>, np: usize) 
    -> Array2<f64> {

    let x = x.into_shared();

    let mut phase = Phase1D::new(x.clone(), y, np);
    let p_slice = p.slice(s![0 .. np]).to_vec(); 
    phase.eval(&p_slice);

    phase.y.clone() 
}

//==============================================================================
pub fn fit_2d(x_direct: Array1<f64>, x_indirect: Array1<f64>, y: Array2<f64>, 
              resonances: Array1<usize>, dimensions: Array1<usize>, knots: Array1<f64>,
              p: Array1<f64>, lb: Array1<f64>, ub: Array1<f64>,
              nl: usize, nb: usize, np: usize,
              eq: Option<Vec<(usize, f64, Vec<usize>, Vec<usize>)>>,
              iq: Option<Vec<(usize, f64, Vec<usize>, Vec<usize>)>>)

    -> (Array<f64, Ix1>, Result<(SuccessState, f64), (FailState, f64)>) {

    // Generate Fit object
    let fit = Fit2D::new(x_direct, x_indirect, y, resonances.clone(), dimensions.clone(), knots, nl, nb, np);

    // Initializing nlopt object
    let n_par: usize = nl + nb*4 + np;
    let mut opt = Nlopt::new(Algorithm::Slsqp, n_par, Fit2D::obj, Target::Minimize, fit);
    //let mut opt = Nlopt::new(Algorithm::Cobyla, n_par, Fit1D::obj, Target::Minimize, fit);

    // Common nlopt setup
    add_constraints(&mut opt, lb, ub, eq, iq);

    // Fixing height of direct and indirect components 
    let ranges = Lineshape2D::gen_ranges(resonances, dimensions);

    for i in 0 .. ranges.len() {
        
        let lhs: usize = ranges[i][0].start; 
        let rhs: usize = ranges[i][1].start; 

        let con = Constraint::new(vec![lhs/4], vec![rhs/4], p[lhs+2]/p[rhs+2]);
        opt.add_equality_constraint(Constraint::height, con, 1e-8).unwrap();

    }

    // Run the optimization
    let mut p = p.to_vec();
    let out = opt.optimize(&mut p[..]);

    // Reforming array from vector
    let p = Array::from_shape_vec((p.len(), ), p).unwrap();
    
    (p, out)
}

//==============================================================================
pub fn eval_2d(x_direct: Array1<f64>, x_indirect: Array1<f64>, 
               resonances: Array1<usize>, dimensions: Array1<usize>, knots: Array1<f64>,
               p: Array1<f64>,  nl: usize, nb: usize, np: usize)
    -> Array2<f64> {

    let x_direct = x_direct.into_shared();
    let x_indirect = x_indirect.into_shared();

    // First lineshape
    let mut lineshape = Lineshape2D::new(x_direct.clone(), x_indirect.clone(),
                                         resonances, dimensions);
    let p_slice = p.slice(s![0 .. nl]).to_vec(); 
    lineshape.eval(&p_slice);

    // Then baseline
    let mut baseline = Baseline2D::new(x_direct.clone(), x_indirect.clone(),
                                       knots, nb);
    let p_slice = p.slice(s![nl .. (nl + nb*4)]).to_vec(); 
    baseline.eval(&p_slice);

    // Output is the sum of lineshape and baseline components
    let y = &lineshape.y + &baseline.y;

    // Adding phase (only useful for simulating phase errors)
    let mut phase = Phase2D::new(x_direct.clone(), x_indirect.clone(), y, np);
    let p_slice = p.slice(s![(nl + nb*4) .. (nl + nb*4 + np)]).to_vec(); 
    phase.eval(&p_slice);

    phase.y.clone() 
}

//==============================================================================
pub fn baseline_2d(x_direct: Array1<f64>, x_indirect: Array1<f64>, 
                   knots: Array1<f64>, p: Array1<f64>,  nb: usize)
    -> Array2<f64> {

    let x_direct = x_direct.into_shared();
    let x_indirect = x_indirect.into_shared();

    let mut baseline = Baseline2D::new(x_direct.clone(), x_indirect.clone(),
                                       knots, nb);
    let p_slice = p.slice(s![0 .. (nb*4)]).to_vec(); 
    baseline.eval(&p_slice);

    baseline.y.clone() 
}

//==============================================================================
pub fn phase_2d(x_direct: Array1<f64>, x_indirect: Array1<f64>, y: Array2<f64>,
               p: Array1<f64>,  np: usize)
    -> Array2<f64> {

    let x_direct = x_direct.into_shared();
    let x_indirect = x_indirect.into_shared();

    let mut phase = Phase2D::new(x_direct.clone(), x_indirect.clone(), y, np);
    let p_slice = p.slice(s![0 .. np]).to_vec(); 
    phase.eval(&p_slice);

    phase.y.clone() 
}

