/*
extern crate rnmrfit;

#[macro_use(concatenate)]
extern crate ndarray;

use std::f64::consts::PI;
use ndarray::prelude::*;
use ndarray_rand::RandomExt;
use ndarray_rand::rand::SeedableRng;
use ndarray_rand::rand_distr::{Normal, Distribution};
use rand_isaac::isaac64::Isaac64Rng;

use rnmrfit::{eval_1d, fit_1d, eval_2d, fit_2d};

//--------------------------------------
fn test_fit_1d(nl: usize, nb: usize, np: usize, offset: f64, tol: f64) {

    // Setting up ideal fit
    let x: Array1<f64> = Array::linspace(0.0, 1.0, 50);
    let p = [&(vec![0.7, 0.05, 0.8, 0.0, 
                    0.3, 0.05, 0.6, 0.5])[0 .. nl],
             &(vec![0.05, 0.1, 0.1, 0.05][0 .. nb]),
             &(vec![0.05, 0.1, 0.1, 0.05][0 .. nb]),
             &(vec![0.05, 0.05])[0 .. np]].concat();

    // Basic boundary knots
    let knots: Array1<f64> = Array::from_shape_vec(
        (2,), vec![0.0, 1.0]
    ).unwrap();

    // Initialize work arrays
    let mut y: Array2<f64> = Array::zeros((n, 2));
    let mut dy: Array3<f64> = Array::zeros((nl, n, 2));

    // Generate peaks
    let y_slice = y.slice_mut(s![.., ..]);
    let dy_slice = dy.slice_mut(s![.., .., ..]);

    let mut lineshape = Lineshape1D::new(x.clone(), np);
    lineshape.eval(&p[..nl], y_slice, dy_slice);

    // Adding some noise for a bit of realism
    let seed = 1111;
    let mut rng = Isaac64Rng::seed_from_u64(seed);
    let noise = Array::random_using(
        (n, 2), Normal::new(0., 0.01).unwrap(), &mut rng
    );
    y += &noise;

    let mut lb = Array::zeros((p.len(),));
    let mut ub = Array::zeros((p.len(),));

    // Defining reasonable lower and upper bounds for parameters
    for i in (0 .. nl).step_by(4) {
        lb[i] = p[i] - 0.05;
        lb[i+1] = 1e-8;
        lb[i+2] = 0.0;
        lb[i+3] = 0.0;

        ub[i] = p[i] + 0.05; 
        ub[i+1] = 1.0;
        ub[i+2] = 1.0;
        ub[i+3] = 0.9;
    }

    // Defining reasonable lower and upper bounds for baseline
    for i in 0 .. nb {
        lb[nl + i] = -1.0;
        lb[nl + nb + i] = -1.0;
        ub[nl + i] = 1.0; 
        ub[nl + nb + i] = 1.0; 
    }

    // Defining reasonable lower and upper bounds for phase
    for i in 0 .. np {
        lb[nl + nb*2 + i] = -PI/8.0;
        ub[nl + nb*2 + i] = PI/8.0; 
    }

    // Generating set of initial guesses (baseline and phase guess is 0)
    let mut p0 = Array::zeros(p.raw_dim());
    for i in 0 .. nl {
        p0[i] = p[i] + offset;
    }

    // Fitting
    let (par, result) = fit_1d(x, y, knots, p0.clone(), lb, ub, nl, nb, np, None, None, 0, 1e-4, 0.0);

    // Testing optimization state
    let success = match result {
        Ok(_) => true,
        Err(_) => false,
    };
    
    assert!(success, "optimization failed with: {:?}", result);

    // Testing deviation from parameters
    for i in 0 .. nl {
        // Fractions are ignored as they are quite sensitive
        if (i + 1) % 4 != 0 {
            assert!( (par[i] - p[i]).abs() < tol,
                    "parameter {} did not meet tolerance ({} vs {})", 
                    i, par[i], p[i]);
        }
    }
}

#[test]
fn test_single_peak_fit_1d() {
    test_fit_1d(4, 0, 0, 0.05, 0.02);
}

#[test]
fn test_multiple_peak_fit_1d() {
    test_fit_1d(8, 0, 0, 0.05, 0.02);
}

#[test]
fn test_baseline_fit_1d() {
    test_fit_1d(8, 4, 0, 0.05, 0.02);
}

#[test]
fn test_phase_fit_1d() {
    test_fit_1d(8, 4, 2, 0.05, 0.02);
}

//--------------------------------------
fn test_fit_2d(nl: usize, nb: usize, np: usize, offset: f64, tol: f64) {

    // Setting up ideal fit
    let x: Array1<f64> = Array::linspace(0.0, 1.0, 50);
    let p = [&(vec![0.7, 0.05, 0.8, 0.0, 
                    0.7, 0.05, 0.6, 0.5,
                    0.3, 0.05, 0.8, 0.0, 
                    0.3, 0.05, 0.6, 0.5])[0 .. nl],
             &(vec![0.05, 0.1, 0.1, 0.05][0 .. nb]),
             &(vec![0.05, 0.1, 0.1, 0.05][0 .. nb]),
             &(vec![0.05, 0.1, 0.1, 0.05][0 .. nb]),
             &(vec![0.05, 0.1, 0.1, 0.05][0 .. nb]),
             &(vec![0.05, 0.05, 0.05])[0 .. np]].concat();

    let resonances: Array1<usize> = concatenate![Axis(0), Array::zeros((8,)), Array::from_elem((8,), 1)];
    let dimensions: Array1<usize> = concatenate![Axis(0), Array::zeros((4,)), Array::from_elem((4,), 1),
                                                    Array::zeros((4,)), Array::from_elem((4,), 1)];

    let resonances = resonances.slice(s![0 .. nl]).to_owned();
    let dimensions = dimensions.slice(s![0 .. nl]).to_owned();

    // Basic boundary knots
    let knots: Array1<f64> = Array::from_shape_vec((2,), vec![0.0, 1.0]).unwrap();

    // Building grid from x
    let n = x.len();
    let mut x1: Array1<f64> = Array::zeros((n*n,));
    let mut x2 = x1.clone();

    for i in 0 .. n {
        for j in 0 .. n {
            x1[i+n*j] = x[j];
            x2[i+n*j] = x[i];
        }
    }

    // Generating ideal y
    let p = Array::from_shape_vec((p.len(),), p).unwrap();
    let mut y = eval_2d(x1.clone(), x2.clone(), resonances.clone(), dimensions.clone(), 
                        knots.clone(), p.clone(), nl, nb, np); 

    // Adding some noise for a bit of realism
    let seed = 1111;
    let mut rng = Isaac64Rng::seed_from_u64(seed);
    let noise = Array::random_using(y.raw_dim(), Normal::new(0., 0.01).unwrap(), &mut rng);
    y += &noise;

    let mut lb = Array::zeros((p.len(),));
    let mut ub = Array::zeros((p.len(),));

    // Defining reasonable lower and upper bounds for parameters
    for i in (0 .. nl).step_by(4) {
        lb[i] = p[i] - 0.05;
        lb[i+1] = 1e-8;
        lb[i+2] = 0.0;
        lb[i+3] = 0.0;

        ub[i] = p[i] + 0.05; 
        ub[i+1] = 1.0;
        ub[i+2] = 1.0;
        ub[i+3] = 0.9;
    }

    // Defining reasonable lower and upper bounds for baseline
    for i in 0 .. nb {
        for j in 0 .. 4 {
            lb[nl + nb*j + i] = -1.0;
            ub[nl + nb*j + i] = 1.0; 
        }
    }

    // Defining reasonable lower and upper bounds for phase
    for i in 0 .. np {
        lb[nl + nb*4 + i] = -PI/8.0;
        ub[nl + nb*4 + i] = PI/8.0; 
    }

    // Generating set of initial guesses (baseline and phase guess is 0)
    let mut p0 = Array::zeros(p.raw_dim());
    for i in 0 .. nl {
        p0[i] = p[i] + offset;
    }

    // Generating lineshape
    let (par, result) = fit_2d(x1, x2, y, resonances, dimensions, knots, 
                               p0.clone(), lb, ub, nl, nb, np, None, None);

    // Testing optimization state
    let success = match result {
        Ok(_) => true,
        Err(_) => false,
    };
    
    assert!(success, "optimization failed with: {:?}", result);

    // Testing deviation from parameters
    for i in 0 .. nl {
        // Fractions are ignored as they are quite sensitive
        if (i + 1) % 4 != 0 {
            assert!( (par[i] - p[i]).abs() < tol,
                    "parameter {} did not meet tolerance ({} vs {})", 
                    i, par[i], p[i]);
        }
    }
}

#[test]
fn test_single_peak_fit_2d() {
    test_fit_2d(8, 0, 0, 0.05, 0.02);
}

#[test]
fn test_multiple_peak_fit_2d() {
    test_fit_2d(16, 0, 0, 0.05, 0.02);
}

#[test]
fn test_baseline_fit_2d() {
    test_fit_2d(16, 4, 0, 0.05, 0.02);
}

#[test]
fn test_phase_fit_2d() {
    test_fit_2d(16, 4, 2, 0.05, 0.02);
}
*/
