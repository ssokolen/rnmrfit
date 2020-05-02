extern crate rustnmrfit;

use std::f64::consts::PI;
use ndarray::prelude::*;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::BufReader;

use rustnmrfit::fit_1d;

#[derive(Serialize, Deserialize)]
struct JSONVec {
    par: Vec<f64>,
    x: Vec<f64>,
    y_re: Vec<f64>,
    y_im: Vec<f64>,
    basis: Vec<f64>,
    nb: Vec<usize>,
    np: Vec<usize>,
}

struct JSONArray {
    par: Array<f64, Ix1>,
    x: Array<f64, Ix1>,
    y: Array<f64, Ix2>,
    basis: Array<f64, Ix2>,
    nb: usize,
    np: usize,
}

//--------------------------------------
fn read_json(path: &str) -> JSONArray {

    let file = File::open(path).unwrap();
    let reader = BufReader::new(file);

    let json: JSONVec = serde_json::from_reader(reader).unwrap();

    let n_par: usize = json.par.len();
    let n: usize = json.x.len();

    let par: Array<f64, Ix1> = Array::from_shape_vec((n_par, ), json.par).unwrap(); 
    let x: Array<f64, Ix1> = Array::from_shape_vec((n, ), json.x).unwrap(); 

    let y_re: Array<f64, Ix1> = Array::from_shape_vec((n, ), json.y_re).unwrap(); 
    let y_im: Array<f64, Ix1> = Array::from_shape_vec((n, ), json.y_im).unwrap(); 
    
    let mut y: Array<f64, Ix2> = Array::zeros((2, n));

    let mut view = y.slice_mut(s![0, .. ]);
    view.assign(&y_re);
    let mut view = y.slice_mut(s![1, .. ]);
    view.assign(&y_im);

    let basis: Array<f64, Ix2> = Array::from_shape_vec(( n, json.nb[0].max(1) ), json.basis).unwrap(); 

    JSONArray {
        par: par,
        x: x,
        y: y,
        basis: basis,
        nb: json.nb[0],
        np: json.np[0],
    }
}

//--------------------------------------
fn approx_eq(lhs: f64, rhs: f64, tol: f64) -> bool {
    (lhs - rhs).abs() < tol
}

//--------------------------------------
fn test_fit(path: &str, offset: f64, tol: f64) {

    // Reading JSON data
    let json = read_json(path);
    
    let n_par = json.par.len();
    let nb = json.nb;
    let np = json.np;
    let nl = n_par - nb*2 - np;

    let mut lb = Array::zeros((n_par,));
    let mut ub = Array::zeros((n_par,));

    // Defining reasonable lower and upper bounds for parameters
    for i in (0 .. nl).step_by(4) {
        lb[i] = 0.0;
        lb[i+1] = 1e-8;
        lb[i+2] = 0.0;
        lb[i+3] = 0.0;

        ub[i] = 1.0; 
        ub[i+1] = 1.0;
        ub[i+2] = 1.0;
        ub[i+3] = 1.0-1e-8;
    }

    // Defining reasonable lower and upper bounds for baseline
    for i in 0 .. (nb) {
        lb[nl + i] = -1.0;
        lb[nl + nb + i] = -1.0;
        ub[nl + i] = 1.0; 
        ub[nl + nb + i] = 1.0; 
    }

    // Defining reasonable lower and upper bounds for phase
    for i in 0 .. (np) {
        lb[nl + nb*2 + i] = -PI/8.0;
        ub[nl + nb*2 + i] = PI/8.0; 
    }

    // Generating set of initial guesses
    let mut par = &json.par + offset;

    // Initial guess for phase is always 0
    for i in 0 .. (np) {
        par[nl + nb*2 + i] = 0.0;
    }

    // Generating lineshape
    // println!("Before: [{}]", par.iter().fold(String::new(), |acc, &num| acc + &format!("{:.3}", &num) + ", "));
    let (par, result) = fit_1d(json.x, json.y, par, lb, ub, nl, json.nb, json.np, Some(json.basis), None);
    // println!("After: [{}]", par.iter().fold(String::new(), |acc, &num| acc + &format!("{:.3}", &num) + ", "));

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
            assert!(approx_eq(json.par[i], par[i], tol),
                    "parameter {} did not meet tolerance ({} vs {})", 
                    i+1, json.par[i], par[i]);
        }
    }
}

#[test]
fn test_lorentz_fit() {
    test_fit("./tests/data/lorentz_fit.json", 0.05, 2e-2);
}

#[test]
fn test_voigt_fit() {
    test_fit("./tests/data/voigt_fit.json", 0.05, 2e-2);
}

#[test]
fn test_combined_fit() {
    test_fit("./tests/data/combined_fit.json", 0.05, 2e-2);
}

#[test]
fn test_baseline_fit() {
    test_fit("./tests/data/baseline_fit.json", 0.05, 2e-2);
}

#[test]
fn test_phase_fit() {
    test_fit("./tests/data/phase_fit.json", 0.05, 2e-2);
}


