use libc::c_double;
use ndarray::prelude::*;
use std::slice;

use rustnmrfit;

fn parse_constraints(eq: &[f64]) 
    -> Vec<(usize, f64, Vec<usize>, Vec<usize>)> {
    
    let mut eq_in: Vec<(usize, f64, Vec<usize>, Vec<usize>)> = Vec::new();

    let mut flag: usize;
    let mut offset: f64;
    let mut lhs: Vec<usize> = Vec::new();
    let mut rhs: Vec<usize> = Vec::new();

    let neq = eq.len();
    
    // Starting parsing
    if neq > 2 {
        flag = eq[0].round() as usize;
        offset = eq[1] as f64;
    } else {
        panic!("Constraints too short to be parsed, aborting.")
    }

    let mut i: usize = 2;
    while i < neq {
    
        // NAN signifies new set of constraints
        if eq[i].is_nan() {

            // Append previous set
            eq_in.push( (flag, offset, lhs.clone(), rhs.clone()) );
            lhs = Vec::new();
            rhs = Vec::new();

            flag = eq[i+1].round() as usize;
            offset = eq[i+2] as f64;
            i += 2;
            
        } else if eq[i] > 0.0 {
            
            lhs.push( (eq[i].round() as usize) - 1 );
            
        } else if eq[i] < 0.0 {
        
            rhs.push( (-eq[i].round() as usize) - 1 );
            
        }
        
        i += 1;
    }
    
    // Pushing last created constraint
    eq_in.push( (flag, offset, lhs.clone(), rhs.clone()) );
    
    eq_in
}

#[no_mangle]
pub extern fn fit_1d(x: *const c_double, y: *const c_double, knots: *const c_double,
                     p: *mut c_double, lb: *const c_double, ub: *const c_double,  
                     n: i32, nl: i32, nb: i32, np: i32, nk: i32,
                     eq: *const c_double, iq: *const c_double,
                     neq: i32, niq: i32) {

    let n = n as usize;
    let nl = nl as usize;
    let nb = nb as usize;
    let np = np as usize;
    let nk = nk as usize;
    let neq = neq as usize;
    let niq = niq as usize;

    // Parsing equality constraint tuples from a flat list
    let eq_in: Option<Vec<(usize, f64, Vec<usize>, Vec<usize>)>> = if neq > 0 {
        let eq: &[f64] = unsafe { slice::from_raw_parts(eq, neq) };
        Some(parse_constraints(eq))
    } else {
        None
    };

    // Parsing inequality constraint tuples from a flat list
    let iq_in: Option<Vec<(usize, f64, Vec<usize>, Vec<usize>)>> = if niq > 0 {
        let iq: &[f64] = unsafe { slice::from_raw_parts(iq, niq) };
        Some(parse_constraints(iq))
    } else {
        None
    };

    // Converting rest into slices
    let x: &[f64] = unsafe { slice::from_raw_parts(x, n) };
    let y: &[f64] = unsafe { slice::from_raw_parts(y, n*2) };
    let knots: &[f64] = unsafe { slice::from_raw_parts(knots, nk) };
    let p: &mut[f64] = unsafe { slice::from_raw_parts_mut(p, nl + nb*2 + np) };   
    let lb: &[f64] = unsafe { slice::from_raw_parts(lb, nl + nb*2 + np) };   
    let ub: &[f64] = unsafe { slice::from_raw_parts(ub, nl + nb*2 + np) };   

    // Converting rest into arrays
    let x = Array::from_shape_vec((n,), x.to_vec()).unwrap();
    let y = Array::from_shape_vec((2, n), y.to_vec()).unwrap();
    let knots = Array::from_shape_vec((nk,), knots.to_vec()).unwrap();
    let p_in = Array::from_shape_vec((nl + nb*2 + np,), p.to_vec()).unwrap();
    let lb = Array::from_shape_vec((nl + nb*2 + np,), lb.to_vec()).unwrap();
    let ub = Array::from_shape_vec((nl + nb*2 + np,), ub.to_vec()).unwrap();

    // Calling fit
    let (p_out, _) = rustnmrfit::fit_1d(x, y, knots, p_in, lb, ub, nl, nb, np, eq_in, iq_in);

    // Copying over new parameters
    for i in 0 .. (nl + nb*2 + np) {
        p[i] = p_out[i];
    }
}

#[no_mangle]
pub extern fn eval_1d(x: *const c_double, y: *mut c_double, knots: *const c_double,
                      p: *const c_double, n: i32, nl: i32, nb: i32, np: i32, nk: i32) {

    let n = n as usize;
    let nl = nl as usize;
    let nb = nb as usize;
    let np = np as usize;
    let nk = nk as usize;

    // Converting rest into slices
    let x: &[f64] = unsafe { slice::from_raw_parts(x, n) };
    let y: &mut[f64] = unsafe { slice::from_raw_parts_mut(y, n*2) };
    let p: &[f64] = unsafe { slice::from_raw_parts(p, nl + nb*2 + np) };   
    let knots: &[f64] = unsafe { slice::from_raw_parts(knots, nk) };

    // Converting rest into arrays
    let x = Array::from_shape_vec((n,), x.to_vec()).unwrap();
    let p = Array::from_shape_vec((nl + nb*2 + np,), p.to_vec()).unwrap();
    let knots = Array::from_shape_vec((nk,), knots.to_vec()).unwrap();

    // Calling eval
    let y_out = rustnmrfit::eval_1d(x, knots, p, nl, nb, np);

    // Copying over y values
    for i in 0 .. n {
        y[i] = y_out[[0, i]];
        y[i + n] = y_out[[1, i]];
    }
}
