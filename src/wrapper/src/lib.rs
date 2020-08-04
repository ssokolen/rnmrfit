use libc::c_double;
use ndarray::prelude::*;
use std::slice;
use nlopt::{SuccessState, FailState};

use rnmrfit;

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
                     neq: i32, niq: i32, alg: i32, xtr: f64, mxt: f64,
                     out: *mut i32) {

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
    let out: &mut[i32] = unsafe { slice::from_raw_parts_mut(out, 1) };

    // Converting rest into arrays
    let x = Array::from_shape_vec((n,), x.to_vec()).unwrap();
    let y = Array::from_shape_vec((2, n), y.to_vec()).unwrap();
    let knots = Array::from_shape_vec((nk,), knots.to_vec()).unwrap();
    let p_in = Array::from_shape_vec((nl + nb*2 + np,), p.to_vec()).unwrap();
    let lb = Array::from_shape_vec((nl + nb*2 + np,), lb.to_vec()).unwrap();
    let ub = Array::from_shape_vec((nl + nb*2 + np,), ub.to_vec()).unwrap();

    // Calling fit
    let (p_out, result) = rnmrfit::fit_1d(x, y, knots, p_in, lb, ub, nl, nb, np, eq_in, iq_in,
                                          alg, xtr, mxt);

    // Copying over new parameters
    for i in 0 .. (nl + nb*2 + np) {
        p[i] = p_out[i];
    }

    // Encoding status
    out[0] = match result {
        Ok((code, _value)) => { 
            match code {
                SuccessState::Success => 1,
                SuccessState::StopValReached => 2,
                SuccessState::FtolReached => 3,
                SuccessState::XtolReached => 4,
                SuccessState::MaxEvalReached => 5,
                SuccessState::MaxTimeReached => 6,
            }
        },
        Err((code, _value)) => { 
            match code {
                FailState::Failure => -1,
                FailState::InvalidArgs => -2,
                FailState::OutOfMemory => -3,
                FailState::RoundoffLimited => -4,
                FailState::ForcedStop => -5,
            }
        }
    };
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
    let y_out = rnmrfit::eval_1d(x, knots, p, nl, nb, np);

    // Copying over y values
    for i in 0 .. n {
        y[i] = y_out[[0, i]];
        y[i + n] = y_out[[1, i]];
    }
}

#[no_mangle]
pub extern fn baseline_1d(x: *const c_double, y: *mut c_double, knots: *const c_double,
                          p: *const c_double, n: i32, nb: i32, nk: i32) {

    let n = n as usize;
    let nb = nb as usize;
    let nk = nk as usize;

    // Converting rest into slices
    let x: &[f64] = unsafe { slice::from_raw_parts(x, n) };
    let y: &mut[f64] = unsafe { slice::from_raw_parts_mut(y, n*2) };
    let p: &[f64] = unsafe { slice::from_raw_parts(p, nb*2) };   
    let knots: &[f64] = unsafe { slice::from_raw_parts(knots, nk) };

    // Converting rest into arrays
    let x = Array::from_shape_vec((n,), x.to_vec()).unwrap();
    let p = Array::from_shape_vec((nb*2,), p.to_vec()).unwrap();
    let knots = Array::from_shape_vec((nk,), knots.to_vec()).unwrap();

    // Calling baseline
    let y_out = rnmrfit::baseline_1d(x, knots, p, nb);

    // Copying over y values
    for i in 0 .. n {
        y[i] = y_out[[0, i]];
        y[i + n] = y_out[[1, i]];
    }
}

#[no_mangle]
pub extern fn phase_1d(x: *const c_double, y: *mut c_double,
                       p: *const c_double, n: i32, np: i32) {

    let n = n as usize;
    let np = np as usize;

    // Converting rest into slices
    let x: &[f64] = unsafe { slice::from_raw_parts(x, n) };
    let y: &mut[f64] = unsafe { slice::from_raw_parts_mut(y, n*2) };
    let p: &[f64] = unsafe { slice::from_raw_parts(p, np) };   

    // Converting rest into arrays
    let x = Array::from_shape_vec((n,), x.to_vec()).unwrap();
    let y_array = Array::from_shape_vec((2,n), y.to_vec()).unwrap();
    let p = Array::from_shape_vec((np,), p.to_vec()).unwrap();

    // Calling baseline
    let y_out = rnmrfit::phase_1d(x, y_array, p, np);

    // Copying over y values
    for i in 0 .. n {
        y[i] = y_out[[0, i]];
        y[i + n] = y_out[[1, i]];
    }
}


#[no_mangle]
pub extern fn fit_2d(x_direct: *const c_double, x_indirect: *const c_double, y: *const c_double, 
                     resonances: *const i32, dimensions: *const i32, knots: *const c_double,
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
    let x_direct: &[f64] = unsafe { slice::from_raw_parts(x_direct, n) };
    let x_indirect: &[f64] = unsafe { slice::from_raw_parts(x_indirect, n) };
    let y: &[f64] = unsafe { slice::from_raw_parts(y, n*4) };
    let resonances: &[i32] = unsafe { slice::from_raw_parts(resonances, nl) };
    let dimensions: &[i32] = unsafe { slice::from_raw_parts(dimensions, nl) };
    let knots: &[f64] = unsafe { slice::from_raw_parts(knots, nk) };
    let p: &mut[f64] = unsafe { slice::from_raw_parts_mut(p, nl + nb*4 + np) };   
    let lb: &[f64] = unsafe { slice::from_raw_parts(lb, nl + nb*4 + np) };   
    let ub: &[f64] = unsafe { slice::from_raw_parts(ub, nl + nb*4 + np) };   

    // Converting rest into arrays
    let x_direct = Array::from_shape_vec((n,), x_direct.to_vec()).unwrap();
    let x_indirect = Array::from_shape_vec((n,), x_indirect.to_vec()).unwrap();
    let y = Array::from_shape_vec((4, n), y.to_vec()).unwrap();
    let resonances = Array::from_shape_vec((nl,), resonances.to_vec()).unwrap();
    let resonances = resonances.mapv(|elem| elem as usize);
    let dimensions = Array::from_shape_vec((nl,), dimensions.to_vec()).unwrap();
    let dimensions = dimensions.mapv(|elem| elem as usize);
    let knots = Array::from_shape_vec((nk,), knots.to_vec()).unwrap();
    let p_in = Array::from_shape_vec((nl + nb*4 + np,), p.to_vec()).unwrap();
    let lb = Array::from_shape_vec((nl + nb*4 + np,), lb.to_vec()).unwrap();
    let ub = Array::from_shape_vec((nl + nb*4 + np,), ub.to_vec()).unwrap();

    // Calling fit
    let (p_out, _) = rnmrfit::fit_2d(x_direct, x_indirect, y, resonances, dimensions, 
                                     knots, p_in, lb, ub, nl, nb, np, eq_in, iq_in);

    // Copying over new parameters
    for i in 0 .. (nl + nb*4 + np) {
        p[i] = p_out[i];
    }
}

#[no_mangle]
pub extern fn eval_2d(x_direct: *const c_double, x_indirect: *const c_double, y: *mut c_double, 
                      resonances: *const i32, dimensions: *const i32, knots: *const c_double,
                      p: *const c_double, n: i32, nl: i32, nb: i32, np: i32, nk: i32) {

    let n = n as usize;
    let nl = nl as usize;
    let nb = nb as usize;
    let np = np as usize;
    let nk = nk as usize;

    // Converting rest into slices
    let x_direct: &[f64] = unsafe { slice::from_raw_parts(x_direct, n) };
    let x_indirect: &[f64] = unsafe { slice::from_raw_parts(x_indirect, n) };
    let y: &mut[f64] = unsafe { slice::from_raw_parts_mut(y, n*4) };
    let resonances: &[i32] = unsafe { slice::from_raw_parts(resonances, nl) };
    let dimensions: &[i32] = unsafe { slice::from_raw_parts(dimensions, nl) };
    let p: &[f64] = unsafe { slice::from_raw_parts(p, nl + nb*4 + np) };   
    let knots: &[f64] = unsafe { slice::from_raw_parts(knots, nk) };

    // Converting rest into arrays
    let x_direct = Array::from_shape_vec((n,), x_direct.to_vec()).unwrap();
    let x_indirect = Array::from_shape_vec((n,), x_indirect.to_vec()).unwrap();
    let resonances = Array::from_shape_vec((nl,), resonances.to_vec()).unwrap();
    let resonances = resonances.mapv(|elem| elem as usize);
    let dimensions = Array::from_shape_vec((nl,), dimensions.to_vec()).unwrap();
    let dimensions = dimensions.mapv(|elem| elem as usize);
    let p = Array::from_shape_vec((nl + nb*4+ np,), p.to_vec()).unwrap();
    let knots = Array::from_shape_vec((nk,), knots.to_vec()).unwrap();

    // Calling eval
    let y_out = rnmrfit::eval_2d(x_direct, x_indirect, resonances, dimensions, 
                                 knots, p, nl, nb, np);

    // Copying over y values
    for i in 0 .. n {
        for j in 0 .. 4 {
            y[i + n*j] = y_out[[j, i]];
        }
    }
}

#[no_mangle]
pub extern fn baseline_2d(x_direct: *const c_double, x_indirect: *const c_double, y: *mut c_double, 
                          knots: *const c_double, p: *const c_double, n: i32, nb: i32, nk: i32) {

    let n = n as usize;
    let nb = nb as usize;
    let nk = nk as usize;

    // Converting rest into slices
    let x_direct: &[f64] = unsafe { slice::from_raw_parts(x_direct, n) };
    let x_indirect: &[f64] = unsafe { slice::from_raw_parts(x_indirect, n) };
    let y: &mut[f64] = unsafe { slice::from_raw_parts_mut(y, n*4) };
    let p: &[f64] = unsafe { slice::from_raw_parts(p, nb*4) };   
    let knots: &[f64] = unsafe { slice::from_raw_parts(knots, nk) };

    // Converting rest into arrays
    let x_direct = Array::from_shape_vec((n,), x_direct.to_vec()).unwrap();
    let x_indirect = Array::from_shape_vec((n,), x_indirect.to_vec()).unwrap();
    let p = Array::from_shape_vec((nb*4,), p.to_vec()).unwrap();
    let knots = Array::from_shape_vec((nk,), knots.to_vec()).unwrap();

    // Calling eval
    let y_out = rnmrfit::baseline_2d(x_direct, x_indirect, knots, p, nb);

    // Copying over y values
    for i in 0 .. n {
        for j in 0 .. 4 {
            y[i + n*j] = y_out[[j, i]];
        }
    }
}

#[no_mangle]
pub extern fn phase_2d(x_direct: *const c_double, x_indirect: *const c_double, y: *mut c_double, 
                       p: *const c_double, n: i32, np: i32) {

    let n = n as usize;
    let np = np as usize;

    // Converting rest into slices
    let x_direct: &[f64] = unsafe { slice::from_raw_parts(x_direct, n) };
    let x_indirect: &[f64] = unsafe { slice::from_raw_parts(x_indirect, n) };
    let y: &mut[f64] = unsafe { slice::from_raw_parts_mut(y, n*4) };
    let p: &[f64] = unsafe { slice::from_raw_parts(p, np) };   

    // Converting rest into arrays
    let x_direct = Array::from_shape_vec((n,), x_direct.to_vec()).unwrap();
    let x_indirect = Array::from_shape_vec((n,), x_indirect.to_vec()).unwrap();
    let y_array = Array::from_shape_vec((4,n), y.to_vec()).unwrap();
    let p = Array::from_shape_vec((np,), p.to_vec()).unwrap();

    // Calling eval
    let y_out = rnmrfit::phase_2d(x_direct, x_indirect, y_array, p, np);

    // Copying over y values
    for i in 0 .. n {
        for j in 0 .. 4 {
            y[i + n*j] = y_out[[j, i]];
        }
    }
}
