use libc::c_double;
use ndarray::prelude::*;
use std::slice;

use rustnmrfit;

#[no_mangle]
pub extern fn fit_1d(x: *const c_double, y: *const c_double, 
                     p: *mut c_double, lb: *const c_double, ub: *const c_double,  
                     n: i32, nl: i32, nb: i32, np: i32,
                     basis: *const c_double) {

    let n = n as usize;
    let nl = nl as usize;
    let nb = nb as usize;
    let np = np as usize;

    println!("{:?}", vec![n, nl, nb, np]);

    // Basis is special -- if there are no baseline terms, it is ignored
    let basis_in: Option<Array<f64, Ix2>> = if nb > 0 {
        let basis: &[f64] = unsafe { slice::from_raw_parts(basis, n*nb) };   
        Some(Array::from_shape_vec((n, nb), basis.to_vec()).unwrap())
    } else {
        None
    };

    // Converting rest into slices
    let x: &[f64] = unsafe { slice::from_raw_parts(x, n) };
    let y: &[f64] = unsafe { slice::from_raw_parts(y, n*2) };
    let p: &mut[f64] = unsafe { slice::from_raw_parts_mut(p, nl + nb*2 + np) };   
    let lb: &[f64] = unsafe { slice::from_raw_parts(lb, nl + nb*2 + np) };   
    let ub: &[f64] = unsafe { slice::from_raw_parts(ub, nl + nb*2 + np) };   

    // Converting rest into arrays
    let x = Array::from_shape_vec((n,), x.to_vec()).unwrap();
    let y = Array::from_shape_vec((2, n), y.to_vec()).unwrap();
    let p_in = Array::from_shape_vec((nl + nb*2 + np,), p.to_vec()).unwrap();
    let lb = Array::from_shape_vec((nl + nb*2 + np,), lb.to_vec()).unwrap();
    let ub = Array::from_shape_vec((nl + nb*2 + np,), ub.to_vec()).unwrap();

    // Calling fit
    let (p_out, _) = rustnmrfit::fit_1d(x, y, p_in, lb, ub, nl, nb, np, basis_in);

    // Copying over new parameters
    for i in 0 .. (nl + nb*2 + np) {
        p[i] = p_out[i];
    }
}

#[no_mangle]
pub extern fn eval_1d(x: *const c_double, y: *mut c_double, p: *const c_double, 
                      n: i32, nl: i32, nb: i32, np: i32,
                      basis: *const c_double) {

    let n = n as usize;
    let nl = nl as usize;
    let nb = nb as usize;
    let np = np as usize;

    // Basis is special -- if there are no baseline terms, it is ignored
    let basis_in: Option<Array<f64, Ix2>> = if nb > 0 {
        let basis: &[f64] = unsafe { slice::from_raw_parts(basis, n*nb) };   
        Some(Array::from_shape_vec((n, nb), basis.to_vec()).unwrap())
    } else {
        None
    };

    // Converting rest into slices
    let x: &[f64] = unsafe { slice::from_raw_parts(x, n) };
    let y: &mut[f64] = unsafe { slice::from_raw_parts_mut(y, n*2) };
    let p: &[f64] = unsafe { slice::from_raw_parts(p, nl + nb*2 + np) };   

    // Converting rest into arrays
    let x = Array::from_shape_vec((n,), x.to_vec()).unwrap();
    let p = Array::from_shape_vec((nl + nb*2 + np,), p.to_vec()).unwrap();

    // Calling eval
    let y_out = rustnmrfit::eval_1d(x, p, nl, nb, np, basis_in);

    // Copying over y values
    for i in 0 .. (n) {
        y[i] = y_out[[0, i]];
        y[i + n] = y_out[[1, i]];
    }
}
