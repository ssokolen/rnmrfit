use std::f64::consts::PI;
use std::f64::consts::SQRT_2;
use rustfft::num_complex::Complex;

use crate::faddeeva;

// Some constants
#[allow(dead_code)]
const SQRT_PI: f64 = 1.772453850905515881919427556567825376987457275390625_f64;

//=============================================================================
// Constraint

//-----------------------------------------------------------------------------
#[allow(dead_code)]
#[derive(Debug)]
pub struct Constraint {

    // Two sets of peak indexes
    lhs: Vec<usize>,
    rhs: Vec<usize>,

    // Offset between peak indexes, e.g. lhs - rhs = offset or lhs/rhs = offset
    offset: f64,

}

//-----------------------------------------------------------------------------
#[allow(dead_code)]
impl Constraint {

    pub fn new(lhs: Vec<usize>, rhs: Vec<usize>, offset: f64) 
        -> Constraint {

        Constraint {
            lhs: lhs,
            rhs: rhs,
            offset: offset,
        }
    }
}

//-----------------------------------------------------------------------------
#[allow(dead_code)]
impl Constraint {

    //--------------------------------------
    pub fn position(p: &[f64], grad: Option<&mut [f64]>, obj: &mut Constraint) -> f64 {

        // The typical scenario is a simple position_1 - position_2 = offset
        // But this is generalized to sum(position_i) - sum(position_j) = offset

        let mut value: f64 = 0.0;
        let n_lhs = obj.lhs.len();
        let n_rhs = obj.rhs.len();

        // i -> i*4 to from peak index to position parameter
        for i in 0 .. n_lhs { value += p[ obj.lhs[i]*4 ]; }
        for i in 0 .. n_rhs { value -= p[ obj.rhs[i]*4 ]; }

        // Gradients are simple +1/-1
        if grad.is_some() {
            let grad = grad.unwrap();

            for i in 0 .. grad.len() { grad[i] = 0.0; }
            for i in 0 .. n_lhs { grad[ obj.lhs[i]*4 ] = 1.0; }
            for i in 0 .. n_rhs { grad[ obj.rhs[i]*4 ] = -1.0; }
        }

        value - obj.offset
    }

    //--------------------------------------
    pub fn width(p: &[f64], grad: Option<&mut [f64]>, obj: &mut Constraint) -> f64 {

        // The typical scenario is a simple width_1/width_2 = offset
        // But this is generalized to sum(width_i)/sum(width_j) = offset

        let mut num: f64 = 0.0;
        let mut denom: f64 = 0.0;

        let n_lhs = obj.lhs.len();
        let n_rhs = obj.rhs.len();

        // i -> i*4 + 1 to from peak index to width parameter
        for i in 0 .. n_lhs { num += p[ obj.lhs[i]*4 + 1 ]; }
        for i in 0 .. n_rhs { denom += p[ obj.rhs[i]*4 + 1 ]; }

        let value: f64 = num/denom;

        // Gradients depend on whether term is in numerator or denominator
        if grad.is_some() {
            let grad = grad.unwrap();

            let grad_num: f64 = 1.0/denom;
            let grad_denom: f64 = -value/denom;

            for i in 0 .. grad.len() { grad[i] = 0.0; }
            for i in 0 .. n_lhs { grad[ obj.lhs[i]*4 + 1] = grad_num; }
            for i in 0 .. n_rhs { grad[ obj.rhs[i]*4 + 1] = grad_denom; }
        }

        value - obj.offset
    }

    //--------------------------------------
    pub fn height(p: &[f64], grad: Option<&mut [f64]>, obj: &mut Constraint) -> f64 {

        // The typical scenario is a simple height_1/height_2 = offset
        // But this is generalized to sum(height_i)/sum(height_j) = offset

        let mut num: f64 = 0.0;
        let mut denom: f64 = 0.0;

        let n_lhs = obj.lhs.len();
        let n_rhs = obj.rhs.len();

        // i -> i*4 + 2 to from peak index to height parameter
        for i in 0 .. n_lhs { num += p[ obj.lhs[i]*4 + 2 ]; }
        for i in 0 .. n_rhs { denom += p[ obj.rhs[i]*4 + 2 ]; }

        let value: f64 = num/denom;

        // Gradients depend on whether term is in numerator or denominator
        if grad.is_some() {
            let grad = grad.unwrap();

            let grad_num: f64 = 1.0/denom;
            let grad_denom: f64 = -value/denom;

            for i in 0 .. grad.len() { grad[i] = 0.0; }
            for i in 0 .. n_lhs { grad[ obj.lhs[i]*4 + 2] = grad_num; }
            for i in 0 .. n_rhs { grad[ obj.rhs[i]*4 + 2] = grad_denom; }
        }

        value - obj.offset
    }

    //--------------------------------------
    pub fn fraction(p: &[f64], grad: Option<&mut [f64]>, obj: &mut Constraint) -> f64 {

        // The typical scenario is a simple fraction_1 - fraction_2 = offset
        // But this is generalized to sum(fraction_i) - sum(fraction_j) = offset

        let mut value: f64 = 0.0;
        let n_lhs = obj.lhs.len();
        let n_rhs = obj.rhs.len();

        // i -> i*4 + 3 to from peak index to fraction parameter
        for i in 0 .. n_lhs { value += p[ obj.lhs[i]*4 + 3]; }
        for i in 0 .. n_rhs { value -= p[ obj.rhs[i]*4 + 3]; }

        // Gradients are simple +1/-1
        if grad.is_some() {
            let grad = grad.unwrap();

            for i in 0 .. grad.len() { grad[i] = 0.0; }
            for i in 0 .. n_lhs { grad[ obj.lhs[i]*4 + 3] = 1.0; }
            for i in 0 .. n_rhs { grad[ obj.rhs[i]*4 + 3] = -1.0; }
        }

        value - obj.offset
    }

    //--------------------------------------
    pub fn area(p: &[f64], grad: Option<&mut [f64]>, obj: &mut Constraint) -> f64 {

        // General area function
        let f = | i | {
            let w = p[i*4 + 1];
            let h = p[i*4 + 2];
            let f = p[i*4 + 3];
            let wg = w * f/(1.0 - f);

            if f < 1e-8 { 
              
                PI * w * h 
            
            } else if f < (1.0 - 1e-8) { 
              let z: Complex<f64> = Complex::new(0.0, w/(SQRT_2 * wg));
              let y = faddeeva::w( z, 1e-6 );

              SQRT_2 * SQRT_PI * wg * h / y.re
            
            } else if f <= 1.0 {
              
                SQRT_2 * SQRT_PI * wg * h
            
            } else {
                panic!("Fraction Gauss fell out of 0-1 range, aborting.")
            }        
        };



        // Unpacking parameters

        // The typical scenario is a simple area_1/area_2 = offset
        // But this is generalized to sum(area_i)/sum(area_j) = offset

        let mut num: f64 = 0.0;
        let mut denom: f64 = 0.0;

        let n_lhs = obj.lhs.len();
        let n_rhs = obj.rhs.len();

        for i in 0 .. n_lhs { num += f( obj.lhs[i] ); }
        for i in 0 .. n_rhs { denom += f( obj.rhs[i] ); }

        let value: f64 = num/denom;

        // Gradients depend on whether term is in numerator or denominator
        // The df function is used to calculate one component of the chain
        // which is then scaled by the other fixed component
        if grad.is_some() {
            let grad = grad.unwrap();

            // General gradient function
            let df = | i, grad: &mut[f64] | {

                let w = p[i*4 + 1];
                let h = p[i*4 + 2];
                let f = p[i*4 + 3];
                let wg = w * f/(1.0 - f);

                if f < 1e-8 { 
                  
                    grad[i*4 + 1] = PI * h;
                    grad[i*4 + 2] = PI * w;
                
                } else if f < (1.0 - 1e-8) { 
                    let z: Complex<f64> = Complex::new(0.0, w/(SQRT_2 * wg));
                    let y = faddeeva::w( z, 1e-6 );

                    let a = SQRT_2 * SQRT_PI / y.re;
                    let b = a / PI;

                    grad[i*4 + 1] = a * h * f/(1.0 - f);
                    grad[i*4 + 2] = a * wg;
                    grad[i*4 + 3] = a * w * h * 
                        ( 1.0/( (1.0-f)*(1.0-f) ) + 1.0/( f*f ) - b/( f*(1.0 - f) ) );
                
                } else if f <= 1.0 {
                  
                    let a = SQRT_2 * SQRT_PI;
                    grad[i*4 + 1] = a * h * f/(1.0 - f);
                    grad[i*4 + 2] = a * wg;
                    grad[i*4 + 3] = a * h * w /( (1.0 - f)*(1.0 - f) );
                
                } else {
                    panic!("Fraction Gauss fell out of 0-1 range, aborting.")
                }
            };

            let grad_num: f64 = 1.0/denom;
            let grad_denom: f64 = -value/denom;

            for i in 0 .. grad.len() { grad[i] = 0.0; }
            for i in 0 .. n_lhs { 
                df( obj.lhs[i], grad );
                grad[ obj.lhs[i]*4 + 1] *= grad_num;
                grad[ obj.lhs[i]*4 + 2] *= grad_num;
                grad[ obj.lhs[i]*4 + 3] *= grad_num; 
            };
            for i in 0 .. n_rhs { 
                df( obj.rhs[i], grad );
                grad[ obj.rhs[i]*4 + 1] *= grad_denom; 
                grad[ obj.rhs[i]*4 + 2] *= grad_denom;
                grad[ obj.rhs[i]*4 + 3] *= grad_denom;
            };
        }

        value - obj.offset
    }
}

