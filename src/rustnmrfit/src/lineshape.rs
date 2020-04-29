use std::f64::consts::SQRT_2;
use num::complex::Complex;
use faddeeva;

use ndarray::prelude::*;

// Some constants
const SQRT_PI: f64 = 1.772453850905515881919427556567825376987457275390625_f64;

//==============================================================================
// Main lineshape struct

pub struct Lineshape1D {

    // Number of parameters
    np: usize,

    // Unique set of x values
    pub x: Array<f64, Ix1>,

    // Fit data, grouped by resonance and x value
    pub y: Array<f64, Ix3>,

    // Gradients, grouped by parameter and x value
    pub dydp: Array<f64, Ix3>,

}

//--------------------------------------
impl Lineshape1D {

    //--------------------------------------
    pub fn new(x: Array<f64, Ix1>, np: usize, nr: usize) -> Lineshape1D {

        let n = x.len();

        Lineshape1D {
            np: np,
            x: x,
            y: Array::zeros((2, nr, n)),
            dydp: Array::zeros((2, np, n)),
        }
    }

    //--------------------------------------
    pub fn new_flat(x: Array<f64, Ix1>, np: usize) -> Lineshape1D {
        
        Lineshape1D::new(x, np, 1)
    }

    //--------------------------------------
    pub fn calculate(&mut self, p: &[f64], r: &[usize]) {

        // Initializing y to zero (unlike dydp, y is incremented)
        self.y.fill(0.0);

        // Loop through each peak
        for i in (0 .. self.np).step_by(4) {

            let p_slice = &p[i .. (i+4)];
            let y_slice = self.y.slice_mut(s![.., r[i], ..]);
            let dydp_slice = self.dydp.slice_mut(s![.., i .. (i+4), ..]);

            calculate_peak(&self.x, &p_slice, y_slice, dydp_slice);
        }
    }

    //--------------------------------------
    pub fn calculate_flat(&mut self, p: &[f64]) {

        let r: Vec<usize> = vec![0; self.np]; 

        self.calculate(p, &r[..]);

    }

}

//==============================================================================
// Peak calculating functions

//------------------------------------------------------------------------------
pub fn calculate_peak(x: &Array<f64, Ix1>, p: &[f64], 
                      y_fit: ArrayViewMut<f64, Ix2>, 
                      y_grad: ArrayViewMut<f64, Ix3>) {

    if p[3] <= 1e-8 {
        calculate_lorentz(x, p, y_fit, y_grad) 
    } else if p[3] <= (1.0 - 1e-8) {
        calculate_voigt(x, p, y_fit, y_grad) 
    } else {
        panic!("Gauss peaks not supported yet.")
    }


}

//------------------------------------------------------------------------------
fn calculate_lorentz(x: &Array<f64, Ix1>, param: &[f64],
                         mut y_fit: ArrayViewMut<f64, Ix2>, 
                         mut y_grad: ArrayViewMut<f64, Ix3>) {

    // Unpacking parameters to avoid typos
    let p: f64 = param[0];
    let w: f64 = param[1];
    let h: f64 = param[2];

    let n: usize = x.len();

    // Terms independent of x
    let dzdp = Complex::new(1.0 / w, 0.0);

    // Looping over x
    for i in 0 .. n {

        // Some common terms
        let z = (p - x[i]) / w;
        let z_2 = z * z;

        // Calculating an intermediate y with no h term
        let yo = Complex::new(1.0, z)  / 
                 Complex::new(z_2 + 1.0, 0.0);
           
        // Tacking on the h for the final y value
        let y = Complex::new(h, 0.0) * yo;

        // Storing output
        y_fit[[0, i]] += y.re;
        y_fit[[1, i]] += y.im;

        // Derivative of y (for chain rule)
        let dydz = Complex::new(-2.0 * z, 1.0 - z_2) * 
                   Complex::new(h / ( (z_2 + 1.0) * (z_2 + 1.0) ), 0.0);

        // Derivatives of z (for chain rule)
        let dzdw = Complex::new(-z / w, 0.0);

        // Position gradient
        let dydp = dydz*dzdp;
        y_grad[[0, 0, i]] = dydp.re;
        y_grad[[1, 0, i]] = dydp.im;

        // Lorentz width gradient (corrected dzndw term is 0, cancelling out other terms)
        let dydw = dydz*dzdw;
        y_grad[[0, 1, i]] = dydw.re;
        y_grad[[1, 1, i]] = dydw.im;

        // Height gradient (consistent at yo)
        y_grad[[0, 2, i]] = yo.re;
        y_grad[[1, 2, i]] = yo.im;

        // Fraction gradient ignored (consistent at 0)
    }
}

//------------------------------------------------------------------------------
fn calculate_voigt(x: &Array<f64, Ix1>, param: &[f64], 
                       mut y_fit: ArrayViewMut<f64, Ix2>, 
                       mut y_grad: ArrayViewMut<f64, Ix3>) {

    // Unpacking parameters to avoid typos
    let p: f64 = param[0];
    let w: f64 = param[1];
    let h: f64 = param[2];
    let f: f64 = param[3];

    let wg = w * f / (1.0 - f);

    let n: usize = x.len();

    // Terms independent of x
    let zn = Complex::new( 0.0, w / (SQRT_2 * wg) );
    let yn = faddeeva::w( zn, 1e-6 );

    let a = Complex::new( h, 0.0 ) / yn;

    let dyndzn = Complex::new( -2.0, 0.0 ) * zn * yn + 
                 Complex::new( 0.0, 2.0/SQRT_PI );  

    let b = Complex::new( wg / (w * f * f), 0.0 );
    let dzndf = -zn * b;

    let dzdp = Complex::new( 1.0/( SQRT_2 * wg ), 0.0 );

    // Looping over x
    for i in 0 .. n {

        // Some common terms
        let z = Complex::new( p - x[i] , w )  / 
                Complex::new( SQRT_2 * wg, 0.0 );

        // Calculating an intermediate f with no h term
        let yo = faddeeva::w( z, 1e-6 );

        // Tacking on the h and yn for the final y value
        let y = a * yo;

        // Storing output
        y_fit[[0, i]] += y.re;
        y_fit[[1, i]] += y.im;

        // Derivative of f (for chain rule)
        let dyodz = Complex::new( -2.0, 0.0 ) * z * yo + 
                    Complex::new( 0.0, 2.0/SQRT_PI );  

        // Derivatives of z (for chain rule)
        let dzdw = Complex::new( -(p - x[i])/(SQRT_2 * w * wg), 0.0 );
        let dzdf = -z * b;

        // Position gradient
        let dydp = a * dyodz * dzdp;
        y_grad[[0, 0, i]] = dydp.re;
        y_grad[[1, 0, i]] = dydp.im;

        // Lorentz width gradient (corrected dzndw term is 0, cancelling out other terms)
        let dydw = a * dyodz * dzdw;
        y_grad[[0, 1, i]] = dydw.re;
        y_grad[[1, 1, i]] = dydw.im;

        // Height gradient
        let dydh = yo / yn;
        y_grad[[0, 2, i]] = dydh.re;
        y_grad[[1, 2, i]] = dydh.im;

        // Fraction gradient
        let dydf = a * (dyodz * dzdf * yn - dyndzn * dzndf * yo) / yn;
        y_grad[[0, 3, i]] = dydf.re;
        y_grad[[1, 3, i]] = dydf.im;
    }
}
