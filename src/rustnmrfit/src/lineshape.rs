use ndarray::prelude::*;
use num::complex::Complex;
use std::f64::consts::SQRT_2;

use faddeeva;

// Some constants
const SQRT_PI: f64 = 1.772453850905515881919427556567825376987457275390625_f64;

//==============================================================================
// Lineshape struct where each x is unique

pub struct Lineshape1D {

    // Index of parameters that correspond to starting points of new peaks
    indexes: Array<usize, Ix1>,

    // Unique set of x values
    pub x: Array<f64, Ix1>,

    // Fit data, grouped by x value
    pub y: Array<f64, Ix2>,

    // Gradients, grouped by parameter and x value
    pub dydp: Array<f64, Ix3>,

}

//--------------------------------------
impl Lineshape1D {

    //--------------------------------------
    pub fn new(x: Array<f64, Ix1>, indexes: Array<usize, Ix1>) -> Lineshape1D {

        let n = x.len();
        let np = indexes.len()*4;

        Lineshape1D {
            indexes: indexes,
            x: x,
            y: Array::zeros((2, n)),
            dydp: Array::zeros((2, np, n)),
        }
    }

    //--------------------------------------
    pub fn calculate(&mut self, p: &[f64]) {

        // Initializing y to zero (unlike dydp, y is incremented)
        self.y.fill(0.0);

        // Loop through each peak
        for &i in self.indexes.iter() {

            let p_slice = &p[i .. (i+4)];
            let y_slice = self.y.slice_mut(s![.., ..]);
            let dydp_slice = self.dydp.slice_mut(s![.., i .. (i+4), ..]);

            calculate_peak(&self.x, &p_slice, y_slice, dydp_slice);
        }
    }

}

//==============================================================================
// Lineshape struct where some x are repeated

pub struct LineshapeND {

    // Index of parameters that correspond to starting points of new peaks
    np: usize,
    indexes: Array<usize, Ix1>,

    // Unique set of x values
    pub x: Array<f64, Ix1>,

    // Indexes mapping x values to y values
    xi: Array<usize, Ix1>,

    // Fit data, grouped by x value
    y_temp: Array<f64, Ix2>,
    pub y: Array<f64, Ix2>,

    // Gradients, grouped by parameter and x value
    dydp_temp: Array<f64, Ix3>,
    pub dydp: Array<f64, Ix3>,

}

//--------------------------------------
impl LineshapeND {

    //--------------------------------------
    pub fn new(x: Array<f64, Ix1>, xi: Array<usize, Ix1>, indexes: Array<usize, Ix1>) -> LineshapeND {

        let n = x.len();
        let ni = xi.len();
        let np = indexes.len()*4;

        LineshapeND {
            np: np,
            indexes: indexes,
            x: x,
            xi: xi,
            y_temp: Array::zeros((2, n)),
            y: Array::zeros((2, ni)),
            dydp_temp: Array::zeros((2, np, n)),
            dydp: Array::zeros((2, np, ni)),
        }
    }

    //--------------------------------------
    pub fn calculate(&mut self, p: &[f64]) {

        // Initializing y to zero (unlike dydp, y is incremented)
        self.y_temp.fill(0.0);

        // Loop through each peak using unique x values
        for &i in self.indexes.iter() {

            let p_slice = &p[i .. (i+4)];
            let y_temp_slice = self.y_temp.slice_mut(s![.., ..]);
            let dydp_temp_slice = self.dydp_temp.slice_mut(s![.., i .. (i+4), ..]);

            calculate_peak(&self.x, &p_slice, y_temp_slice, dydp_temp_slice);
        }

        // And then map unique values to repeats
        for i in 0 .. 2 {

            // First, the y values themselves
            let mut y_slice: ArrayViewMut<_, Ix1> = self.y.slice_mut(s![i, ..]);   
            let y_temp_slice: ArrayView<_, Ix1> = self.y_temp.slice(s![i, ..]);   

            for (j, &k) in self.xi.iter().enumerate() {
                y_slice[j] = y_temp_slice[k];
            }

            // Then, the derivatives
            for j in 0 .. self.np {
                let mut dydp_slice: ArrayViewMut<_, Ix1> = self.dydp.slice_mut(s![i, j, ..]);   
                let dydp_temp_slice: ArrayView<_, Ix1> = self.dydp_temp.slice(s![i, j, ..]);   

                for (j, &k) in self.xi.iter().enumerate() {
                    dydp_slice[j] = dydp_temp_slice[k];
                }
            }

        }
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

#[cfg(test)]
mod tests {

    use ndarray::prelude::*;

    #[test]
    fn nd_vs_1d() {

        let x: Array<f64, Ix1> = Array::linspace(0.0, 1.0, 5);

        let xi = vec![0, 0, 1, 1, 2, 2, 3, 3, 4, 4];
        let xi: Array<usize, Ix1> = Array::from_shape_vec((xi.len(),), xi).unwrap();
        let indexes =  vec![0, 4];
        let indexes: Array<usize, Ix1> = Array::from_shape_vec((indexes.len(),), indexes).unwrap();

        let mut l1 = super::Lineshape1D::new(x.clone(), indexes.clone());
        let mut ln = super::LineshapeND::new(x.clone(), xi, indexes.clone());

        let p = vec![0.5, 0.5, 0.5, 0.5, 0.3, 0.3, 0.3, 0.3];

        l1.calculate(&p[..]);
        ln.calculate(&p[..]);

        for i in 0 .. x.len() {

            for j in 0 .. 2 {

                assert!(l1.y[[0, i]] == ln.y[[0, i*2+j]], 
                        "ND vs 1D fit mismatch");
                assert!(l1.y[[1, i]] == ln.y[[1, i*2+j]], 
                        "ND vs 1D fit mismatch");

                for k in 0 .. indexes.len()*4 {

                    assert!(l1.dydp[[0, k, i]] == ln.dydp[[0, k, i*2+j]], 
                            "ND vs 1D gradient mismatch");
                    assert!(l1.dydp[[1, k, i]] == ln.dydp[[1, k, i*2+j]],
                            "ND vs 1D gradient mismatch");

                }
            }
        }
    }
}
