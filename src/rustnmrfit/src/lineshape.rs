use ndarray::{prelude::*, Data, DataMut, Zip};
use num::complex::Complex;

use crate::peak::{Peak, PeakFunctions};

//==============================================================================
// General 1D lineshape

pub struct Lineshape1D {

    // Intermediate values
    y_temp: Array2<f64>,
    dydp_temp: Array3<f64>,

    // Final output
    pub y: Array2<f64>,
    pub dydp: Array3<f64>,

}

//--------------------------------------
impl Lineshape1D {

    //--------------------------------------
    pub fn new(n: usize, nu: usize, np: usize) -> Lineshape1D {

        Lineshape1D {
            y_temp: Array::zeros((2, nu)),
            dydp_temp: Array::zeros((2, np, nu)),
            
            y: Array::zeros((2, n)),
            dydp: Array::zeros((2, np, n)),
        }
    }

    //------------------------------------------------------------------------------
    pub fn eval_peak<S,T>(mut peak: impl PeakFunctions,  x: &ArrayBase<S, Ix1>, 
                          mut y: ArrayBase<T, Ix2>, mut dydp: ArrayBase<T, Ix3>)
    where S: Data<Elem = f64>, T: DataMut<Elem = f64> {

        // Temporary holder for gradient terms
        let mut grad: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); 4];

        // Looping over x
        for i in 0 .. x.len() {

            let fit = peak.gradients(x[i], &mut grad);

            y[[0, i]] += fit.re;
            y[[1, i]] += fit.im;

            for j in 0 .. 4 {

                dydp[[0, j, i]] = grad[j].re;
                dydp[[1, j, i]] = grad[j].im;

            }
        }
    }

    //--------------------------------------
    fn eval_peaks<S,T>(x: &ArrayBase<S, Ix1>, p: &ArrayBase<S, Ix1>,
                        y: &mut ArrayBase<T, Ix2>,  dydp: &mut ArrayBase<T, Ix3>)
    where S: Data<Elem = f64>, T: DataMut<Elem = f64> {

        // Initializing y to zero (unlike dydp, y is incremented)
        y.fill(0.0);

         // Loop through each peak
        for i in (0 .. p.len()).step_by(4) {

            // Generate peak object and match based on result
            let peak = Peak::new(p[i], p[i+1], p[i+2], p[i+3]);

            let y_slice = y.slice_mut(s![.., ..]);
            let dydp_slice = dydp.slice_mut(s![.., i .. (i+4), ..]);

            match peak {
                Peak::Lorentz( lorentz ) => Lineshape1D::eval_peak(lorentz, x, y_slice, dydp_slice),
                Peak::Voigt( voigt ) => Lineshape1D::eval_peak(voigt, x, y_slice, dydp_slice),
            };
        }
    }       

    //--------------------------------------
    pub fn eval<S>(&mut self, x: &ArrayBase<S, Ix1>, p: &ArrayBase<S, Ix1>)
    where S: Data<Elem = f64> {       

        // If x values are not being mapped to duplicates, feed y directly into eval_peaks
        Lineshape1D::eval_peaks(x, p, &mut self.y, &mut self.dydp);

    }

    //--------------------------------------
    pub fn eval_map<S,T>(&mut self, x: &ArrayBase<S, Ix1>, x_map: &ArrayBase<T, Ix1>, p: &ArrayBase<S, Ix1>)
    where S: Data<Elem = f64>, T: Data<Elem = usize> {       

        // If x values are being mapped, then feed temporary arrays first
        Lineshape1D::eval_peaks(x, p, &mut self.y_temp, &mut self.dydp_temp);

        let np = p.len();

        // And then map unique values to repeats
        for i in 0 .. 2 {

            // First, the y values themselves
            let from: ArrayView<_, Ix1> = self.y_temp.slice(s![i, ..]);   
            let to: ArrayViewMut<_, Ix1> = self.y.slice_mut(s![i, ..]);   

            Zip::from(to).and(x_map).apply(|x, &i| *x = from[i]);

            // Then, the derivatives
            for j in 0 .. np {
                let from: ArrayView<_, Ix1> = self.dydp_temp.slice(s![i, j, ..]);   
                let to: ArrayViewMut<_, Ix1> = self.dydp.slice_mut(s![i, j, ..]);   

                Zip::from(to).and(x_map).apply(|x, &i| *x = from[i]);

            }

        }
    }
}

//==============================================================================
// Peak calculating functions

#[cfg(test)]
mod tests {

    use ndarray::prelude::*;

    #[test]
    fn map() {

        let x: Array<f64, Ix1> = Array::linspace(0.0, 1.0, 5);

        let xi = vec![0, 0, 1, 1, 2, 2, 3, 3, 4, 4];
        let xi: Array<usize, Ix1> = Array::from_shape_vec((xi.len(),), xi).unwrap();

        let p = vec![0.5, 0.5, 0.5, 0.5, 0.3, 0.3, 0.3, 0.3];
        let p: Array<f64, Ix1> = Array::from_shape_vec((p.len(),), p).unwrap();

        let mut l1 = super::Lineshape1D::new(x.len(), x.len(), p.len());
        let mut ln = super::Lineshape1D::new(xi.len(), x.len(), p.len());

        l1.eval(&x, &p);
        ln.eval_map(&x, &xi, &p);

        for i in 0 .. x.len() {

            for j in 0 .. 2 {

                assert!(l1.y[[0, i]] == ln.y[[0, i*2+j]], 
                        "ND vs 1D fit mismatch");
                assert!(l1.y[[1, i]] == ln.y[[1, i*2+j]], 
                        "ND vs 1D fit mismatch");

                for k in 0 .. p.len() {

                    assert!(l1.dydp[[0, k, i]] == ln.dydp[[0, k, i*2+j]], 
                            "ND vs 1D gradient mismatch");
                    assert!(l1.dydp[[1, k, i]] == ln.dydp[[1, k, i*2+j]],
                            "ND vs 1D gradient mismatch");

                }
            }
        }
    }
}
