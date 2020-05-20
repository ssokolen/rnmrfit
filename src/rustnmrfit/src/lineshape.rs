use ndarray::{prelude::*, ArcArray1, Zip};
use num::complex::Complex;

use crate::peak::{Peak, PeakFunctions};

//==============================================================================
// General 1D lineshape

pub struct Lineshape1D {

    // Reference to x values for convenience
    x: ArcArray1<f64>,

    // Final output
    pub y: Array2<f64>,
    pub dydp: Array3<f64>,

}

//--------------------------------------
impl Lineshape1D {

    //--------------------------------------
    pub fn new(x: ArcArray1<f64>, np: usize) -> Lineshape1D {

        let n = x.len();

        Lineshape1D {
            x: x,
            
            y: Array::zeros((2, n)),
            dydp: Array::zeros((2, np, n)),
        }
    }

    //--------------------------------------
    pub fn eval(&mut self, p: &Array1<f64>) {

        // Initializing y to zero (unlike dydp, y is incremented)
        self.y.fill(0.0);

         // Loop through each peak
        for i in (0 .. p.len()).step_by(4) {

            // Generate peak object and match based on result
            let peak = Peak::new(p[i], p[i+1], p[i+2], p[i+3]);

            match peak {
                Peak::Lorentz( lorentz ) => self.eval_peak(lorentz, i),
                Peak::Voigt( voigt ) => self.eval_peak(voigt, i),
            };
        }
    }   

    //------------------------------------------------------------------------------
    fn eval_peak(&mut self, mut peak: impl PeakFunctions, i: usize) {

        // Temporary holder for gradient terms
        let mut grad: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); 4];

        let mut y = self.y.slice_mut(s![.., ..]);
        let mut dydp = self.dydp.slice_mut(s![.., i .. (i+4), ..]);

        // Looping over x
        for j in 0 .. self.x.len() {

            let fit = peak.gradients(self.x[j], &mut grad);

            y[[0, j]] += fit.re;
            y[[1, j]] += fit.im;

            for k in 0 .. 4 {

                dydp[[0, k, j]] = grad[k].re;
                dydp[[1, k, j]] = grad[k].im;

            }
        }
    }
}

//==============================================================================
// Container for 1D lineshape that maps repeated values

pub struct Lineshape1DMap {

    // Reference to x indices for convenience
    x_map: ArcArray1<usize>,

    // Lineshape object for unique values
    lineshape: Lineshape1D,

    // Final output
    pub y: Array2<f64>,
    pub dydp: Array3<f64>,

}

//--------------------------------------
impl Lineshape1DMap {

    //--------------------------------------
    pub fn new(x: ArcArray1<f64>, x_map: ArcArray1<usize>, np: usize) -> Lineshape1DMap {

        let n = x_map.len();

        Lineshape1DMap {
            x_map: x_map,

            lineshape: Lineshape1D::new(x, np),
            
            y: Array::zeros((2, n)),
            dydp: Array::zeros((2, np, n)),
        }
    }

    //--------------------------------------
    pub fn eval(&mut self, p: &Array1<f64>) {       

        // First, evaluate unique values
        self.lineshape.eval(p);

        // And then map unique values to repeats
        let np = p.len();

        for i in 0 .. 2 {

            // First, the y values themselves
            let from: ArrayView<_, Ix1> = self.lineshape.y.slice(s![i, ..]);   
            let to: ArrayViewMut<_, Ix1> = self.y.slice_mut(s![i, ..]);   

            Zip::from(to).and(&self.x_map).apply(|x, &i| *x = from[i]);

            // Then, the derivatives
            for j in 0 .. np {
                let from: ArrayView<_, Ix1> = self.lineshape.dydp.slice(s![i, j, ..]);   
                let to: ArrayViewMut<_, Ix1> = self.dydp.slice_mut(s![i, j, ..]);   

                Zip::from(to).and(&self.x_map).apply(|x, &i| *x = from[i]);

            }

        }
    }
}

//==============================================================================
// Double checking that map works correctly 

#[cfg(test)]
mod tests {

    use ndarray::prelude::*;

    #[test]
    fn map() {

        let x  = Array::linspace(0.0, 1.0, 5).into_shared();

        let xi = vec![0, 0, 1, 1, 2, 2, 3, 3, 4, 4];
        let xi = Array::from_shape_vec((xi.len(),), xi).unwrap().into_shared();

        let p = vec![0.5, 0.5, 0.5, 0.5, 0.3, 0.3, 0.3, 0.3];
        let p = Array::from_shape_vec((p.len(),), p).unwrap();

        let mut l = super::Lineshape1D::new(x.clone(), p.len());
        let mut lmap = super::Lineshape1DMap::new(x.clone(), xi.clone(), p.len());

        l.eval(&p);
        lmap.eval(&p);

        for i in 0 .. x.len() {

            for j in 0 .. 2 {

                assert!(l.y[[0, i]] == lmap.y[[0, i*2+j]], 
                        "ND vs 1D fit mismatch");
                assert!(l.y[[1, i]] == lmap.y[[1, i*2+j]], 
                        "ND vs 1D fit mismatch");

                for k in 0 .. p.len() {

                    assert!(l.dydp[[0, k, i]] == lmap.dydp[[0, k, i*2+j]], 
                            "ND vs 1D gradient mismatch");
                    assert!(l.dydp[[1, k, i]] == lmap.dydp[[1, k, i*2+j]],
                            "ND vs 1D gradient mismatch");

                }
            }
        }
    }
}
