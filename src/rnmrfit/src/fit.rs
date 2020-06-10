use ndarray::prelude::*;

use crate::common::{NMRFitComponent, NMRFit};
use crate::lineshape::{Lineshape1D, Lineshape2D};
use crate::baseline::{Baseline1D, Baseline2D};
use crate::phase::{Phase1D, Phase2D};

//==============================================================================
// Fit1D 

//------------------------------------------------------------------------------
pub struct Fit1D {

    // Number of parameters (p) dedicated to each component of fit
    nl: usize,
    nb: usize,
    np: usize,

    // Fit components
    lineshape: Lineshape1D,
    baseline: Baseline1D,
    phase: Phase1D,

    // Intermediate matrices
    y_diff: Array2<f64>,

    // Primary output
    pub y_fit: Array2<f64>,
    pub grad: Array2<f64>,

}

//------------------------------------------------------------------------------
impl Fit1D {

    //--------------------------------------
    pub fn new(x: Array1<f64>, y: Array2<f64>, knots: Array1<f64>,
               nl: usize, nb: usize, np: usize) 
        -> Fit1D {

        let n = x.len();
        let x = x.into_shared();

        Fit1D {
            lineshape: Lineshape1D::new(x.clone(), nl),
            baseline: Baseline1D::new(x.clone(), knots, nb),
            phase: Phase1D::new(x.clone(), y, np),

            nl: nl,
            nb: nb,
            np: np,

            y_diff: Array::zeros((2, n)),
            y_fit: Array::zeros((2, n)),
            grad: Array::zeros((2, nl + nb*2 + np)),
        }
    }

    //--------------------------------------
    pub fn obj(p: &[f64], grad: Option<&mut [f64]>, obj: &mut Fit1D) -> f64 {
        obj.eval(p, grad)
    }

}

//--------------------------------------
impl NMRFit for Fit1D {

    //--------------------------------------
    fn eval(&mut self, p: &[f64], grad: Option<&mut [f64]>) -> f64 {

        let nl = self.nl;
        let nb = self.nb;
        let np = self.np;

        // Initializing fit
        self.y_fit.fill(0.0);

        // Phasing
        self.phase.eval(&p[(nl + nb*2) ..]);
            
        // Lineshape
        self.lineshape.eval(&p[.. nl]);

        self.y_fit += &self.lineshape.y;

        // Baseline
        self.baseline.eval(&p[nl .. (nl + nb*2)]);

        if nb > 0 {
            self.y_fit += &self.baseline.y;
        }

        //--------------------------------------
        // Gradients are last because they require the calculation of y-yfit

        // Calculating total difference
        self.y_diff = &self.phase.y - &self.y_fit;

        // Combine real/imaginary components into one loop
        for i in 0 .. 2 {

            let y_diff_slice: ArrayView1<_> = self.y_diff.slice(s![i, ..]);

            // Lineshape
            let dydp_slice: ArrayView2<_> = self.lineshape.dydp.slice(s![i, .., ..]);
            let mut grad_slice: ArrayViewMut1<_> = self.grad.slice_mut(s![i, 0 .. nl]);

            grad_slice.assign(&( -2.0*dydp_slice.dot(&y_diff_slice) ));

            // Baseline
            if nb > 0 {
                let nb1 = nb*i;
                let nb2 = nb*(i+1);

                let dydp_slice: ArrayView2<_> = self.baseline.dydp.slice(s![i, nb1 .. nb2, ..]);
                let mut grad_slice: ArrayViewMut1<_> = self.grad.slice_mut(s![i, (nl + nb1) .. (nl + nb2)]);

                grad_slice.assign(&( -2.0*dydp_slice.dot(&y_diff_slice) ));
            }

            // Phase
            if np > 0 {
                let dydp_slice: ArrayView2<_> = self.phase.dydp.slice(s![i, .., ..]);
                let mut grad_slice: ArrayViewMut1<_> = self.grad.slice_mut(s![i, (nl + nb*2) ..]);

                grad_slice.assign(&( 2.0*dydp_slice.dot(&y_diff_slice) ));
            }

        }

        // Filling in the overall gradient
        if grad.is_some() {
            let grad = grad.unwrap();

            for i in 0 .. grad.len() {
                grad[i] = self.grad.slice(s![.., i]).sum();
            }
        }

        
        
        (&self.y_diff * &self.y_diff).sum()
    }

}

//==============================================================================
// Fit2D 

//------------------------------------------------------------------------------
pub struct Fit2D {

    // Number of parameters (p) dedicated to each component of fit
    nl: usize,
    nb: usize,
    np: usize,

    // Fit components
    lineshape: Lineshape2D,
    baseline: Baseline2D,
    phase: Phase2D,

    // Intermediate matrices
    y_diff: Array2<f64>,

    // Primary output
    pub y_fit: Array2<f64>,
    pub grad: Array2<f64>,

}

//------------------------------------------------------------------------------
impl Fit2D {

    //--------------------------------------
    pub fn new(x_direct: Array1<f64>, x_indirect: Array1<f64>, y: Array2<f64>, 
               resonances: Array1<usize>, dimensions: Array1<usize>, knots: Array1<f64>,
               nl: usize, nb: usize, np: usize) 
        -> Fit2D {

        let n = x_direct.len();
        let x_direct = x_direct.into_shared();
        let x_indirect = x_indirect.into_shared();

        Fit2D {
            lineshape: Lineshape2D::new(x_direct.clone(), x_indirect.clone(), resonances, dimensions),
            baseline: Baseline2D::new(x_direct.clone(), x_indirect.clone(), knots, nb),
            phase: Phase2D::new(x_direct.clone(), x_indirect.clone(), y, np),

            nl: nl,
            nb: nb,
            np: np,

            y_diff: Array::zeros((4, n)),
            y_fit: Array::zeros((4, n)),
            grad: Array::zeros((4, nl + nb*4 + np)),
        }
    }

    //--------------------------------------
    pub fn obj(p: &[f64], grad: Option<&mut [f64]>, obj: &mut Fit2D) -> f64 {
        obj.eval(p, grad)
    }

}

//--------------------------------------
impl NMRFit for Fit2D {

    //--------------------------------------
    fn eval(&mut self, p: &[f64], grad: Option<&mut [f64]>) -> f64 {

        let nl = self.nl;
        let nb = self.nb;
        let np = self.np;

        // Initializing fit
        self.y_fit.fill(0.0);

        // Phasing
        self.phase.eval(&p[(nl + nb*4) ..]);
            
        // Lineshape
        self.lineshape.eval(&p[.. nl]);

        self.y_fit += &self.lineshape.y;

        // Baseline
        self.baseline.eval(&p[nl .. (nl + nb*4)]);

        if nb > 0 {
            self.y_fit += &self.baseline.y;
        }

        //--------------------------------------
        // Gradients are last because they require the calculation of y-yfit

        // Calculating total difference
        self.y_diff = &self.phase.y - &self.y_fit;

        // Combine real/imaginary components into one loop
        for i in 0 .. 4 {

            let y_diff_slice: ArrayView1<_> = self.y_diff.slice(s![i, ..]);

            // Lineshape
            let dydp_slice: ArrayView2<_> = self.lineshape.dydp.slice(s![i, .., ..]);
            let mut grad_slice: ArrayViewMut1<_> = self.grad.slice_mut(s![i, 0 .. nl]);

            grad_slice.assign(&( -2.0*dydp_slice.dot(&y_diff_slice) ));

            // Baseline
            if nb > 0 {
                let nb1 = nb*i;
                let nb2 = nb*(i+1);

                let dydp_slice: ArrayView2<_> = self.baseline.dydp.slice(s![i, nb1 .. nb2, ..]);
                let mut grad_slice: ArrayViewMut1<_> = self.grad.slice_mut(s![i, (nl + nb1) .. (nl + nb2)]);

                grad_slice.assign(&( -2.0*dydp_slice.dot(&y_diff_slice) ));
            }

            // Phase
            if np > 0 {
                let dydp_slice: ArrayView2<_> = self.phase.dydp.slice(s![i, .., ..]);
                let mut grad_slice: ArrayViewMut1<_> = self.grad.slice_mut(s![i, (nl + nb*4) ..]);

                grad_slice.assign(&( 2.0*dydp_slice.dot(&y_diff_slice) ));
            }

        }

        // Filling in the overall gradient
        if grad.is_some() {
            let grad = grad.unwrap();

            for i in 0 .. grad.len() {
                grad[i] = self.grad.slice(s![.., i]).sum();
            }
        }

        
        
        (&self.y_diff * &self.y_diff).sum()
    }

}

//==============================================================================
// Unit tests

#[cfg(test)]
mod tests {

    use ndarray::prelude::*;
    use ndarray_rand::RandomExt;
    use ndarray_rand::rand::SeedableRng;
    use ndarray_rand::rand_distr::{Normal, Distribution};
    use rand_isaac::isaac64::Isaac64Rng;

    use crate::common::{NMRFitComponent, NMRFit};
    use crate::lineshape::{Lineshape1D, Lineshape2D};

    //--------------------------------------
    #[test]
    fn fit_1d_gradient() {

        let x: Array1<f64> = Array::linspace(0.0, 1.0, 20);
        let x = x.into_shared();
        let p = vec![0.5, 0.5, 0.5, 0.5, 0.3, 0.3, 0.3, 0.3, 
                     0.1, 0.2, 0.3, 0.4, 0.5, 0.1, 0.2, 0.3, 0.4, 0.5, 
                     0.1, 0.05];

        let nl = 8;
        let nb = 5;
        let np = 2;

        // Double checking basic math
        assert!((nl + nb*2 + np) == p.len(), "Error setting up test.");

        // Generate realistic peaks
        let mut lineshape = Lineshape1D::new(x.clone(), nl);
        lineshape.eval(&p[0 .. nl]);
        let mut y = lineshape.y.clone();

        // Adding some noise for a bit of realism
        let seed = 1111;
        let mut rng = Isaac64Rng::seed_from_u64(seed);
        let noise = Array::random_using(y.raw_dim(), Normal::new(0., 0.01).unwrap(), &mut rng);
        y += &noise;

        // Basic boundary knots
        let knots: Array1<f64> = Array::from_shape_vec((2,), vec![0.0, 1.0]).unwrap();

        // Initialize fit
        let mut fit = super::Fit1D::new(x.to_owned(), y, knots, nl, nb, np);

        // Summing up errors across multiple points/dimensions increases threshold
        fit.grad_test(&p, 5e-4);

    }

    //--------------------------------------
    #[test]
    fn fit_2d_gradient() {

        let x: Array1<f64> = Array::linspace(0.0, 1.0, 20);
        let n = x.len();
        let p = vec![0.5, 0.5, 0.5, 0.5, 0.3, 0.3, 0.3, 0.3, 
                     0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4, 0.1, 0.2, 0.3, 0.4,
                     0.1, 0.05, 0.05];

        let nl = 8;
        let nb = 4;
        let np = 3;

        // Double checking basic math
        assert!((nl + nb*4 + np) == p.len(), "Error setting up test.");

        // Building grid from x
        let mut x1: Array1<f64> = Array::zeros((n*n,));
        let mut x2 = x1.clone();

        for i in 0 .. n {
            for j in 0 .. n {
                x1[i+n*j] = x[j];
                x2[i+n*j] = x[i];
            }
        }

        let x1 = x1.into_shared();
        let x2 = x2.into_shared();

        // Generate realistic peaks
        let resonances: Array1<usize> = Array::zeros((nl,));
        let mut dimensions: Array1<usize> = Array::zeros((nl,));
        for i in 4 .. 8 {
            dimensions[i] = 1;
        }

        let mut lineshape = Lineshape2D::new(x1.clone(), x2.clone(), 
                                             resonances.clone(), dimensions.clone());
        lineshape.eval(&p[0 .. nl]);
        let mut y = lineshape.y.clone();

        // Adding some noise for a bit of realism
        let seed = 1111;
        let mut rng = Isaac64Rng::seed_from_u64(seed);
        let noise = Array::random_using(y.raw_dim(), Normal::new(0., 0.01).unwrap(), &mut rng);
        y += &noise;

        // Basic boundary knots
        let knots: Array1<f64> = Array::from_shape_vec((2,), vec![0.0, 1.0]).unwrap();

        // Initialize fit
        let mut fit = super::Fit2D::new(x1.to_owned(), x2.to_owned(), y,
                                        resonances, dimensions, knots, 
                                        nl, nb, np);

        // Summing up errors across multiple points/dimensions increases threshold
        fit.grad_test(&p, 5e-4);

    }

}
