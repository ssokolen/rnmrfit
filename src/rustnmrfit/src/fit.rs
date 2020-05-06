use ndarray::prelude::*;
use std::iter::FromIterator;

use super::lineshape::Lineshape1D;
use super::lineshape::LineshapeND;

//==============================================================================
// Fit1D 

//------------------------------------------------------------------------------
pub struct Fit1D {

    // y is the observed spectral data
    y: Array<f64, Ix2>,

    // Number of parameters (p) dedicated to each component of fit
    nl: usize,
    nb: usize,
    np: usize,

    // NMRLineshape object holding gradients and lineshape fit
    lineshape: Lineshape1D,

    // Baseline b-spline basis
    basis: Array<f64, Ix2>,

    // Phase angles and intermediate terms
    theta: Array<f64, Ix1>,
    sin_theta: Array<f64, Ix1>,
    cos_theta: Array<f64, Ix1>,

    // Intermediate matrices
    y_mod: Array<f64, Ix2>,
    y_diff: Array<f64, Ix2>,
    grad_theta: Array<f64, Ix2>,

    // Primary output
    pub y_fit: Array<f64, Ix2>,
    pub grad: Array<f64, Ix2>,

}

//------------------------------------------------------------------------------
impl Fit1D {

    //--------------------------------------
    pub fn new(x: Array<f64, Ix1>, y: Array<f64, Ix2>, nl: usize, nb: usize, np: usize,
               basis: Option<Array<f64, Ix2>>) 
        -> Fit1D {

        let n = x.len();

        let basis: Array<f64, Ix2> = match basis {
            Some(array) => array,
            None => Array::zeros((n, 1)),
        };

        Fit1D {
            y: y,

            nl: nl,
            nb: nb,
            np: np,

            lineshape: Lineshape1D::new(x, Array::from_iter((0 .. nl).step_by(4))),
            basis: basis,

            theta: Array::zeros((n, )),
            sin_theta: Array::zeros((n, )),
            cos_theta: Array::zeros((n, )),

            y_mod: Array::zeros((2, n)),
            y_diff: Array::zeros((2, n)),
            grad_theta: Array::zeros((2, n)),

            y_fit: Array::zeros((2, n)),
            grad: Array::zeros((2, nl + 2*nb + np)),
        }
    }

    //--------------------------------------
    pub fn obj(p: &[f64], grad: Option<&mut [f64]>, obj: &mut Fit1D) -> f64 {
        obj.eval(p, grad)
    }

}

//------------------------------------------------------------------------------
impl Fit1D {

    //--------------------------------------
    pub fn eval(&mut self, p: &[f64], grad: Option<&mut [f64]>) -> f64 {

        let nl = self.nl;
        let nb = self.nb;
        let np = self.np;

        // Initializing fit
        self.y_fit.fill(0.0);

        let y_r: ArrayView<_, Ix1> = self.y.slice(s![0, ..]);
        let y_i: ArrayView<_, Ix1> = self.y.slice(s![1, ..]);

        //--------------------------------------
        // Phasing

        // Modifying y based on phase correction
        if np > 0 {

            if np == 1 {
                let theta_0 = p[nl + nb*2];
                self.theta.fill(theta_0);
                self.sin_theta.fill(theta_0.sin());
                self.cos_theta.fill(theta_0.cos());
            } else {
                let theta_0 = p[nl + nb*2];
                let theta_1 = p[nl + nb*2 + 1];
                self.theta.assign(&( self.lineshape.x.mapv(|x| theta_0 + theta_1*x) ));
                self.sin_theta.assign(&( self.theta.mapv(|x| x.sin()) ));
                self.cos_theta.assign(&( self.theta.mapv(|x| x.cos()) ));
            }

            let mut y_mod_r: ArrayViewMut<_, Ix1> = self.y_mod.slice_mut(s![0, ..]);
            y_mod_r.assign(&( &y_r * &self.cos_theta + &y_i * &self.sin_theta ));
            
            let mut y_mod_i: ArrayViewMut<_, Ix1> = self.y_mod.slice_mut(s![1, ..]);
            y_mod_i.assign(&( -&y_r * &self.sin_theta + &y_i * &self.cos_theta ));

        } else {    
            let mut y_mod_r: ArrayViewMut<_, Ix1> = self.y_mod.slice_mut(s![0, ..]);
            y_mod_r.assign(&y_r);
            
            let mut y_mod_i: ArrayViewMut<_, Ix1> = self.y_mod.slice_mut(s![1, ..]);
            y_mod_i.assign(&y_i);
        }
            
        //--------------------------------------
        // Lineshape

        // Calculating peak fit and gradient terms 
        self.lineshape.calculate(p);

        // Updating with lineshapes
        self.y_fit += &self.lineshape.y.slice(s![.., ..]);

        //--------------------------------------
        // Baseline

        // Tacking on baseline
        if nb > 0 {
            let pb_r = p[nl .. (nl + nb)].to_vec();
            let pb_r: Array<_, Ix1> = Array::from_shape_vec((nb, ), pb_r).unwrap(); 
            let mut y_fit_r: ArrayViewMut<_, Ix1> = self.y_fit.slice_mut(s![0, ..]);

            y_fit_r += &( self.basis.dot(&pb_r) );

            let pb_i = p[(nl + nb) .. (nl + nb*2)].to_vec();
            let pb_i: Array<_, Ix1> = Array::from_shape_vec((nb, ), pb_i).unwrap(); 
            let mut y_fit_i: ArrayViewMut<_, Ix1> = self.y_fit.slice_mut(s![1, ..]);

            y_fit_i += &( self.basis.dot(&pb_i) );
        }

        //--------------------------------------
        // All gradients

        // Calculating total difference
        self.y_diff = &self.y_mod - &self.y_fit;

        let y_diff_r: ArrayView<_, Ix1> = self.y_diff.slice(s![0, ..]);
        let y_diff_i: ArrayView<_, Ix1> = self.y_diff.slice(s![1, ..]);

        // Lineshape
        for i in 0 .. nl {
            let dydp_r: ArrayView<_, Ix1> = self.lineshape.dydp.slice(s![0, i, ..]);
            let dydp_i: ArrayView<_, Ix1> = self.lineshape.dydp.slice(s![1, i, ..]);

            self.grad[[0, i]] = y_diff_r.dot(&dydp_r);
            self.grad[[1, i]] = y_diff_i.dot(&dydp_i);
        } 

        // Baseline 
        if nb > 0 {
            let y_diff_r: ArrayView<_, Ix2> = self.y_diff.slice(s![0 .. 1, ..]);
            let y_diff_i: ArrayView<_, Ix2> = self.y_diff.slice(s![1 .. 2, ..]);

            let grad_r = y_diff_r.dot(&self.basis);
            let grad_i = y_diff_i.dot(&self.basis);

            for i in 0 .. nb {
                self.grad[[0, nl + i]] = grad_r[[0, i]];
                self.grad[[1, nl + nb + i]] = grad_i[[0, i]];
            } 
        }

        // Phase
        if np > 0 {

            let mut grad_theta_r: ArrayViewMut<_, Ix1> = self.grad_theta.slice_mut(s![0, ..]);
            grad_theta_r.assign(&( -&y_r * &self.sin_theta + &y_i * &self.cos_theta ));
            
            let mut grad_theta_i: ArrayViewMut<_, Ix1> = self.grad_theta.slice_mut(s![1, ..]);
            grad_theta_i.assign(&( -&y_r * &self.cos_theta - &y_i * &self.sin_theta ));

            let grad_theta_r: ArrayView<_, Ix1> = self.grad_theta.slice(s![0, ..]);
            let grad_theta_i: ArrayView<_, Ix1> = self.grad_theta.slice(s![1, ..]);

            // The negative is a hack to counteract the general negative below
            if np == 1 {
                self.grad[[0, nl + nb*2]] = -y_diff_r.dot(&grad_theta_r);
                self.grad[[1, nl + nb*2]] = -y_diff_i.dot(&grad_theta_i);
            } else {
                self.grad[[0, nl + nb*2]] = -y_diff_r.dot(&grad_theta_r);
                self.grad[[1, nl + nb*2]] = -y_diff_i.dot(&grad_theta_i);

                self.grad[[0, nl + nb*2 + 1]] = -y_diff_r.dot(
                    &( &grad_theta_r * &self.lineshape.x )
                );
                self.grad[[1, nl + nb*2 + 1]] = -y_diff_i.dot(
                    &( &grad_theta_i * &self.lineshape.x )
                );
            }

        }

        // Tacking on the -2 scalar term to the gradient
        self.grad *= -2.0;

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

    // y is the observed spectral data
    y: Array<f64, Ix2>,

    // Number of parameters (p) dedicated to each component of fit
    nl: usize,
    nb: usize,
    np: usize,

    // Lineshapes are held as pairs of direct and indirect components,
    // within vectors of resonances
    direct: Vec<LineshapeND>,
    indirect: Vec<LineshapeND>,

    // Baseline b-spline basis
    //basis: Array<f64, Ix2>,

    // Phase angles and intermediate terms
    //theta: Array<f64, Ix1>,
    //sin_theta: Array<f64, Ix1>,
    //cos_theta: Array<f64, Ix1>,

    // Intermediate matrices
    y_mod: Array<f64, Ix2>,
    y_diff: Array<f64, Ix2>,
    grad_theta: Array<f64, Ix2>,

    // Primary output
    pub y_fit: Array<f64, Ix2>,
    pub grad: Array<f64, Ix2>,

}
/*
//------------------------------------------------------------------------------
impl Fit2D {

    //--------------------------------------
    pub fn new(x_direct: Array<f64, Ix1>, x_indirect: Array<f64, Ix1>, 
               xi_direct: Array<usize, Ix1>, xi_indirect: Array<usize, Ix1>,
               
               y: Array<f64, Ix2>,  
               nl: usize, nb: usize, np: usize,
               _: Option<Array<f64, Ix2>>) 
        -> Fit2D {

        let n = x.len();

        let basis: Array<f64, Ix2> = match basis {
            Some(array) => array,
            None => Array::zeros((n, 1)),
        };

        Fit1D {
            y: y,

            nl: nl,
            nb: nb,
            np: np,

            lineshape: Lineshape1D::new(x, Array::from_iter((0 .. nl).step_by(4))),
            basis: basis,

            theta: Array::zeros((n, )),
            sin_theta: Array::zeros((n, )),
            cos_theta: Array::zeros((n, )),

            y_mod: Array::zeros((2, n)),
            y_diff: Array::zeros((2, n)),
            grad_theta: Array::zeros((2, n)),

            y_fit: Array::zeros((2, n)),
            grad: Array::zeros((2, nl + 2*nb + np)),
        }
    }

    //--------------------------------------
    pub fn obj(p: &[f64], grad: Option<&mut [f64]>, obj: &mut Fit1D) -> f64 {
        obj.eval(p, grad)
    }

}

//------------------------------------------------------------------------------
impl Fit2D {

    //--------------------------------------
    pub fn eval(&mut self, p: &[f64], grad: Option<&mut [f64]>) -> f64 {

        let nl = self.nl;
        let nb = self.nb;
        let np = self.np;

        // Initializing fit
        self.y_fit.fill(0.0);

        let y_r: ArrayView<_, Ix1> = self.y.slice(s![0, ..]);
        let y_i: ArrayView<_, Ix1> = self.y.slice(s![1, ..]);

        //--------------------------------------
        // Phasing

        // Modifying y based on phase correction
        if np > 0 {

            if np == 1 {
                let theta_0 = p[nl + nb*2];
                self.theta.fill(theta_0);
                self.sin_theta.fill(theta_0.sin());
                self.cos_theta.fill(theta_0.cos());
            } else {
                let theta_0 = p[nl + nb*2];
                let theta_1 = p[nl + nb*2 + 1];
                self.theta.assign(&( self.lineshape.x.mapv(|x| theta_0 + theta_1*x) ));
                self.sin_theta.assign(&( self.theta.mapv(|x| x.sin()) ));
                self.cos_theta.assign(&( self.theta.mapv(|x| x.cos()) ));
            }

            let mut y_mod_r: ArrayViewMut<_, Ix1> = self.y_mod.slice_mut(s![0, ..]);
            y_mod_r.assign(&( &y_r * &self.cos_theta + &y_i * &self.sin_theta ));
            
            let mut y_mod_i: ArrayViewMut<_, Ix1> = self.y_mod.slice_mut(s![1, ..]);
            y_mod_i.assign(&( -&y_r * &self.sin_theta + &y_i * &self.cos_theta ));

        } else {    
            let mut y_mod_r: ArrayViewMut<_, Ix1> = self.y_mod.slice_mut(s![0, ..]);
            y_mod_r.assign(&y_r);
            
            let mut y_mod_i: ArrayViewMut<_, Ix1> = self.y_mod.slice_mut(s![1, ..]);
            y_mod_i.assign(&y_i);
        }
            
        //--------------------------------------
        // Lineshape

        // Calculating peak fit and gradient terms 
        self.lineshape.calculate(p);

        // Updating with lineshapes
        self.y_fit += &self.lineshape.y.slice(s![.., ..]);

        //--------------------------------------
        // Baseline

        // Tacking on baseline
        if nb > 0 {
            let pb_r = p[nl .. (nl + nb)].to_vec();
            let pb_r: Array<_, Ix1> = Array::from_shape_vec((nb, ), pb_r).unwrap(); 
            let mut y_fit_r: ArrayViewMut<_, Ix1> = self.y_fit.slice_mut(s![0, ..]);

            y_fit_r += &( self.basis.dot(&pb_r) );

            let pb_i = p[(nl + nb) .. (nl + nb*2)].to_vec();
            let pb_i: Array<_, Ix1> = Array::from_shape_vec((nb, ), pb_i).unwrap(); 
            let mut y_fit_i: ArrayViewMut<_, Ix1> = self.y_fit.slice_mut(s![1, ..]);

            y_fit_i += &( self.basis.dot(&pb_i) );
        }

        //--------------------------------------
        // All gradients

        // Calculating total difference
        self.y_diff = &self.y_mod - &self.y_fit;

        let y_diff_r: ArrayView<_, Ix1> = self.y_diff.slice(s![0, ..]);
        let y_diff_i: ArrayView<_, Ix1> = self.y_diff.slice(s![1, ..]);

        // Lineshape
        for i in 0 .. nl {
            let dydp_r: ArrayView<_, Ix1> = self.lineshape.dydp.slice(s![0, i, ..]);
            let dydp_i: ArrayView<_, Ix1> = self.lineshape.dydp.slice(s![1, i, ..]);

            self.grad[[0, i]] = y_diff_r.dot(&dydp_r);
            self.grad[[1, i]] = y_diff_i.dot(&dydp_i);
        } 

        // Baseline 
        if nb > 0 {
            let y_diff_r: ArrayView<_, Ix2> = self.y_diff.slice(s![0 .. 1, ..]);
            let y_diff_i: ArrayView<_, Ix2> = self.y_diff.slice(s![1 .. 2, ..]);

            let grad_r = y_diff_r.dot(&self.basis);
            let grad_i = y_diff_i.dot(&self.basis);

            for i in 0 .. nb {
                self.grad[[0, nl + i]] = grad_r[[0, i]];
                self.grad[[1, nl + nb + i]] = grad_i[[0, i]];
            } 
        }

        // Phase
        if np > 0 {

            let mut grad_theta_r: ArrayViewMut<_, Ix1> = self.grad_theta.slice_mut(s![0, ..]);
            grad_theta_r.assign(&( -&y_r * &self.sin_theta + &y_i * &self.cos_theta ));
            
            let mut grad_theta_i: ArrayViewMut<_, Ix1> = self.grad_theta.slice_mut(s![1, ..]);
            grad_theta_i.assign(&( -&y_r * &self.cos_theta - &y_i * &self.sin_theta ));

            let grad_theta_r: ArrayView<_, Ix1> = self.grad_theta.slice(s![0, ..]);
            let grad_theta_i: ArrayView<_, Ix1> = self.grad_theta.slice(s![1, ..]);

            // The negative is a hack to counteract the general negative below
            if np == 1 {
                self.grad[[0, nl + nb*2]] = -y_diff_r.dot(&grad_theta_r);
                self.grad[[1, nl + nb*2]] = -y_diff_i.dot(&grad_theta_i);
            } else {
                self.grad[[0, nl + nb*2]] = -y_diff_r.dot(&grad_theta_r);
                self.grad[[1, nl + nb*2]] = -y_diff_i.dot(&grad_theta_i);

                self.grad[[0, nl + nb*2 + 1]] = -y_diff_r.dot(
                    &( &grad_theta_r * &self.lineshape.x )
                );
                self.grad[[1, nl + nb*2 + 1]] = -y_diff_i.dot(
                    &( &grad_theta_i * &self.lineshape.x )
                );
            }

        }

        // Tacking on the -2 scalar term to the gradient
        self.grad *= -2.0;

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
*/
