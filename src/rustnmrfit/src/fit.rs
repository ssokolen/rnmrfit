use ndarray::prelude::*;

use crate::testing::Eval;
use crate::lineshape::Lineshape1D;
use crate::baseline::Baseline1D;
use crate::phase::Phase1D;

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
            grad: Array::zeros((2, nl + 2*nb + np)),
        }
    }

    //--------------------------------------
    pub fn obj(p: &[f64], grad: Option<&mut [f64]>, obj: &mut Fit1D) -> f64 {
        obj.eval(p, grad)
    }

    //--------------------------------------
    pub fn eval(&mut self, p: &[f64], grad: Option<&mut [f64]>) -> f64 {

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


/*
//==============================================================================
// Fit2D 

//------------------------------------------------------------------------------
pub struct Fit2D {

    // x, y spectral data (where the x values are unique)
    x_direct: Array<f64, Ix1>,
    x_indirect: Array<f64, Ix1>,
    y: Array<f64, Ix2>,

    // Mapping of unique x values to y values
    xi_direct: Array<usize, Ix1>,
    xi_indirect: Array<usize, Ix1>,

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
    //grad_theta: Array<f64, Ix2>,

    // Primary output
    pub y_fit: Array<f64, Ix2>,
    pub grad: Array<f64, Ix2>,

}

//------------------------------------------------------------------------------
impl Fit2D {

    //--------------------------------------
    pub fn new(x_direct: Array<f64, Ix1>, x_indirect: Array<f64, Ix1>, y: Array<f64, Ix2>,  
               xi_direct: Array<usize, Ix1>, xi_indirect: Array<usize, Ix1>,
               peak_resonances: Array<usize, Ix1>, peak_dimensions: Array<usize, Ix1>,
               nl: usize, nb: usize, np: usize,
               _: Option<Array<f64, Ix2>>) 
        -> Fit2D {

        let ni = xi_direct.len();

        // Generating vectors of lineshapes
        let mut direct: Vec<Vec<usize>> = Vec::new();
        let mut indirect: Vec<Vec<usize>> = Vec::new();
        
        let mut ir: usize;
        let mut id: usize;
        
        // Looping through
        for i in 0 .. peak_resonances.len() {
            
            ir = peak_resonances[i];
            id = peak_dimensions[i];
            
            // Ensure correct vector length
            while direct.len() <= ir {
                direct.push(Vec::new());
                indirect.push(Vec::new());
            }
            
            // Adding peak index
            if id == 0 {
                direct[ir].push(i * 4);
                direct[ir].push(i * 4 + 1);
                direct[ir].push(i * 4 + 2);
                direct[ir].push(i * 4 + 3);
            } else {
                indirect[ir].push(i * 4);
                indirect[ir].push(i * 4 + 1);
                indirect[ir].push(i * 4 + 2);
                indirect[ir].push(i * 4 + 3);
            }
            
        }
        
        let f_direct = |x: &Vec<usize>| -> LineshapeND {
            let indexes = Array::from_shape_vec((x.len(),), x.to_vec()).unwrap();
            LineshapeND::new(x_direct.len(), ni, indexes)
        };

        let f_indirect = |x: &Vec<usize>| -> LineshapeND {
            let indexes = Array::from_shape_vec((x.len(),), x.to_vec()).unwrap();
            LineshapeND::new(x_indirect.len(), ni, indexes)
        };
        
        // Converting into Arrays
        let direct: Vec<LineshapeND> = direct.iter().map(f_direct).collect();
        let indirect: Vec<LineshapeND> = indirect.iter().map(f_indirect).collect();

        Fit2D {
            x_direct: x_direct,
            x_indirect: x_indirect,
            y: y,
            xi_direct: xi_direct,
            xi_indirect: xi_indirect,

            nl: nl,
            nb: nb,
            np: np,

            direct: direct,
            indirect: indirect,

            //basis: basis,

            //theta: Array::zeros((n, )),
            //sin_theta: Array::zeros((n, )),
            //cos_theta: Array::zeros((n, )),

            y_mod: Array::zeros((4, ni)),
            y_diff: Array::zeros((4, ni)),
            //grad_theta: Array::zeros((2, ni)),

            y_fit: Array::zeros((4, ni)),
            grad: Array::zeros((4, nl + 2*nb + np)),
        }
    }

    //--------------------------------------
    pub fn obj(p: &[f64], grad: Option<&mut [f64]>, obj: &mut Fit2D) -> f64 {
        obj.eval(p, grad)
    }

}

//------------------------------------------------------------------------------
impl Fit2D {

    //--------------------------------------
    pub fn eval(&mut self, p: &[f64], grad: Option<&mut [f64]>) -> f64 {

        let _nl = self.nl;
        let _nb = self.nb;
        let _np = self.np;

        // Initializing fit
        self.y_fit.fill(0.0);

        let _y_rr: ArrayView<_, Ix1> = self.y.slice(s![0, ..]);
        let _y_ri: ArrayView<_, Ix1> = self.y.slice(s![1, ..]);
        let _y_ir: ArrayView<_, Ix1> = self.y.slice(s![2, ..]);
        let _y_ii: ArrayView<_, Ix1> = self.y.slice(s![3, ..]);

        //--------------------------------------
        // Phasing

        /*
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
        */

        self.y_mod.assign(&self.y);
            
        //--------------------------------------
        // Lineshape

        // Calculating peak fit and gradient terms 
        for i in 0 .. self.direct.len() {
            let x_direct = self.x_direct.slice(s![..]);
            let x_indirect = self.x_indirect.slice(s![..]);

            let xi_direct = self.xi_direct.slice(s![..]);
            let xi_indirect = self.xi_indirect.slice(s![..]);

            self.direct[i].eval(x_direct, xi_direct, p);
            self.indirect[i].eval(x_indirect, xi_indirect, p);

            // Updating with lineshapes
            let mut y_fit_rr: ArrayViewMut<_, Ix1> = self.y_fit.slice_mut(s![0, ..]);
            y_fit_rr += &( &self.direct[i].y.slice(s![0, ..]) *
                           &self.direct[i].y.slice(s![0, ..]) );

            let mut y_fit_ri: ArrayViewMut<_, Ix1> = self.y_fit.slice_mut(s![1, ..]);
            y_fit_ri += &( &self.direct[i].y.slice(s![0, ..]) *
                           &self.direct[i].y.slice(s![1, ..]) );

            let mut y_fit_ir: ArrayViewMut<_, Ix1> = self.y_fit.slice_mut(s![2, ..]);
            y_fit_ir += &( &self.direct[i].y.slice(s![1, ..]) *
                           &self.direct[i].y.slice(s![0, ..]) );

            let mut y_fit_ii: ArrayViewMut<_, Ix1> = self.y_fit.slice_mut(s![3, ..]);
            y_fit_ii += &( &self.direct[i].y.slice(s![1, ..]) *
                           &self.direct[i].y.slice(s![1, ..]) );
        }

        //--------------------------------------
        // Baseline

        /*
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
        */

        //--------------------------------------
        // All gradients

        // Calculating total difference
        self.y_diff = &self.y_mod - &self.y_fit;

        let y_diff_rr: ArrayView<_, Ix1> = self.y_diff.slice(s![0, ..]);
        let y_diff_ri: ArrayView<_, Ix1> = self.y_diff.slice(s![1, ..]);
        let y_diff_ir: ArrayView<_, Ix1> = self.y_diff.slice(s![2, ..]);
        let y_diff_ii: ArrayView<_, Ix1> = self.y_diff.slice(s![3, ..]);

        // Looping over direct dimension peaks
        for i in 0 .. self.direct.len() {

            // Direct dimension i is multiplied by indirect dimension i
            let y_indirect_r: ArrayView<_, Ix1> = self.indirect[i].y.slice(s![0, ..]);
            let y_indirect_i: ArrayView<_, Ix1> = self.indirect[i].y.slice(s![1, ..]);
            
            // Index j corresponds to parameter j of np accounted for by each lineshape
            for j in 0 .. self.direct[i].np {

                let dydp_r: ArrayView<_, Ix1> = self.direct[i].dydp.slice(s![0, j, ..]);
                let dydp_i: ArrayView<_, Ix1> = self.direct[i].dydp.slice(s![1, j, ..]);

                // Use index array to figure out which parameter j actually represents
                let k = self.direct[i].indexes[j];

                self.grad[[0, k]] = y_diff_rr.dot(&(&dydp_r * &y_indirect_r));
                self.grad[[1, k]] = y_diff_ri.dot(&(&dydp_r * &y_indirect_i));
                self.grad[[2, k]] = y_diff_ir.dot(&(&dydp_i * &y_indirect_r));
                self.grad[[3, k]] = y_diff_ii.dot(&(&dydp_i * &y_indirect_i));

            }
        }

        // Repeating the same procedure for indirect dimension peaks
        for i in 0 .. self.indirect.len() {

            // Direct dimension i is multiplied by indirect dimension i
            let y_direct_r: ArrayView<_, Ix1> = self.direct[i].y.slice(s![0, ..]);
            let y_direct_i: ArrayView<_, Ix1> = self.direct[i].y.slice(s![1, ..]);
            
            // Index j corresponds to parameter j of np accounted for by each lineshape
            for j in 0 .. self.indirect[i].np {

                let dydp_r: ArrayView<_, Ix1> = self.indirect[i].dydp.slice(s![0, j, ..]);
                let dydp_i: ArrayView<_, Ix1> = self.indirect[i].dydp.slice(s![1, j, ..]);

                // Use index array to figure out which parameter j actually represents
                let k = self.direct[i].indexes[j];

                self.grad[[0, k]] = y_diff_rr.dot(&(&dydp_r * &y_direct_r));
                self.grad[[1, k]] = y_diff_ri.dot(&(&dydp_i * &y_direct_r));
                self.grad[[2, k]] = y_diff_ir.dot(&(&dydp_r * &y_direct_i));
                self.grad[[3, k]] = y_diff_ii.dot(&(&dydp_i * &y_direct_i));

            }
        } 

        /*
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
        */

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
