use core::ops::Range;
use ndarray::prelude::*;

use crate::fitting::{
    common::{ComponentEval, FitEval},
    lineshape::{Lineshape1D},//, Lineshape2D},
    baseline::{Baseline1D},
    phase::{Phase1D}//, Phase2D},
};

//=============================================================================
// Fit1D 

//-----------------------------------------------------------------------------
pub struct Fit1D {

    // Number of data points
    n: usize,

    // Number of parameters (p) dedicated to each component of fit
    nl: usize,
    nb: usize,
    np: usize,

    // Fit components
    lineshape: Lineshape1D,
    baseline: Baseline1D,
    phase: Phase1D,

    // Work arrays 
    y_diff: Array2<f64>,
    y: Array2<f64>,
    dy: Array3<f64>,
}

//-----------------------------------------------------------------------------
impl Fit1D {

    //--------------------------------------
    pub fn new(
            x: Array1<f64>, 
            y: Array2<f64>,
            ranges: Vec<Range<usize>>,
            assignments: Vec<usize>,
            basis: Array2<f64>,
            nl: usize, 
            nb: usize, 
            np: usize, 
            tol: f64
        ) -> Fit1D {

        let n = x.len();
        let x = x.into_shared();

        Fit1D {
            n: n,

            nl: nl,
            nb: nb,
            np: np,

            lineshape: Lineshape1D::new(x.clone(), ranges, assignments, tol),
            baseline: Baseline1D::new(basis),
            phase: Phase1D::new(x.clone(), y, np),
            y_diff: Array::zeros((n, 2)),
            y: Array::zeros((n, 2)),
            dy: Array::zeros((nl + nb*2 + np, n, 2)),
        }
    }

    //--------------------------------------
    pub fn obj(p: &[f64], grad: Option<&mut [f64]>, obj: &mut Fit1D) -> f64 {
        obj.eval(p, grad)
    }

}

//--------------------------------------
impl FitEval for Fit1D {

    //--------------------------------------
    fn eval(&mut self, p: &[f64], grad: Option<&mut [f64]>) -> f64 {

        let n = self.n;

        let nl = self.nl;
        let nb = self.nb;
        let np = self.np;

        // Starting y_diff from phased y value
        let range = (nl + nb*2) .. (nl + nb*2 + np);
        let y_obs = self.y_diff.slice_mut(s![.., ..]);
        let dy = self.dy.slice_mut(s![range.clone(), .., ..]);
        self.phase.eval(&p[range], y_obs, dy);

        // Lineshape
        let range = 0 .. nl;
        let y = self.y.slice_mut(s![.., ..]);
        let dy = self.dy.slice_mut(s![range.clone(), .., ..]);
        self.lineshape.eval(&p[range], y, dy);

        // Updating y_diff with fit
        self.y_diff -= &self.y;

        // Baseline
        let range = nl .. (nl + nb*2);
        let y = self.y.slice_mut(s![.., ..]);
        let dy = self.dy.slice_mut(s![range.clone(), .., ..]);
        self.baseline.eval(&p[range], y, dy);

        // Updating y_diff with baseline
        self.y_diff -= &self.y;

        // Overall d(y_diff^2)/dp is just the dot product of y_diff and dy/dp
        // (but the dot product requires matrix reshaping)
        let dy = self.dy.slice(s![..,..,..])
            .into_shape((nl + nb*2 + np, n*2, ))
            .unwrap();

        let y_diff = self.y_diff.slice(s![..,..])
            .into_shape((n*2, ))
            .unwrap();

        let mut grad_mat = dy.dot(&y_diff);

        // Lineshape and baseline terms are multiplied by -2 (from (y-y_fit)^2)
        // and phase terms are multiplies by 2.
        let mut terms = grad_mat.slice_mut(s![0 .. (nl + nb*2)]);
        terms *= -2.0;

        let mut terms = grad_mat.slice_mut(
            s![(nl + nb*2) .. (nl + nb*2 + np)]
        );
        terms *= 2.0;

        // Once gradient is calculated, y_diff is no longer required and
        // can be squared
        for y in self.y_diff.iter_mut() {
            *y = *y * *y;
        }

        // Filling in output
        if let Some(grad) = grad {
            grad.copy_from_slice(&(grad_mat.to_vec())[..]);
        }
        
        self.y_diff.sum()
    }
}


//=============================================================================
// Fit2D 

/*
//-----------------------------------------------------------------------------
#[allow(dead_code)]
pub struct Fit2D {

    // Number of data points
    n: usize,

    // Number of parameters (p) dedicated to each component of fit
    nl: usize,
    nb: usize,
    np: usize,

    // Fit components
    lineshape: Lineshape2D,
    //baseline: Baseline2D,
    phase: Phase2D,

    // Work arrays
    y_diff: Array2<f64>,
    y: Array2<f64>,
    dy: Array3<f64>,
}

//-----------------------------------------------------------------------------
#[allow(dead_code)]
impl Fit2D {

    //--------------------------------------
    pub fn new(
            x_direct: Array1<f64>, 
            x_indirect: Array1<f64>,
            y: Array2<f64>,
            parameter_map: Vec<(usize, usize)>,
            _knots: Array1<f64>,
            nl: usize,
            nb: usize,
            np: usize,
            tol: f64
        ) -> Fit2D {

        let n = x_direct.len();
        let x_direct = x_direct.into_shared();
        let x_indirect = x_indirect.into_shared();

        Fit2D {
            n: n,

            nl: nl,
            nb: nb,
            np: np,

            lineshape: Lineshape2D::new(
                x_direct.clone(), x_indirect.clone(), parameter_map, tol
            ),
            //baseline: Baseline2D::new(
            //    x_direct.clone(), x_indirect.clone(), knots, nb
            //),
            phase: Phase2D::new(x_direct.clone(), x_indirect.clone(), y, np),

            y_diff: Array::zeros((n, 4)),
            y: Array::zeros((n, 4)),
            dy: Array::zeros((nl + nb*4 + np, n, 4)),
        }
    }


    //--------------------------------------
    pub fn obj(p: &[f64], grad: Option<&mut [f64]>, obj: &mut Fit2D) -> f64 {
        obj.eval(p, grad)
    }
}


//--------------------------------------
impl FitEval for Fit2D {

    //--------------------------------------
    fn eval(&mut self, p: &[f64], grad: Option<&mut [f64]>) -> f64 {

        let n = self.n;

        let nl = self.nl;
        let nb = self.nb;
        let np = self.np;

        // Starting y_diff from phased y value
        let range = (nl + nb*4) .. (nl + nb*4 + np);
        let y_obs = self.y_diff.slice_mut(s![.., ..]);
        let dy = self.dy.slice_mut(s![range.clone(), .., ..]);
        self.phase.eval(&p[range], y_obs, dy);

        // Initializing fit
        self.y.fill(0.0);

        // Lineshape
        let range = 0 .. nl;
        let y = self.y.slice_mut(s![.., ..]);
        let dy = self.dy.slice_mut(s![range.clone(), .., ..]);
        self.lineshape.eval(&p[range], y, dy);

        // Updating y_diff with fit
        self.y_diff -= &self.y;

        // Baseline
        //let range = nl .. (nl + nb*2);
        //let y = self.y.slice_mut(s![.., ..]);
        //let dy = self.dy.slice_mut(s![range.clone(), .., ..]);
        //self.baseline.eval(&p[range], &y, &dy);

        // Updating y_diff with baseline
        //self.y_diff -= &self.y;
            
        // Overall d(y_diff^2)/dp is just the dot product of y_diff and dy/dp
        // (but the dot product requires matrix reshaping)
        let dy = self.dy.slice(s![..,..,..])
            .into_shape((nl + nb*4 + np, n*4, ))
            .unwrap();

        let y_diff = self.y_diff.slice(s![..,..])
            .into_shape((n*4, ))
            .unwrap();

        let mut grad_mat = dy.dot(&y_diff);

        // Lineshape and baseline terms are multiplied by -2 (from (y-y_fit)^2)
        // and phase terms are multiplies by 2.
        let mut terms = grad_mat.slice_mut(s![0 .. (nl + nb*4)]);
        terms *= -2.0;

        let mut terms = grad_mat.slice_mut(
            s![(nl + nb*4) .. (nl + nb*4 + np)]
        );
        terms *= 2.0;

        // Once gradient is calculated, y_diff is no longer required and
        // can be squared
        for y in self.y_diff.iter_mut() {
            *y = *y * *y;
        }

        // Filling in output
        if let Some(grad) = grad {
            grad.copy_from_slice(&(grad_mat.to_vec())[..]);
        }
        
        self.y_diff.sum()
    }
}
*/


//=============================================================================
// Unit tests

#[cfg(test)]
mod tests {

    use ndarray::prelude::*;
    use ndarray_rand::{
        RandomExt,
        rand::SeedableRng,
        rand_distr::Normal,
    };
    use rand_isaac::isaac64::Isaac64Rng;
    use std::iter;

    use crate::fitting::{
        common::{ComponentEval, FitEval},
        lineshape::{Lineshape1D},//, Lineshape2D},
    };

    /*
    //--------------------------------------
    #[test]
    fn fit_1d_gradient() {

        let n = 20;
        let x: Array1<f64> = Array::linspace(0.0, 1.0, n);
        let x = x.into_shared();

        let p = vec![
            // peak
            0.5, 0.5, 0.5, 0.5, 
            // peak
            0.3, 0.3, 0.3, 0.3, 
            // baseline         
            0.1, 0.2, 0.3, 0.4, 0.5, 
            0.1, 0.2, 0.3, 0.4, 0.5, 
            // phase         
            0.1, 0.05
        ];

        let nl = 8;
        let nb = 5;
        let np = 2;

        // Double checking basic math
        assert!((nl + nb*2 + np) == p.len(), "Error setting up test.");

        // Initialize work arrays
        let mut y: Array2<f64> = Array::zeros((n, 2));
        let mut dy: Array3<f64> = Array::zeros((nl, n, 2));

        // Generate peaks
        let y_slice = y.slice_mut(s![.., ..]);
        let dy_slice = dy.slice_mut(s![.., .., ..]);

        let mut lineshape = Lineshape1D::new(x.clone(), np, 1e-10);
        lineshape.eval(&p[..nl], y_slice, dy_slice);

        // Adding some noise for a bit of realism
        let seed = 1111;
        let mut rng = Isaac64Rng::seed_from_u64(seed);
        let noise = Array::random_using(
            (n, 2), Normal::new(0., 0.01).unwrap(), &mut rng
        );
        y += &noise;

        // Basic boundary knots
        let knots: Array1<f64> = Array::from_shape_vec(
            (2,), vec![0.0, 1.0]
        ).unwrap();

        // Initialize fit
        let mut fit = super::Fit1D::new(
            x.to_owned(), y, knots, nl, nb, np, 1e-10
        );

        // Summing up errors increases threshold
        fit.grad_test(&p, 5e-4);
    }


    //--------------------------------------
    #[test]
    fn fit_2d_gradient() {

        let n = 20;
        let x: Array1<f64> = Array::linspace(0.0, 1.0, n);

        let p = vec![
            // peak
            0.5, 0.5, 0.5, 0.5, 
            // peak
            0.3, 0.3, 0.3, 0.3, 
            // baseline 
            0.1, 0.2, 0.3, 0.4, 
            0.1, 0.2, 0.3, 0.4, 
            0.1, 0.2, 0.3, 0.4, 
            0.1, 0.2, 0.3, 0.4,
            //phase 
            0.1, 0.05, 0.05
        ];

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

        // Place 2nd peak in indirect dimension
        let parameter_map: Vec<(usize, usize)> = iter::repeat((1,1)).take(4)
            .chain(iter::repeat((1,2)).take(4))
            .collect();

        // Initialize work arrays
        let mut y: Array2<f64> = Array::zeros((n*n, 4));
        let mut dy: Array3<f64> = Array::zeros((nl, n*n, 4));

        // Generate peaks
        let y_slice = y.slice_mut(s![.., ..]);
        let dy_slice = dy.slice_mut(s![.., .., ..]);

        let mut lineshape = Lineshape2D::new(
            x1.clone(), x2.clone(), parameter_map.clone(), 1e-10
        );
        lineshape.eval(&p[..nl], y_slice, dy_slice);

        // Adding some noise for a bit of realism
        let seed = 1111;
        let mut rng = Isaac64Rng::seed_from_u64(seed);
        let noise = Array::random_using(
            (n*n, 4), Normal::new(0., 0.01).unwrap(), &mut rng
        );
        y += &noise;

        // Basic boundary knots
        let knots: Array1<f64> = Array::from_shape_vec(
            (2,), vec![0.0, 1.0]
        ).unwrap();

        // Initialize fit
        let mut fit = super::Fit2D::new(
            x1.to_owned(), x2.to_owned(), y,
            parameter_map, knots, nl, nb, np, 1e-10
        );

        // Summing up errors increases threshold
        fit.grad_test(&p, 5e-4);
    }
    */
}
