use ndarray::{prelude::*, ArcArray1};
use itertools::{
    izip,
    Itertools
};

use crate::fitting::common::ComponentEval;


//=============================================================================
// Commonly used types

type YSlice<'a> = ArrayViewMut<'a, f64, Ix2>;
type DySlice<'a> = ArrayViewMut<'a, f64, Ix3>;


//=============================================================================
// General 1D phase correction

pub struct Phase1D {

    // Reference to x values for convenience
    x: ArcArray1<f64>,

    // y data used to generate phase
    y: Array2<f64>,

    // Phase function based on parameter number
    _np: usize,
    f: fn (&mut Phase1D, &[f64], y: YSlice, dy: DySlice),
}


//--------------------------------------
impl Phase1D {

    //--------------------------------------
    pub fn new(x: ArcArray1<f64>, y: Array2<f64>, np: usize) -> Phase1D {

        // Choose f based on number of parameters
        let f = match np {
            0 => Phase1D::eval_0,
            1 => Phase1D::eval_1,
            2 => Phase1D::eval_2,
            _ => panic!("Nonlinear phase correction is not supported.")
        };

        Phase1D {
            x: x,
            y: y,
            _np: np,
            f: f,
        }
    }


    //--------------------------------------
    pub fn eval_0(&mut self, _p: &[f64], mut y: YSlice, _dy: DySlice) {       

        // If no phasing terms, output stored y
        y.assign(&self.y);

    }


    //--------------------------------------
    pub fn eval_1(&mut self, p: &[f64], y: YSlice, dy: DySlice) {       

        // With a single phase term, generate single value intermediates
        let theta = p[0];
        let s = theta.sin();
        let c = theta.cos();

        let yo = self.y.slice(s![..,..]);

        let iterator = izip!(
            yo.into_iter().tuples::<(_,_)>(),
            y.into_iter().tuples::<(_,_)>(),
            dy.into_iter().tuples::<(_,_)>(),
        );

        for (yo, y, dy) in iterator {

            let (yo_r, yo_i) = yo;
            let (y_r, y_i) = y;
            let (dy_r, dy_i) = dy;

            *y_r =  yo_r*c + yo_i*s;
            *y_i = -yo_r*s + yo_i*c;

            *dy_r = -yo_r*s + yo_i*c;
            *dy_i = -yo_r*c - yo_i*s;
        }
    }


    //--------------------------------------
    pub fn eval_2(&mut self, p: &[f64], y: YSlice, mut dy: DySlice) {       

        // With two phase terms, intermediates must be calculated within zip
        let theta_0 = p[0];
        let theta_1 = p[1];

        let x = self.x.slice(s![..]);
        let yo = self.y.slice(s![..,..]);

        let (dy0, dy1) = dy.multi_slice_mut((
                s![0, .., ..],
                s![1, .., ..]
        ));

        let iterator = izip!(
            x.into_iter(),
            yo.into_iter().tuples::<(_,_)>(),
            y.into_iter().tuples::<(_,_)>(),
            dy0.into_iter().tuples::<(_, _)>(),
            dy1.into_iter().tuples::<(_, _)>(),
        );

        for (x, yo, y, dy0, dy1) in iterator {

            let s = (theta_0 + theta_1*x).sin(); 
            let c = (theta_0 + theta_1*x).cos(); 

            let (yo_r, yo_i) = yo;
            let (y_r, y_i) = y;
            let (dy0_r, dy0_i) = dy0;
            let (dy1_r, dy1_i) = dy1;

            *y_r =  yo_r*c + yo_i*s;
            *y_i = -yo_r*s + yo_i*c;

            *dy0_r = -yo_r*s + yo_i*c;
            *dy0_i = -yo_r*c - yo_i*s;

            *dy1_r = *dy0_r * x;
            *dy1_i = *dy0_i * x;
        }
    }
}


impl ComponentEval for Phase1D {

    fn nd(&self) -> usize {
        2 
    }

    fn n(&self) -> usize {
        self.x.len()
    }


    //--------------------------------------
    fn eval(&mut self, p: &[f64], y: YSlice, dy: DySlice) {       

        (self.f)(self, p, y, dy);
    }
}


//=============================================================================
// General 2D phase correction
#[allow(dead_code)]
pub struct Phase2D {

    // Reference to x values for convenience
    x_direct: ArcArray1<f64>,
    x_indirect: ArcArray1<f64>,

    // y data used to generate phase
    y: Array2<f64>,

    // Phase function based on parameter number
    np: usize,
    f: fn (&mut Phase2D, &[f64], y: YSlice, dy: DySlice),
}


//--------------------------------------
#[allow(dead_code)]
impl Phase2D {

    //--------------------------------------
    pub fn new(x_direct: ArcArray1<f64>, x_indirect: ArcArray1<f64>, 
               y: Array2<f64>, np: usize) -> Phase2D {

        // Choose f based on number of parameters
        let f = match np {
            0 => Phase2D::eval_0,
            2 => Phase2D::eval_2,
            3 => Phase2D::eval_3,
            _ => {
                panic!(
                    "2D phase correction may only use 2 terms \
                     (0 order in two dimensions) or 3 (for 1st order)."
                )
            }
        };

        Phase2D {
            x_direct: x_direct,
            x_indirect: x_indirect,
            y: y,
            np: np,
            f: f,
        }
    }


    //--------------------------------------
    pub fn eval_0(&mut self, _p: &[f64], mut y: YSlice, _dy: DySlice) {       

        // If no phasing terms, output stored y
        y.assign(&self.y);
    }


    //--------------------------------------
    pub fn eval_2(&mut self, p: &[f64], y: YSlice, mut dy: DySlice) {       

        // With a single phase term, generate single value intermediates
        let theta = p[0];
        let sd = theta.sin();
        let cd = theta.cos();

        let theta = p[1];
        let si = theta.sin();
        let ci = theta.cos();

        let yo = self.y.slice(s![..,..]);

        let (dy0, dy1) = dy.multi_slice_mut((
                s![0, .., ..],
                s![1, .., ..]
        ));

        let iterator = izip!(
            yo.into_iter().tuples::<(_,_,_,_)>(),
            y.into_iter().tuples::<(_,_,_,_)>(),
            dy0.into_iter().tuples::<(_,_,_,_)>(),
            dy1.into_iter().tuples::<(_,_,_,_)>(),
        );

        for (yo, y, dy0, dy1) in iterator {

            let (yo_rr, yo_ri, yo_ir, yo_ii) = yo;
            let (y_rr, y_ri, y_ir, y_ii) = y;
            let (dy0_rr, dy0_ri, dy0_ir, dy0_ii) = dy0;
            let (dy1_rr, dy1_ri, dy1_ir, dy1_ii) = dy1;

            *y_rr =  yo_rr*cd*ci + yo_ri*cd*si + yo_ir*sd*ci + yo_ii*sd*si;
            *y_ri = -yo_rr*cd*si + yo_ri*cd*ci - yo_ir*sd*si + yo_ii*sd*ci;
            *y_ir = -yo_rr*sd*ci - yo_ri*sd*si + yo_ir*cd*ci + yo_ii*cd*si;
            *y_ii =  yo_rr*sd*si - yo_ri*sd*ci - yo_ir*cd*si + yo_ii*cd*ci;

            *dy0_rr = -yo_rr*sd*ci - yo_ri*sd*si + yo_ir*cd*ci + yo_ii*cd*si;
            *dy0_ri =  yo_rr*sd*si - yo_ri*sd*ci - yo_ir*cd*si + yo_ii*cd*ci;
            *dy0_ir = -yo_rr*cd*ci - yo_ri*cd*si - yo_ir*sd*ci - yo_ii*sd*si;
            *dy0_ii =  yo_rr*cd*si - yo_ri*cd*ci + yo_ir*sd*si - yo_ii*sd*ci;

            *dy1_rr = -yo_rr*cd*si + yo_ri*cd*ci - yo_ir*sd*si + yo_ii*sd*ci;
            *dy1_ri = -yo_rr*cd*ci - yo_ri*cd*si - yo_ir*sd*ci - yo_ii*sd*si;
            *dy1_ir =  yo_rr*sd*si - yo_ri*sd*ci - yo_ir*cd*si + yo_ii*cd*ci;
            *dy1_ii =  yo_rr*sd*ci + yo_ri*sd*si - yo_ir*cd*ci - yo_ii*cd*si;
        }
    }


    //--------------------------------------
    pub fn eval_3(&mut self, p: &[f64], y: YSlice, mut dy: DySlice) {       

        // With multiple terms, intermediates must be calculated within zip
        let theta_0 = p[0];
        let theta_d = p[1];
        let theta_i = p[2];

        let xd = self.x_direct.slice(s![..]);
        let xi = self.x_indirect.slice(s![..]);

        let yo = self.y.slice(s![..,..]);

        let (dy0, dy1, dy2) = dy.multi_slice_mut((
                s![0, .., ..],
                s![1, .., ..],
                s![2, .., ..],
        ));

        let iterator = izip!(
            xd.into_iter(),
            xi.into_iter(),
            yo.into_iter().tuples::<(_,_,_,_)>(),
            y.into_iter().tuples::<(_,_,_,_)>(),
            dy0.into_iter().tuples::<(_,_,_,_)>(),
            dy1.into_iter().tuples::<(_,_,_,_)>(),
            dy2.into_iter().tuples::<(_,_,_,_)>(),
        );

        for (x_d, x_i, yo, y, dy0, dy1, dy2) in iterator {

            let sd = (theta_0 + theta_d*x_d).sin(); 
            let cd = (theta_0 + theta_d*x_d).cos(); 
            let si = (theta_0 + theta_i*x_i).sin(); 
            let ci = (theta_0 + theta_i*x_i).cos(); 

            let (yo_rr, yo_ri, yo_ir, yo_ii) = yo;
            let (y_rr, y_ri, y_ir, y_ii) = y;
            let (dy0_rr, dy0_ri, dy0_ir, dy0_ii) = dy0;
            let (dy1_rr, dy1_ri, dy1_ir, dy1_ii) = dy1;
            let (dy2_rr, dy2_ri, dy2_ir, dy2_ii) = dy2;

            *y_rr =  yo_rr*cd*ci + yo_ri*cd*si + yo_ir*sd*ci + yo_ii*sd*si;
            *y_ri = -yo_rr*cd*si + yo_ri*cd*ci - yo_ir*sd*si + yo_ii*sd*ci;
            *y_ir = -yo_rr*sd*ci - yo_ri*sd*si + yo_ir*cd*ci + yo_ii*cd*si;  
            *y_ii =  yo_rr*sd*si - yo_ri*sd*ci - yo_ir*cd*si + yo_ii*cd*ci;

            *dy0_rr =  yo_rr*(-sd*ci - cd*si) + yo_ri*(-sd*si + cd*ci) 
                     + yo_ir*( cd*ci - sd*si) + yo_ii*( cd*si + sd*ci);
            *dy0_ri = -yo_rr*(-sd*si + cd*ci) + yo_ri*(-sd*ci - cd*si) 
                     - yo_ir*( cd*si + sd*ci) + yo_ii*( cd*ci - sd*si);
            *dy0_ir = -yo_rr*( cd*ci - sd*si) - yo_ri*( cd*si + sd*ci) 
                     + yo_ir*(-sd*ci - cd*si) + yo_ii*(-sd*si + cd*ci);
            *dy0_ii =  yo_rr*( cd*si + sd*ci) - yo_ri*( cd*ci - sd*si) 
                     - yo_ir*(-sd*si + cd*ci) + yo_ii*(-sd*ci - cd*si);

            *dy1_rr = (-yo_rr*sd*ci - yo_ri*sd*si 
                      + yo_ir*cd*ci + yo_ii*cd*si)*x_d;
            *dy1_ri = ( yo_rr*sd*si - yo_ri*sd*ci 
                      - yo_ir*cd*si + yo_ii*cd*ci)*x_d;
            *dy1_ir = (-yo_rr*cd*ci - yo_ri*cd*si 
                      - yo_ir*sd*ci - yo_ii*sd*si)*x_d;               
            *dy1_ii = ( yo_rr*cd*si - yo_ri*cd*ci 
                      + yo_ir*sd*si - yo_ii*sd*ci)*x_d;

            *dy2_rr = (-yo_rr*cd*si + yo_ri*cd*ci 
                      - yo_ir*sd*si + yo_ii*sd*ci)*x_i;
            *dy2_ri = (-yo_rr*cd*ci - yo_ri*cd*si 
                      - yo_ir*sd*ci - yo_ii*sd*si)*x_i;
            *dy2_ir = ( yo_rr*sd*si - yo_ri*sd*ci 
                      - yo_ir*cd*si + yo_ii*cd*ci)*x_i;               
            *dy2_ii = ( yo_rr*sd*ci + yo_ri*sd*si
                      - yo_ir*cd*ci - yo_ii*cd*si)*x_i;

        }
    }
}


impl ComponentEval for Phase2D {

    fn nd(&self) -> usize {
       4 
    }

    fn n(&self) -> usize {
        self.x_direct.len()
    }


    //--------------------------------------
    fn eval(&mut self, p: &[f64], y: YSlice, dy: DySlice) {       

        (self.f)(self, p, y, dy);

    }
}

/*
#[cfg(test)]
mod tests {
    
    use ndarray::prelude::*;
    use std::iter;

    use crate::fitting::{
        common::ComponentEval,
        lineshape::{Lineshape1D, Lineshape2D},
    };

    use super::{Phase1D, Phase2D};

    #[test]
    fn phase_1d_gradient() {

        let n = 20;
        let x: Array1<f64> = Array::linspace(0.0, 1.0, n);
        let x = x.into_shared();
        
        let p = vec![
            // peak
            0.5, 0.5, 0.5, 0.5, 
            // peak
            0.3, 0.3, 0.3, 0.3
        ];
        let np = p.len();

        // Initialize work arrays
        let mut y: Array2<f64> = Array::zeros((n, 2));
        let mut dy: Array3<f64> = Array::zeros((np, n, 2));

        // Generate peaks
        let y_slice = y.slice_mut(s![.., ..]);
        let dy_slice = dy.slice_mut(s![.., .., ..]);

        let mut lineshape = Lineshape1D::new(x.clone(), np, 1e-10);
        lineshape.eval(&p, y_slice, dy_slice);

        // Then use a separate set of parameters for the phase
        let p = vec![0.15, 0.05];

        println!("0-order phase first");
        let mut phase = Phase1D::new(x.clone(), y.clone(), 1);
        phase.grad_test(&p[..1], 1e-8);

        println!("1st-order phase second");
        let mut phase = Phase1D::new(x.clone(), y.clone(), 2);
        phase.grad_test(&p, 1e-8);

    }

    #[test]
    fn phase_2d_gradient() {

        let n = 20;
        let x: Array1<f64> = Array::linspace(0.0, 1.0, n);
        
        let p = vec![
            //peak
            0.5, 0.5, 0.5, 0.5, 
            //peak
            0.3, 0.3, 0.3, 0.3
        ];
        let np = p.len();

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
        let mut dy: Array3<f64> = Array::zeros((np, n*n, 4));

        // Generate peaks
        let y_slice = y.slice_mut(s![.., ..]);
        let dy_slice = dy.slice_mut(s![.., .., ..]);

        let mut lineshape = Lineshape2D::new(
            x1.clone(), x2.clone(), parameter_map, 1e-10
        );
        lineshape.eval(&p, y_slice, dy_slice);

        // Then use a separate set of parameters for the phase
        let p = vec![0.15, 0.1, 0.05];

        println!("0-order phase first");
        let mut phase = Phase2D::new(x1.clone(), x2.clone(), y.clone(), 2);
        phase.grad_test(&p[..2], 1e-8);

        println!("1st-order phase second");
        let mut phase = Phase2D::new(x1.clone(), x2.clone(), y.clone(), 3);
        phase.grad_test(&p, 1e-8);
    }
}

*/
