use ndarray::{prelude::*, ArcArray1, Zip};

use crate::common::NMRFitComponent;

//==============================================================================
// General 1D phase correction

pub struct Phase1D {

    // Reference to x values for convenience
    x: ArcArray1<f64>,

    // Initial values
    pub y_in: Array2<f64>,

    // Intermediate values
    grad: Array2<f64>,

    // Phase function
    f: fn (&mut Phase1D, &[f64]),

    // Final output
    pub y: Array2<f64>,
    pub dydp: Array3<f64>,

}

//--------------------------------------
impl Phase1D {

    //--------------------------------------
    pub fn new(x: ArcArray1<f64>, y: Array2<f64>, np: usize) -> Phase1D {

        let n = x.len();
        
        // Choose f based on number of parameters
        let f = match np {
            0 => Phase1D::eval_0,
            1 => Phase1D::eval_1,
            2 => Phase1D::eval_2,
            _ => panic!("Nonlinear phase correction is not supported.")
        };

        Phase1D {
            x: x,
            y_in: y.clone(),
            
            grad: Array::zeros((2, n)),

            f: f,

            y: y,
            dydp: Array::zeros((2, np, n)),
        }
    }

    //--------------------------------------
    pub fn eval_0(&mut self, _p: &[f64]) {       

        // Nothing to do -- if there are no phasing terms, then y is unchanged

    }

    //--------------------------------------
    pub fn eval_1(&mut self, p: &[f64]) {       

        // With a single phase term, generate single value intermediates
        let theta = p[0];
        let s = theta.sin();
        let c = theta.cos();

        Zip::from(self.y.lanes_mut(Axis(0)))
            .and(self.dydp.slice_mut(s![.., 0, ..]).lanes_mut(Axis(0)))
            .and(self.y_in.lanes(Axis(0)))
            .apply(|mut y, mut dy, yo| {
                let yo_r = yo[0];
                let yo_i = yo[1];

                y[0] =   yo_r*c + yo_i*s;
                y[1] =  -yo_r*s + yo_i*c;

                dy[0] = -yo_r*s + yo_i*c;
                dy[1] = -yo_r*c - yo_i*s;
            });
    }

    //--------------------------------------
    pub fn eval_2(&mut self, p: &[f64]) {       

        // With two phase terms, intermediates must be calculated within zip
        let theta_0 = p[0];
        let theta_1 = p[1];

        // Calculating derivative with respect to (theta_0 + theta_1*x) first
        Zip::from(self.y.lanes_mut(Axis(0)))
            .and(self.grad.lanes_mut(Axis(0)))
            .and(self.y_in.lanes(Axis(0)))
            .and(&self.x)
            .apply(|mut y, mut dy, yo, &x| {
                let s = (theta_0 + theta_1*x).sin(); 
                let c = (theta_0 + theta_1*x).cos(); 
                
                let yo_r = yo[0];
                let yo_i = yo[1];

                y[0] =   yo_r*c + yo_i*s;
                y[1] =  -yo_r*s + yo_i*c;

                dy[0] = -yo_r*s + yo_i*c;
                dy[1] = -yo_r*c - yo_i*s;
            });
    
        // And then applying chain rule for separate theta_0 and theta_1 terms
        let mut dydp_0 = self.dydp.slice_mut(s![.., 0, ..]);       
        dydp_0.assign(&self.grad);

        Zip::from(&mut self.dydp.slice_mut(s![.., 1, ..]))
            .and(&self.grad)
            .and_broadcast(&self.x)
            .apply(|dy1, &dy, &x| {
                *dy1 = dy*x;
        });

    }

}

impl NMRFitComponent for Phase1D {

    fn get_y(&self) -> Array2<f64> {
        self.y.clone()
    }

    fn get_dydp(&self) -> Array3<f64> {
        self.dydp.clone()
    }

    //--------------------------------------
    fn eval(&mut self, p: &[f64]) {       

        (self.f)(self, p);

    }
}

//==============================================================================
// General 2D phase correction

pub struct Phase2D {

    // Reference to x values for convenience
    x_direct: ArcArray1<f64>,
    x_indirect: ArcArray1<f64>,

    // Initial values
    pub y_in: Array2<f64>,

    // Phase function
    f: fn (&mut Phase2D, &[f64]),

    // Final output
    pub y: Array2<f64>,
    pub dydp: Array3<f64>,

}

//--------------------------------------
impl Phase2D {

    //--------------------------------------
    pub fn new(x_direct: ArcArray1<f64>, x_indirect: ArcArray1<f64>, 
               y: Array2<f64>, np: usize) -> Phase2D {

        let n = x_direct.len();
        
        // Choose f based on number of parameters
        let f = match np {
            0 => Phase2D::eval_0,
            2 => Phase2D::eval_2,
            3 => Phase2D::eval_3,
            _ => panic!("2D phase correction may only use 2 terms (0 order in two dimensions) or 3 (for 1st order).")
        };

        Phase2D {
            x_direct: x_direct,
            x_indirect: x_indirect,
            y_in: y.clone(),
            
            f: f,

            y: y,
            dydp: Array::zeros((4, np, n)),
        }
    }

    //--------------------------------------
    pub fn eval_0(&mut self, _p: &[f64]) {       

        // Nothing to do -- if there are no phasing terms, then y is unchanged

    }

    //--------------------------------------
    pub fn eval_2(&mut self, p: &[f64]) {       

        // With a single phase term per dimension, generate single value intermediates
        let theta = p[0];
        let sd = theta.sin();
        let cd = theta.cos();

        let theta = p[1];
        let si = theta.sin();
        let ci = theta.cos();

        Zip::from(self.y.lanes_mut(Axis(0)))
            .and(self.dydp.axis_iter_mut(Axis(2)))
            .and(self.y_in.lanes(Axis(0)))
            .apply(|mut y, mut dy, yo| {
                let yo_rr = yo[0];
                let yo_ri = yo[1];
                let yo_ir = yo[2];
                let yo_ii = yo[3];

                y[0] =  yo_rr*cd*ci + yo_ri*cd*si + yo_ir*sd*ci + yo_ii*sd*si;
                y[1] = -yo_rr*cd*si + yo_ri*cd*ci - yo_ir*sd*si + yo_ii*sd*ci;
                y[2] = -yo_rr*sd*ci - yo_ri*sd*si + yo_ir*cd*ci + yo_ii*cd*si;               
                y[3] =  yo_rr*sd*si - yo_ri*sd*ci - yo_ir*cd*si + yo_ii*cd*ci;

                dy[[0,0]] = -yo_rr*sd*ci - yo_ri*sd*si + yo_ir*cd*ci + yo_ii*cd*si;
                dy[[1,0]] =  yo_rr*sd*si - yo_ri*sd*ci - yo_ir*cd*si + yo_ii*cd*ci;
                dy[[2,0]] = -yo_rr*cd*ci - yo_ri*cd*si - yo_ir*sd*ci - yo_ii*sd*si;               
                dy[[3,0]] =  yo_rr*cd*si - yo_ri*cd*ci + yo_ir*sd*si - yo_ii*sd*ci;

                dy[[0,1]] = -yo_rr*cd*si + yo_ri*cd*ci - yo_ir*sd*si + yo_ii*sd*ci;
                dy[[1,1]] = -yo_rr*cd*ci - yo_ri*cd*si - yo_ir*sd*ci - yo_ii*sd*si;
                dy[[2,1]] =  yo_rr*sd*si - yo_ri*sd*ci - yo_ir*cd*si + yo_ii*cd*ci;               
                dy[[3,1]] =  yo_rr*sd*ci + yo_ri*sd*si - yo_ir*cd*ci - yo_ii*cd*si;

            });

    }

    //--------------------------------------
    pub fn eval_3(&mut self, p: &[f64]) {       

        // With multiple phase terms, intermediates must be calculated within zip
        let theta_0 = p[0];
        let theta_d = p[1];
        let theta_i = p[2];

        Zip::from(self.y.lanes_mut(Axis(0)))
            .and(self.dydp.axis_iter_mut(Axis(2)))
            .and(self.y_in.lanes(Axis(0)))
            .and(&self.x_direct)
            .and(&self.x_indirect)
            .apply(|mut y, mut dy, yo, &x_d, &x_i| {
                let sd = (theta_0 + theta_d*x_d).sin(); 
                let cd = (theta_0 + theta_d*x_d).cos(); 
                let si = (theta_0 + theta_i*x_i).sin(); 
                let ci = (theta_0 + theta_i*x_i).cos(); 

                let yo_rr = yo[0];
                let yo_ri = yo[1];
                let yo_ir = yo[2];
                let yo_ii = yo[3];

                y[0] =  yo_rr*cd*ci + yo_ri*cd*si + yo_ir*sd*ci + yo_ii*sd*si;
                y[1] = -yo_rr*cd*si + yo_ri*cd*ci - yo_ir*sd*si + yo_ii*sd*ci;
                y[2] = -yo_rr*sd*ci - yo_ri*sd*si + yo_ir*cd*ci + yo_ii*cd*si;               
                y[3] =  yo_rr*sd*si - yo_ri*sd*ci - yo_ir*cd*si + yo_ii*cd*ci;

                dy[[0,0]] =  yo_rr*(-sd*ci - cd*si) + yo_ri*(-sd*si + cd*ci) 
                           + yo_ir*( cd*ci - sd*si) + yo_ii*( cd*si + sd*ci);
                dy[[1,0]] = -yo_rr*(-sd*si + cd*ci) + yo_ri*(-sd*ci - cd*si) 
                           - yo_ir*( cd*si + sd*ci) + yo_ii*( cd*ci - sd*si);
                dy[[2,0]] = -yo_rr*( cd*ci - sd*si) - yo_ri*( cd*si + sd*ci) 
                           + yo_ir*(-sd*ci - cd*si) + yo_ii*(-sd*si + cd*ci);               
                dy[[3,0]] =  yo_rr*( cd*si + sd*ci) - yo_ri*( cd*ci - sd*si) 
                           - yo_ir*(-sd*si + cd*ci) + yo_ii*(-sd*ci - cd*si);

                dy[[0,1]] = (-yo_rr*sd*ci - yo_ri*sd*si + yo_ir*cd*ci + yo_ii*cd*si)*x_d;
                dy[[1,1]] = ( yo_rr*sd*si - yo_ri*sd*ci - yo_ir*cd*si + yo_ii*cd*ci)*x_d;
                dy[[2,1]] = (-yo_rr*cd*ci - yo_ri*cd*si - yo_ir*sd*ci - yo_ii*sd*si)*x_d;               
                dy[[3,1]] = ( yo_rr*cd*si - yo_ri*cd*ci + yo_ir*sd*si - yo_ii*sd*ci)*x_d;

                dy[[0,2]] = (-yo_rr*cd*si + yo_ri*cd*ci - yo_ir*sd*si + yo_ii*sd*ci)*x_i;
                dy[[1,2]] = (-yo_rr*cd*ci - yo_ri*cd*si - yo_ir*sd*ci - yo_ii*sd*si)*x_i;
                dy[[2,2]] = ( yo_rr*sd*si - yo_ri*sd*ci - yo_ir*cd*si + yo_ii*cd*ci)*x_i;               
                dy[[3,2]] = ( yo_rr*sd*ci + yo_ri*sd*si - yo_ir*cd*ci - yo_ii*cd*si)*x_i;

            });

    }

}

impl NMRFitComponent for Phase2D {

    fn get_y(&self) -> Array2<f64> {
        self.y.clone()
    }

    fn get_dydp(&self) -> Array3<f64> {
        self.dydp.clone()
    }

    //--------------------------------------
    fn eval(&mut self, p: &[f64]) {       

        (self.f)(self, p);

    }
}

#[cfg(test)]
mod tests {
    
    use ndarray::prelude::*;

    use crate::common::NMRFitComponent;
    use crate::lineshape::{Lineshape1D, Lineshape2D};

    use super::{Phase1D, Phase2D};

    #[test]
    fn phase_1d_gradient() {

        let x: Array1<f64> = Array::linspace(0.0, 1.0, 20);
        let x = x.into_shared();
        let p = vec![0.5, 0.5, 0.5, 0.5, 0.3, 0.3, 0.3, 0.3];

        // First, generate realistic peaks
        let mut lineshape = Lineshape1D::new(x.clone(), p.len());
        lineshape.eval(&p);
        let y = lineshape.y.clone();

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

        let x: Array1<f64> = Array::linspace(0.0, 1.0, 20);
        let n = x.len();
        let p = vec![0.5, 0.5, 0.5, 0.5, 0.3, 0.3, 0.3, 0.3];

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

        // First, generate realistic peaks
        let resonances: Array1<usize> = Array::zeros((p.len(),));
        let mut dimensions: Array1<usize> = Array::zeros((p.len(),));
        for i in 4 .. 8 {
            dimensions[i] = 1;
        }

        let mut lineshape = Lineshape2D::new(x1.clone(), x2.clone(), resonances, dimensions);
        lineshape.eval(&p);
        let y = lineshape.y.clone();

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


