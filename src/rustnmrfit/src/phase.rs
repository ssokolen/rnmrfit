use ndarray::{prelude::*, Zip};

//==============================================================================
// General 1D phase correction

pub struct Phase1D {

    // Initial values
    pub y_in: Array2<f64>,

    // Intermediate values
    sin_theta: Array1<f64>,
    cos_theta: Array1<f64>,
    grad_r: Array1<f64>,
    grad_i: Array1<f64>,

    // Phase function
    f: fn (&mut Phase1D, &Array<f64, Ix1>, &Array<f64, Ix1>),

    // Final output
    pub y: Array2<f64>,
    pub dydp: Array3<f64>,

}

//--------------------------------------
impl Phase1D {

    //--------------------------------------
    pub fn new(y: &Array2<f64>, n:usize, np: usize) -> Phase1D {
        
        let mut y_copy: Array2<f64> = Array::zeros((2,n));
        y_copy.assign(y);

        // Choose f based on number of parameters
        let f = match np {
            0 => Phase1D::eval_0,
            1 => Phase1D::eval_1,
            2 => Phase1D::eval_2,
            _ => panic!("Nonlinear phase correction is not supported.")
        };

        Phase1D {
            y_in: y_copy.clone(),
            
            sin_theta: Array::zeros((n,)),
            cos_theta: Array::zeros((n,)),
            grad_r: Array::zeros((n,)),
            grad_i: Array::zeros((n,)),

            f: f,

            y: y_copy,
            dydp: Array::zeros((2, np, n)),
        }
    }

    //--------------------------------------
    pub fn eval(&mut self, x: &Array<f64, Ix1>, p: &Array<f64, Ix1>) {       

        (self.f)(self, x, p);

    }

    //--------------------------------------
    pub fn eval_0(&mut self, _x: &Array<f64, Ix1>, _p: &Array<f64, Ix1>) {       

        // Nothing to do -- if there are no phasing terms, then y is unchanged

    }

    //--------------------------------------
    pub fn eval_1(&mut self, _x: &Array<f64, Ix1>, p: &Array<f64, Ix1>) {       

        // With a single phase term, generate single value intermediates
        let theta = p[0];
        let sin_theta = theta.sin();
        let cos_theta = theta.cos();

        // Real terms
        Zip::from(&mut self.y.slice_mut(s![0, ..]))
            .and(&mut self.dydp.slice_mut(s![0, 0, ..]))
            .and(&self.y_in.slice(s![0, ..]))
            .and(&self.y_in.slice(s![1, ..]))
            .apply(|y_r, dy_r, &yo_r, &yo_i| {
                *y_r =   yo_r*cos_theta + yo_i*sin_theta;
                *dy_r = -yo_r*sin_theta + yo_i*cos_theta;
            });

        // Imaginary terms
        Zip::from(&mut self.y.slice_mut(s![1, ..]))
            .and(&mut self.dydp.slice_mut(s![1, 0, ..]))
            .and(&self.y_in.slice(s![0, ..]))
            .and(&self.y_in.slice(s![1, ..]))
            .apply(|y_i, dy_i, &yo_r, &yo_i| {
                *y_i =  -yo_r*sin_theta + yo_i*cos_theta;
                *dy_i = -yo_r*cos_theta - yo_i*sin_theta;
            });
    }

    //--------------------------------------
    pub fn eval_2(&mut self, x: &Array<f64, Ix1>, p: &Array<f64, Ix1>) {       

        // With two phase terms, generate array intermediates
        let theta_0 = p[0];
        let theta_1 = p[1];

        Zip::from(x)
            .and(&mut self.sin_theta)
            .and(&mut self.cos_theta)
            .apply(|&x, s, c| {
                *s = (theta_0 + theta_1*x).sin(); 
                *c = (theta_0 + theta_1*x).cos(); 
            });

        // Real terms (1st gradient only)
        Zip::from(&mut self.y.slice_mut(s![0, ..]))
            .and(&mut self.grad_r)
            .and(&self.sin_theta)
            .and(&self.cos_theta)
            .and(&self.y_in.slice(s![0, ..]))
            .and(&self.y_in.slice(s![1, ..]))
            .apply(|y_r, dy_r, &s, &c, &yo_r, &yo_i| {
                *y_r =   yo_r*c + yo_i*s;
                *dy_r = -yo_r*s + yo_i*c;
            });

        // Imaginary terms (1st gradient only)
        Zip::from(&mut self.y.slice_mut(s![1, ..]))
            .and(&mut self.grad_i)
            .and(&self.sin_theta)
            .and(&self.cos_theta)
            .and(&self.y_in.slice(s![0, ..]))
            .and(&self.y_in.slice(s![1, ..]))
            .apply(|y_i, dy_i, &s, &c, &yo_r, &yo_i| {
                *y_i =  -yo_r*s + yo_i*c;
                *dy_i = -yo_r*c - yo_i*s;
            });

        // Slices must be filled separately to avoid borrow clashes 
        let mut dydp_r = self.dydp.slice_mut(s![0, 0, ..]);
        dydp_r.assign( &self.grad_r );
        let mut dydp_r = self.dydp.slice_mut(s![0, 1, ..]);
        dydp_r.assign( &(&self.grad_r * x) );

        let mut dydp_i = self.dydp.slice_mut(s![1, 0, ..]);
        dydp_i.assign( &self.grad_i );
        let mut dydp_i = self.dydp.slice_mut(s![1, 1, ..]);
        dydp_i.assign( &(&self.grad_i * x) );
    }

}


#[cfg(test)]
mod tests {
    
    use ndarray::prelude::*;
    use nlopt;

    use super::Phase1D;

    fn check_gradient(p: Array1<f64>) {

        // x/y must be hard coded due to nlopt
        let x = Array::from_shape_vec((2,), vec![0.3, 0.7]).unwrap();
        let y = Array::from_shape_vec((2,2), vec![0.5, 0.4, 0.3, 0.2]).unwrap();

        // Analytical gradients
        let mut phase = Phase1D::new(&y, 2, p.len());
        phase.eval(&x, &p);

        // Comparing to real numerical gradients
        fn f_real(p: &[f64]) -> f64 {
            let x = Array::from_shape_vec((2,), vec![0.3, 0.7]).unwrap();
            let y = Array::from_shape_vec((2,2), vec![0.5, 0.4, 0.3, 0.2]).unwrap();
            let p = Array::from_shape_vec((p.len(),), p.to_vec()).unwrap();
            let mut phase = Phase1D::new(&y, 2, p.len());
            phase.eval(&x, &p);
            phase.y[[0, 0]] 
        };

        let mut grad_real = vec![0.0; 2];
        nlopt::approximate_gradient(&(p.to_vec()), f_real, &mut grad_real);

        for i in 0 .. p.len() {
            assert!((phase.dydp[[0,i,0]] - grad_real[i]).abs() < 1e-5, 
                    "Phase gradient error -- real domain, term {}", i+1);
        }

        // Comparing to imaginary numerical gradients
        fn f_imag(p: &[f64]) -> f64 {
            let x = Array::from_shape_vec((2,), vec![0.3, 0.7]).unwrap();
            let y = Array::from_shape_vec((2,2), vec![0.5, 0.4, 0.3, 0.2]).unwrap();
            let p = Array::from_shape_vec((p.len(),), p.to_vec()).unwrap();
            let mut phase = Phase1D::new(&y, 2, p.len());
            phase.eval(&x, &p);
            phase.y[[1, 0]] 
        };

        let mut grad_imag = vec![0.0; 2];
        nlopt::approximate_gradient(&(p.to_vec()), f_imag, &mut grad_imag);

        for i in 0 .. p.len() {
            assert!((phase.dydp[[1,i,0]] - grad_imag[i]).abs() < 1e-5, 
                    "Phase gradient error -- imaginary domain, term {}", i+1);
        }

    }

    #[test]
    fn phase_gradient_0_order() {

        let p = Array::from_shape_vec((1,), vec![0.3]).unwrap();
        check_gradient(p);

    }

    #[test]
    fn phase_gradient_1_order() {

        let p = Array::from_shape_vec((2,), vec![0.3, 0.1]).unwrap();
        check_gradient(p);

    }
}


