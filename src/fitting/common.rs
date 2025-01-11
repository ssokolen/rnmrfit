use ndarray::prelude::*;

//=============================================================================
// Commonly used types

type YSlice<'a> = ArrayViewMut<'a, f64, Ix2>;
type DySlice<'a> = ArrayViewMut<'a, f64, Ix3>;


//=============================================================================
// Single eval trait to define a generic testing function

pub trait ComponentEval {

    // Number of complex dimensions
    fn nd(&self) -> usize;

    // Number of data points in x
    fn n(&self) -> usize; 

    // Lineshape calculation and its derivative given a set of parameters
    fn eval(&mut self, p: &[f64], y: YSlice, dy: DySlice);

    // Gradient test
    fn grad_test(&mut self, p: &[f64], tol: f64) {

        let n = self.n();
        let nd = self.nd();
        let np = p.len();

        let mut p_mod = p.to_vec();

        // Analytical numerical gradients and placeholders
        let mut y: Array2<f64> = Array::zeros((n, nd));
        let mut dy: Array3<f64> = Array::zeros((np, n, nd));

        // Calculate numerical gradients
        // (numerical gradient estimate based on NLOPT approximate_gradient)
        let mut dynm: Array3<f64> = Array::zeros((np, n, nd));

        let eps = std::f64::EPSILON.powf(1.0 / 3.0);

        for (i, p0) in p_mod.clone().into_iter().enumerate() {

            let mut slice = dynm.slice_mut(s![i, .., ..]);
        
            p_mod[i] = p0 + eps;

            let y_slice = y.slice_mut(s![.., ..]);
            let dy_slice = dy.slice_mut(s![.., .., ..]);
            self.eval(&p_mod[..], y_slice, dy_slice);
            slice.assign(&y);
            
            p_mod[i] = p0 - eps;

            let y_slice = y.slice_mut(s![.., ..]);
            let dy_slice = dy.slice_mut(s![.., .., ..]);
            self.eval(&p_mod[..], y_slice, dy_slice);
            slice -= &y;

            slice /= 2.0 * eps;
            p_mod[i] = p0;
        }

        // Calculate analytical gradients
        let y_slice = y.slice_mut(s![.., ..]);
        let dy_slice = dy.slice_mut(s![.., .., ..]);
        self.eval(&p, y_slice, dy_slice);

        // Taking absolute difference and summing up y values
        dy -= &dynm;
        dy.map_mut(|x| *x = (*x).abs() );

        let differences = dy.sum_axis(Axis(1));

        println!("{:?}", differences);

        for i in 0 .. differences.len_of(Axis(0)) {
            for j in 0 .. differences.len_of(Axis(1)) {
                assert!(differences[[i,j]] < tol, 
                    "Problem with parameter {} in dimension {} ({})",
                        i, j, differences[[i,j]])
            }
        }
    }
}

//-----------------------------------------------------------------------------
// Separate trait for testing least squares fit

pub trait FitEval {
    fn eval(&mut self, p: &[f64], grad: Option<&mut [f64]>) -> f64;

    fn grad_test(&mut self, p: &[f64], tol: f64) {

        let mut p = p.to_vec();

        // Calculate analytical gradients
        let mut analytical = p.clone();
        self.eval(&p, Some(&mut analytical));
        
        // Calculate numerical gradients
        // (numerical gradient estimate based on NLOPT approximate_gradient)
        let mut numerical = p.clone();
        let mut differences = p.clone();
        
        let eps = std::f64::EPSILON.powf(1.0 / 3.0);

        for i in 0 .. p.len() {

            let p0 = p[i];

            p[i] = p0 + eps;
            numerical[i] = self.eval(&p, None);
            
            p[i] = p0 - eps;
            numerical[i] -= self.eval(&p, None);
            
            numerical[i] /= 2.0 * eps;
            p[i] = p0;

            // Converting analytical into absolute difference
            differences[i] = (analytical[i] - numerical[i]).abs();
        }

        for i in 0 .. differences.len() {
            assert!(differences[i] < tol, 
                    "Problem with parameter {}, difference: {}",
                    i, differences[i])
        }
    }
}
