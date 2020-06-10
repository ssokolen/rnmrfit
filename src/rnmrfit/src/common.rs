use ndarray::prelude::*;

//------------------------------------------------------------------------------
// Single eval trait to define a generic testing function

pub trait NMRFitComponent {
    fn get_y(&self) -> Array2<f64>;
    fn get_dydp(&self) -> Array3<f64>;

    fn eval(&mut self, p: &[f64]);

    fn grad_test(&mut self, p: &[f64], tol: f64) {

        let mut p = p.to_vec();

        // Calculate analytical gradients
        self.eval(&p);
        let mut analytical = self.get_dydp();
        
        // Calculate numerical gradients
        // (numerical gradient estimate based on rust NLOPT approximate_gradient)
        let mut numerical = Array::zeros(analytical.raw_dim());

        let eps = std::f64::EPSILON.powf(1.0 / 3.0);

        for i in 0 .. p.len() {

            let mut slice: ArrayViewMut2<f64> = numerical.slice_mut(s![.., i, ..]);
        
            let p0 = p[i];

            p[i] = p0 + eps;
            self.eval(&p);
            slice.assign(&self.get_y());
            
            p[i] = p0 - eps;
            self.eval(&p);
            slice -= &self.get_y();

            slice /= 2.0 * eps;
            p[i] = p0;

        }

        // Taking absolute difference and summing up y values
        analytical -= &numerical;
        analytical.map_mut(|x| *x = (*x).abs() );

        let differences = analytical.sum_axis(Axis(2));

        for i in 0 .. differences.len_of(Axis(0)) {
            for j in 0 .. differences.len_of(Axis(1)) {
                assert!(differences[[i,j]] < tol, 
                        "Problem with parameter {} in dimension {}, difference: {}",
                        j, i, differences[[i,j]])
            }
        }
    }
}

//------------------------------------------------------------------------------
// Separate trait for testing least squares fit

pub trait NMRFit {
    fn eval(&mut self, p: &[f64], grad: Option<&mut [f64]>) -> f64;

    fn grad_test(&mut self, p: &[f64], tol: f64) {

        let mut p = p.to_vec();

        // Calculate analytical gradients
        let mut analytical = p.clone();
        self.eval(&p, Some(&mut analytical));
        
        // Calculate numerical gradients
        // (numerical gradient estimate based on rust NLOPT approximate_gradient)
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
