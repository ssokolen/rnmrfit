use ndarray::{prelude::*, ArcArray1, stack};
use num_integer::Roots;

use crate::common::Eval;

//==============================================================================
// Basis function

//--------------------------------------
pub fn gen_basis(x: ArcArray1<f64>, knots: Array1<f64>, np: usize) -> Array2<f64> { 

    // Assuming that the knots vector includes boundary knots
    let n_knots: usize = knots.len();
    
    // The baseline degree is based on the total number of parameters and knots
    // (assuming that n_knots includes both interior and boundary knots)
    let degree: usize = np - n_knots + 1;

    let n = x.len();

    // Expanding the knots based on degree
    let left = vec![knots[0]; degree];
    let left = Array::from_shape_vec((degree,), left).unwrap();
    
    let right = vec![knots[knots.len() - 1]; degree];
    let right = Array::from_shape_vec((degree,), right).unwrap();

    let knots = stack![Axis(0), left, knots.clone(), right]; 

    fn alpha(x: f64, i: usize, j: usize, knots: &Array1<f64>) -> f64 {
        if knots[i+j] != knots[i] {
            (x - knots[i])/(knots[i+j] - knots[i])
        } else {
            0.0
        }
    };

    fn beta(x: f64, i: usize, j: usize, knots: &Array1<f64>) -> f64 {
        if j == 0 {
            if (x >= knots[i]) && (x < knots[i+1]) {
                1.0
            } else {
                0.0
            }
        } else {
            alpha(x, i, j, knots)*beta(x, i, j-1, knots) +
            ( 1.0 - alpha(x, i+1, j, knots) )*beta(x, i+1, j-1, knots)
        }
    };

    // Formulating basis
    let mut basis = Array::zeros((n, np));
    
    for i in 0 .. n {
        for j in 0 .. np {
            basis[[i,j]] = beta(x[i], j, degree, &knots);
        }
    }

    // Fudging last value by assuming boundary knots correspond to end points
    basis[[x.len() - 1, np - 1]] = 1.0;

    basis
}

//==============================================================================
// General 1D baseline

pub struct Baseline1D {

    // Number of parameters per domain
    np: usize,

    // Basis
    basis: Array2<f64>,

    // Baseline function
    f: fn (&mut Baseline1D, &[f64]),

    // Final output
    pub y: Array2<f64>,
    pub dydp: Array3<f64>,

}

//--------------------------------------
impl Baseline1D {

    //--------------------------------------
    pub fn new(x: ArcArray1<f64>, knots: Array1<f64>, np: usize) -> Baseline1D {

        let n: usize = x.len();

        // For baseline terms np represents only one of real/imaginary domains
        // Meaning that it has to be multiplied for the dydp term

        let mut basis: Array2<f64> = Array::zeros((n, np));
        let mut dydp: Array3<f64> = Array::zeros((2, np*2, n));

        if np > 0 {
            let basis_0: Array2<f64> = Array::zeros((n, np));
            basis = gen_basis(x, knots, np);

            let mut dydp_r = dydp.slice_mut(s![0, .., ..]);
            dydp_r.assign(&( stack![Axis(0), basis.t(), basis_0.t()] ));

            let mut dydp_i = dydp.slice_mut(s![1, .., ..]);
            dydp_i.assign(&( stack![Axis(0), basis_0.t(), basis.t()] ));
        }

        let f = match np {
            0 => Baseline1D::eval_0,
            _ => Baseline1D::eval_n
        };

        Baseline1D {
            np: np,
            basis: basis,

            f: f,

            y: Array::zeros((2, n)),
            dydp: dydp,
        }
    }

    //--------------------------------------
    pub fn eval_0(&mut self, _p: &[f64]) {       

        // Nothing to do -- if there are no baseline terms, then y is unchanged

    }

    //--------------------------------------
    pub fn eval_n(&mut self, p: &[f64]) {       

        // Assume that p is divided in half -- real parameters then imaginary ones,
        // so loop over the two halves
        let np = self.np;
        for i in 0 .. 2 {
            let p_slice = Array::from_shape_vec((np,), p[(np*i) .. (np*(i+1))].to_vec()).unwrap();
            let mut y_slice = self.y.slice_mut(s![i, ..]);
            y_slice.assign( &(self.basis.dot(&p_slice)) );
        }

    }
}

impl Eval for Baseline1D {

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
// General 2D baseline

pub struct Baseline2D {

    // Number of parameters per domain
    np: usize,

    // Basis
    basis: Array2<f64>,

    // Baseline function
    f: fn (&mut Baseline2D, &[f64]),

    // Final output
    pub y: Array2<f64>,
    pub dydp: Array3<f64>,

}

//--------------------------------------
impl Baseline2D {

    //--------------------------------------
    pub fn new(x_direct: ArcArray1<f64>, x_indirect: ArcArray1<f64>, 
               knots: Array1<f64>, np: usize) -> Baseline2D {

        let n: usize = x_direct.len();
        
        // The total number of parameter must be split between direct and indirect grids
        // (assume square grid so take square root)
        let np0 = np.sqrt();
        assert!(np0*np0 == np, "The number of 2D baseline parameters must be square.");

        // For baseline terms np represents only one of real/imaginary domains
        // Meaning that it has to be multiplied for the dydp term

        let mut basis: Array2<f64> = Array::zeros((n, np));
        let mut dydp: Array3<f64> = Array::zeros((4, np*4, n));

        if np > 0 {
            let basis_0: Array2<f64> = Array::zeros((n, np));
            
            // Although it would probably be more efficient to construct the bases
            // from unique x values... that optimization can come later
            let basis_direct = gen_basis(x_direct, knots.clone(), np0);
            let basis_indirect = gen_basis(x_indirect, knots.clone(), np0);

            let mut temp_col = Array::zeros((n,));

            // Filling 2D matrix as product of 1D ones
            for i in 0 .. np0 {
                let direct_slice = basis_direct.slice(s![.., i]);
                for j in 0 .. np0 {
                    let mut basis_slice = basis.slice_mut(s![.., i*np0 + j]);
                    let indirect_slice = basis_indirect.slice(s![.., j]);

                    temp_col.assign(&direct_slice);
                    temp_col *= &indirect_slice;

                    basis_slice.assign(&temp_col);
                }
            }

            let mut dydp_rr = dydp.slice_mut(s![0, .., ..]);
            dydp_rr.assign(&( stack![Axis(0), basis.t(), basis_0.t(), basis_0.t(), basis_0.t()] ));

            let mut dydp_ri = dydp.slice_mut(s![1, .., ..]);
            dydp_ri.assign(&( stack![Axis(0), basis_0.t(), basis.t(), basis_0.t(), basis_0.t()] ));

            let mut dydp_ir = dydp.slice_mut(s![2, .., ..]);
            dydp_ir.assign(&( stack![Axis(0), basis_0.t(), basis_0.t(), basis.t(), basis_0.t()] ));

            let mut dydp_ii = dydp.slice_mut(s![3, .., ..]);
            dydp_ii.assign(&( stack![Axis(0), basis_0.t(), basis_0.t(), basis_0.t(), basis.t()] ));
        }

        let f = match np {
            0 => Baseline2D::eval_0,
            _ => Baseline2D::eval_n
        };

        Baseline2D {
            np: np,
            basis: basis,

            f: f,

            y: Array::zeros((4, n)),
            dydp: dydp,
        }
    }

    //--------------------------------------
    pub fn eval_0(&mut self, _p: &[f64]) {       

        // Nothing to do -- if there are no baseline terms, then y is unchanged

    }

    //--------------------------------------
    pub fn eval_n(&mut self, p: &[f64]) {       

        // Assume that p is divided in four -- rr, ri, ir, ii
        // so loop over the quarters
        let np = self.np;
        for i in 0 .. 4 {
            let p_slice = Array::from_shape_vec((np,), p[(np*i) .. (np*(i+1))].to_vec()).unwrap();
            let mut y_slice = self.y.slice_mut(s![i, ..]);
            y_slice.assign( &(self.basis.dot(&p_slice)) );
        }

    }
}

impl Eval for Baseline2D {

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
// Unit tests

#[cfg(test)]
mod tests {

    use ndarray::prelude::*;

    use crate::common::Eval;

    //--------------------------------------
    #[test]
    fn baseline_1d_gradient() {

        let x: Array1<f64> = Array::linspace(0.0, 1.0, 20);
        let x = x.into_shared();
        let p = Array::linspace(0.0, 1.0, 16).to_vec();
        let knots = Array::linspace(0.0, 1.0, 3);

        let mut baseline = super::Baseline1D::new(x, knots, p.len()/2);

        // Summing up errors across one dimensions increases threshold
        baseline.grad_test(&p, 5e-4);

    }

    //--------------------------------------
    #[test]
    fn baseline_2d_gradient() {

        let x: Array1<f64> = Array::linspace(0.0, 1.0, 20);
        let n = x.len();

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
        let p = Array::linspace(0.0, 1.0, 64).to_vec();
        let knots = Array::linspace(0.0, 1.0, 3);

        let mut baseline = super::Baseline2D::new(x1, x2, knots, p.len()/4);

        // Summing up errors across two dimensions increases threshold
        baseline.grad_test(&p, 5e-3);

    }

}
