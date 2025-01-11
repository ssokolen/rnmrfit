use itertools::izip;
use ndarray::{prelude::*, ArcArray1, concatenate};
use ndarray_stats::QuantileExt;
use num_integer::Roots;

use crate::fitting::common::ComponentEval;


//=============================================================================
// Commonly used types

type YSlice<'a> = ArrayViewMut<'a, f64, Ix2>;
type DySlice<'a> = ArrayViewMut<'a, f64, Ix3>;



//=============================================================================
// General 1D baseline

pub struct Baseline1D {

    // Vector of spline baseis
    basis: Array2<f64>,

    // Work arrays
    p_r: Array1<f64>,
    p_i: Array1<f64>,

    // Derivatives are constant and can be generated in advance
    dy: Array3<f64>
}

//--------------------------------------
impl Baseline1D {

    //--------------------------------------
    pub fn new(basis: Array2<f64>) -> Baseline1D {

        // Since derivatives are constant, they can be generated in advance.
        let n = basis.nrows();
        let np = basis.ncols();

        let mut dy: Array3<f64> = Array::zeros((np*2, n, 2));

        for j in 0..np {

            let slice = basis.slice(s![.., j]);

            // Real components are only influenced by even terms,
            // (and imaginary components by odd terms)
            let (mut dy_r, mut dy_i) = dy.multi_slice_mut((
                s![j*2, .., 0],
                s![j*2 + 1, .., 1]
            ));

            let iterator = izip!(
                dy_r.iter_mut(),
                dy_i.iter_mut(), 
                slice.iter()
            );

            for (r, i, value) in iterator {
                *r = *value;
                *i = *value;
            }
        }

        let mut p: Array1<f64> = Array::zeros((np));

        Baseline1D {
            basis: basis,
            p_r: p.clone(),
            p_i: p,
            dy: dy,
        }
    }
}


impl ComponentEval for Baseline1D {

    fn nd(&self) -> usize {
       2 
    }

    fn n(&self) -> usize {
        let (_, n, _) = self.dy.dim();
        n
    }

    //--------------------------------------
    fn eval(&mut self, p: &[f64], mut y: YSlice, mut dy: DySlice) { 

        // Transfer parameters to work arrays
        let iterator = self.p_r.iter_mut()
            .zip(p.iter().step_by(2));

        for (r, p) in iterator {
            *r = *p;
        }

        let iterator = self.p_i.iter_mut()
            .zip(p.iter().skip(1).step_by(2));

        for (i, p) in iterator {
            *i = *p;
        }

        // First the real
        let y_temp = self.basis.dot(&self.p_r);
        let mut y_slice = y.slice_mut(s![.., 0]);
        y_slice.assign(&y_temp);

        // Then the imaginary
        let y_temp = self.basis.dot(&self.p_i);
        let mut y_slice = y.slice_mut(s![.., 1]);
        y_slice.assign(&y_temp);

        dy.assign(&self.dy);
    }
}

/*
//=============================================================================
// General 2D baseline
#[allow(dead_code)]
pub struct Baseline2D {

    // Basis (effectively capturing x)
    basis: Array2<f64>,

    // Baseline function based on parameter number
    np: usize,
    f: fn (&mut Baseline2D, &[f64], YSlice, DySlice),

    // Derivatives are constant and can be generated in advance
    dy: Array3<f64>,
}

//--------------------------------------
#[allow(dead_code)]
impl Baseline2D {

    //--------------------------------------
    pub fn new(
            x_direct: ArcArray1<f64>, 
            x_indirect: ArcArray1<f64>, 
            knots: Array1<f64>, 
            np: usize
        ) -> Baseline2D {

        let n: usize = x_direct.len();
        
        // The total number of parameter must be split between direct and 
        // indirect grids (assume square grid so take square root)
        let np0 = np.sqrt();
        assert!(
            np0*np0 == np,
            "The number of 2D baseline parameters must be square"
        );

        let mut basis: Array2<f64> = Array::zeros((n, np));
        let mut dy: Array3<f64> = Array::zeros((np*4, n, 4));

        // Basis function only has to be built for non-zero np
        if np > 0 {
            let basis_0: Array2<f64> = Array::zeros((n, np));
            
            // TODO: consider constructing basis from unique values 
            let basis_direct = gen_basis(x_direct, knots.clone(), np0);
            let basis_indirect = gen_basis(x_indirect, knots.clone(), np0);

            // Filling 2D matrix as product of 1D ones
            for i in 0 .. np0 {
                let direct_slice = basis_direct.slice(s![.., i]);
                for j in 0 .. np0 {
                    let mut basis_slice = basis.slice_mut(s![.., i + j*np0]);
                    let indirect_slice = basis_indirect.slice(s![.., j]);

                    basis_slice.assign(&direct_slice);
                    basis_slice *= &indirect_slice;
                }
            }

            let mut dy_rr = dy.slice_mut(s![.., .., 0]);
            dy_rr.assign(&( concatenate![
                Axis(0), basis.t(), basis_0.t(), basis_0.t(), basis_0.t()
            ] ));

            let mut dy_ri = dy.slice_mut(s![.., .., 1]);
            dy_ri.assign(&( concatenate![
                Axis(0), basis_0.t(), basis.t(), basis_0.t(), basis_0.t()
            ] ));

            let mut dy_ir = dy.slice_mut(s![.., .., 2]);
            dy_ir.assign(&( concatenate![
                Axis(0), basis_0.t(), basis_0.t(), basis.t(), basis_0.t()
            ] ));

            let mut dy_ii = dy.slice_mut(s![.., .., 3]);
            dy_ii.assign(&( concatenate![
                Axis(0), basis_0.t(), basis_0.t(), basis_0.t(), basis.t()
            ] ));
        }

        let f = match np {
            0 => Baseline2D::eval_0,
            _ => Baseline2D::eval_n
        };

        Baseline2D {
            basis: basis,
            np: np,
            f: f,
            dy: dy,
        }
    }

    //--------------------------------------
    pub fn eval_0(&mut self, _p: &[f64], _y: YSlice, _dy: DySlice) {       

        // Nothing to do -- if there are no baseline terms, then y is unchanged

    }

    //--------------------------------------
    pub fn eval_n(&mut self, p: &[f64], mut y: YSlice, mut dy: DySlice) {

        // Assume that p is divided in four -- rr, ri, ir, ii
        // so loop over the quarters
        let np = self.np;
        for i in 0 .. 4 {
            let p_slice = Array::from_shape_vec(
                (np,), p[(np*i) .. (np*(i+1))].to_vec()
            ).unwrap();
            let mut y_slice = y.slice_mut(s![.., i]);
            y_slice.assign( &(self.basis.dot(&p_slice)) );
        }

        dy.assign( &(self.dy) );
    }
}


impl ComponentEval for Baseline2D {

    fn nd(&self) -> usize {
        4
    }

    fn n(&self) -> usize {
        let (n, _) = self.basis.dim();
        n
    }


    //--------------------------------------
    fn eval(&mut self, p: &[f64], y: YSlice, dy: DySlice) {       

        (self.f)(self, p, y, dy);
    }
}
*/

//=============================================================================
// Unit tests

#[cfg(test)]
mod tests {

    use ndarray::prelude::*;

    use crate::fitting::common::ComponentEval;

    /*
    //--------------------------------------
    #[test]
    fn baseline_1d_gradient() {

        let n = 20;
        let x: Array1<f64> = Array::linspace(0.0, 1.0, n);
        let x = x.into_shared();

        let p = Array::linspace(0.0, 1.0, 16).to_vec();
        let knots = Array::linspace(0.0, 1.0, 3);

        let mut baseline = super::Baseline1D::new(x, knots, p.len()/2);

        // Summing up errors increases threshold
        baseline.grad_test(&p, 5e-4);
    }

    //--------------------------------------
    #[test]
    fn baseline_2d_gradient() {

        let n = 20;
        let x: Array1<f64> = Array::linspace(0.0, 1.0, n);

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

        // Summing up errors increases threshold
        baseline.grad_test(&p, 5e-4);
    }
    */
}
