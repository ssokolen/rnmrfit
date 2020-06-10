use core::ops::Range;
use ndarray::{prelude::*, ArcArray1, Zip};
use num::complex::Complex;
use std::collections::HashMap;

use crate::peak::{Peak, PeakFunctions};
use crate::common::NMRFitComponent;

//==============================================================================
// General 1D lineshape

pub struct Lineshape1D {

    // Reference to x values for convenience
    x: ArcArray1<f64>,

    // Final output
    pub y: Array2<f64>,
    pub dydp: Array3<f64>,

}

//--------------------------------------
impl Lineshape1D {

    //--------------------------------------
    pub fn new(x: ArcArray1<f64>, np: usize) -> Lineshape1D {

        let n = x.len();

        Lineshape1D {
            x: x,
            
            y: Array::zeros((2, n)),
            dydp: Array::zeros((2, np, n)),
        }
    }

    //------------------------------------------------------------------------------
    fn eval_peak(&mut self, mut peak: impl PeakFunctions, i: usize) {

        // Temporary holder for gradient terms
        let mut grad: Vec<Complex<f64>> = vec![Complex::new(0.0, 0.0); 4];

        let mut y = self.y.slice_mut(s![.., ..]);
        let mut dydp = self.dydp.slice_mut(s![.., i .. (i+4), ..]);

        // Looping over x
        for j in 0 .. self.x.len() {

            let fit = peak.gradients(self.x[j], &mut grad);

            y[[0, j]] += fit.re;
            y[[1, j]] += fit.im;

            for k in 0 .. 4 {

                dydp[[0, k, j]] = grad[k].re;
                dydp[[1, k, j]] = grad[k].im;

            }
        }
    }
}

//--------------------------------------
impl NMRFitComponent for Lineshape1D {

    fn get_y(&self) -> Array2<f64> {
        self.y.clone()
    }

    fn get_dydp(&self) -> Array3<f64> {
        self.dydp.clone()
    }

    //--------------------------------------
    fn eval(&mut self, p: &[f64]) {

        // Initializing y to zero (unlike dydp, y is incremented)
        self.y.fill(0.0);

         // Loop through each peak
        for i in (0 .. p.len()).step_by(4) {

            // Generate peak object and match based on result
            let peak = Peak::new(p[i], p[i+1], p[i+2], p[i+3]);

            match peak {
                Peak::Lorentz( lorentz ) => self.eval_peak(lorentz, i),
                Peak::Voigt( voigt ) => self.eval_peak(voigt, i),
            };
        }
    }   
}

//==============================================================================
// Container for 1D lineshape that maps repeated values

pub struct Lineshape1DMap {

    // Reference to x indices for convenience
    x_map: ArcArray1<usize>,

    // Lineshape object for unique values
    lineshape: Lineshape1D,

    // Final output
    pub y: Array2<f64>,
    pub dydp: Array3<f64>,

}

//--------------------------------------
impl Lineshape1DMap {

    //--------------------------------------
    pub fn new(x: ArcArray1<f64>, x_map: ArcArray1<usize>, np: usize) -> Lineshape1DMap {

        let n = x_map.len();

        Lineshape1DMap {
            x_map: x_map,

            lineshape: Lineshape1D::new(x, np),
            
            y: Array::zeros((2, n)),
            dydp: Array::zeros((2, np, n)),
        }
    }

    //--------------------------------------
    pub fn eval(&mut self, p: &[f64]) {       

        // First, evaluate unique values
        self.lineshape.eval(p);

        // And then map unique values to repeats
        let np = p.len();

        for i in 0 .. 2 {

            // First, the y values themselves
            let from: ArrayView<_, Ix1> = self.lineshape.y.slice(s![i, ..]);   
            let to: ArrayViewMut<_, Ix1> = self.y.slice_mut(s![i, ..]);   

            Zip::from(to).and(&self.x_map).apply(|x, &i| *x = from[i]);

            // Then, the derivatives
            for j in 0 .. np {
                let from: ArrayView<_, Ix1> = self.lineshape.dydp.slice(s![i, j, ..]);   
                let to: ArrayViewMut<_, Ix1> = self.dydp.slice_mut(s![i, j, ..]);   

                Zip::from(to).and(&self.x_map).apply(|x, &i| *x = from[i]);

            }

        }
    }
}


//==============================================================================
// General 2D lineshape

pub struct Lineshape2D {

    // Pairs of 1D lineshapes and ranges mapping parameters to specific resonances
    lineshapes: Vec<Vec<Lineshape1DMap>>,
    ranges: Vec<Vec<Range<usize>>>,

    // Intermediate value
    y_temp:Array1<f64>,

    // Final output
    pub y: Array2<f64>,
    pub dydp: Array3<f64>,

}

//--------------------------------------
impl Lineshape2D {

    //--------------------------------------
    pub fn new(x_direct: ArcArray1<f64>, x_indirect: ArcArray1<f64>, 
               resonances: Array1<usize>, dimensions: Array1<usize>) 
        -> Lineshape2D {

        let n = x_direct.len();

        // Breaking down x values into a subset of unique values and maps
        let (x_direct, x_direct_map) = Lineshape2D::map_unique(x_direct);
        let (x_indirect, x_indirect_map) = Lineshape2D::map_unique(x_indirect);

        // Generating individual ranges for overall parameter list
        let np = resonances.len();

        let ranges = Lineshape2D::gen_ranges(resonances, dimensions);

        // Generating a vector of lineshapes from the ranges
        let mut lineshapes: Vec<Vec<Lineshape1DMap>> = Vec::new();

        for i in 0 .. ranges.len() {
            lineshapes.push(vec![
                Lineshape1DMap::new(x_direct.clone(), x_direct_map.clone(), ranges[i][0].len()),
                Lineshape1DMap::new(x_indirect.clone(), x_indirect_map.clone(), ranges[i][1].len())
            ]);
        }

        Lineshape2D {
            lineshapes: lineshapes,
            ranges: ranges,

            y_temp: Array::zeros((n,)),
            
            y: Array::zeros((4, n)),
            dydp: Array::zeros((4, np, n)),
        }
    }


    //--------------------------------------
    fn map_unique(x: ArcArray1<f64>) -> (ArcArray1<f64>, ArcArray1<usize>) {
        
        let mut hash_map = HashMap::new();
        let mut index_map: Vec<usize> = Vec::new();
        let mut x_unique: Vec<f64> = Vec::new();
        
        for i in 0 .. x.len() {
        
            // x-values are assumed to be normalized 0-1
            // and imprecision of 1e-10 should be acceptable
            let x_key = (x[i]*1e10) as i64;
            
            match hash_map.get(&x_key) {
                Some(j) => {
                    index_map.push(*j);  
                },
                None => {
                    let j = x_unique.len();
                    hash_map.insert(x_key, j);
                    x_unique.push(x[i]);
                    index_map.push(j);
                }
            }
            
        }
        
        (Array::from_shape_vec((x_unique.len(),), x_unique).unwrap().into_shared(),
         Array::from_shape_vec((index_map.len(),), index_map).unwrap().into_shared())
    }

    //--------------------------------------
    pub fn gen_ranges(resonances: Array1<usize>, dimensions: Array1<usize>) 
        -> Vec<Vec<Range<usize>>> {
        
        let mut indexes: Vec<Vec<Vec<usize>>> = Vec::new();

        let mut ir: usize;
        let mut id: usize;

        // First, splitting up all the indexes
        for i in 0 .. resonances.len() {
        
            ir = resonances[i];
            id = dimensions[i];

            // Ensure that there are sufficient resonances
            while indexes.len() <= ir {
                indexes.push(vec![Vec::new(), Vec::new()]);
            }

            // Adding indexes
            indexes[ir][id].push(i);
        }

        // Then converting them to ranges if they are contiguous
        let mut ranges: Vec<Vec<Range<usize>>> = Vec::new();

        let mut lengths: Vec<usize> = vec![0, 0];

        for i in 0 .. indexes.len() {

            // Checking if contiguous
            for j in 0 .. 2 {

                lengths[j] = indexes[i][j].len();
                for k in 0 .. (lengths[j] - 1) {
                    if (indexes[i][j][k+1] - indexes[i][j][k]) > 1 {
                        panic!("Resonances within the parameter vector must be contiguous.");
                    }
                }
            }

            // If no issues, then generate ranges
            ranges.push(vec![indexes[i][0][0] .. (indexes[i][0][lengths[0]-1] + 1),
                             indexes[i][1][0] .. (indexes[i][1][lengths[1]-1] + 1)]);
        }
        ranges
    }
}

impl NMRFitComponent for Lineshape2D {

    fn get_y(&self) -> Array2<f64> {
        self.y.clone()
    }

    fn get_dydp(&self) -> Array3<f64> {
        self.dydp.clone()
    }

    //--------------------------------------
    fn eval(&mut self, p: &[f64]) {       

        self.y.fill(0.0);

        // Loop over the resonances
        for i in 0 .. self.lineshapes.len() {

            // First, perform individual evaluations
            let r_direct = self.ranges[i][0].clone();
            let r_indirect = self.ranges[i][1].clone();

            self.lineshapes[i][0].eval(&p[r_direct.clone()]);
            self.lineshapes[i][1].eval(&p[r_indirect.clone()]);

            // Then perform cross product
            let mut j = 0;
            for j1 in 0 .. 2 {
                for j2 in 0 .. 2 {

                    let y_direct = self.lineshapes[i][0].y.slice(s![j1, ..]);
                    let y_indirect = self.lineshapes[i][1].y.slice(s![j2, ..]);

                    let dydp_direct = self.lineshapes[i][0].dydp.slice(s![j1, .., ..]);
                    let dydp_indirect = self.lineshapes[i][1].dydp.slice(s![j2, .., ..]);

                    // Trying to minimize temporary array creation
                    self.y_temp.assign( &y_direct );
                    self.y_temp *= &y_indirect;

                    let mut y = self.y.slice_mut(s![j, ..]);
                    y += &self.y_temp;

                    // Gradients correspond to one set of components -- direct or indirect

                    // First, the direct gradients
                    let mut dydp = self.dydp.slice_mut(s![j, r_direct.clone(), ..]);
                    dydp.assign( &( &dydp_direct * &y_indirect ) );

                    // And the indirect
                    let mut dydp = self.dydp.slice_mut(s![j, r_indirect.clone(), ..]);
                    dydp.assign( &( &dydp_indirect * &y_direct ) );

                    j += 1;
                }
            }
        }
    }
}

//==============================================================================
// Unit tests

#[cfg(test)]
mod tests {

    use ndarray::prelude::*;

    use crate::common::NMRFitComponent;

    //--------------------------------------
    #[test]
    fn map() {

        let x  = Array::linspace(0.0, 1.0, 5).into_shared();

        let xi = vec![0, 0, 1, 1, 2, 2, 3, 3, 4, 4];
        let xi = Array::from_shape_vec((xi.len(),), xi).unwrap().into_shared();

        let p = vec![0.5, 0.5, 0.5, 0.5, 0.3, 0.3, 0.3, 0.3];

        let mut l = super::Lineshape1D::new(x.clone(), p.len());
        let mut lmap = super::Lineshape1DMap::new(x.clone(), xi.clone(), p.len());

        l.eval(&p);
        lmap.eval(&p);

        for i in 0 .. x.len() {

            for j in 0 .. 2 {

                assert!(l.y[[0, i]] == lmap.y[[0, i*2+j]], 
                        "map y mismatch");
                assert!(l.y[[1, i]] == lmap.y[[1, i*2+j]], 
                        "map y mismatch");

                for k in 0 .. p.len() {

                    assert!(l.dydp[[0, k, i]] == lmap.dydp[[0, k, i*2+j]], 
                            "map dydp mismatch");
                    assert!(l.dydp[[1, k, i]] == lmap.dydp[[1, k, i*2+j]],
                            "map dydp mismatch");

                }
            }
        }
    }

    //--------------------------------------
    #[test]
    fn lineshape_1d_gradient() {

        let x: Array1<f64> = Array::linspace(0.0, 1.0, 20);
        let x = x.into_shared();
        let p = vec![0.5, 0.5, 0.5, 0.5, 0.3, 0.3, 0.3, 0.3];

        let mut lineshape = super::Lineshape1D::new(x, p.len());

        // Summing up errors across one dimensions increases threshold
        lineshape.grad_test(&p, 5e-4);

    }

    //--------------------------------------
    #[test]
    fn lineshape_2d_gradient() {

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

        let p = vec![0.5, 0.5, 0.5, 0.5, 0.3, 0.3, 0.3, 0.3];
        let resonances: Array1<usize> = Array::zeros((p.len(),));
        let mut dimensions: Array1<usize> = Array::zeros((p.len(),));
        for i in 4 .. 8 {
            dimensions[i] = 1;
        }

        let mut lineshape = super::Lineshape2D::new(x1, x2, resonances, dimensions);

        // Summing up errors across two dimensions increases threshold
        lineshape.grad_test(&p, 5e-3);

    }

}
