use core::ops::Range;
use itertools::{
    izip, 
    Itertools
};
use ndarray::{prelude::*, ArcArray1};
use std::{
    clone::Clone,
    collections::HashMap,
};

use crate::fitting::{
    peak::f_peak_gradients,
    common::ComponentEval,
};


//=============================================================================
// Commonly used types

type YSlice<'a> = ArrayViewMut<'a, f64, Ix2>;
type DySlice<'a> = ArrayViewMut<'a, f64, Ix3>;


//=============================================================================
// General 1D lineshape

pub struct Lineshape1D {

    // Reference to x values for convenience
    x: ArcArray1<f64>,

    // Work array
    y: Array2<f64>,

    // Groups linking peaks to ranges
    ranges: Vec<Range<usize>>,
    groups: Vec<usize>,

    // Error tolerance for Fadeeva function
    tol: f64
}


//--------------------------------------
impl Lineshape1D {

    //--------------------------------------
    pub fn new(
            x: ArcArray1<f64>, 
            ranges: Vec<Range<usize>>, 
            groups: Vec<usize>, 
            tol: f64
        ) -> Lineshape1D {

        let n = x.len();
        
        Lineshape1D {
            x: x,
            y: Array::zeros((n, 2)),
            ranges: ranges,
            groups: groups,
            tol: tol,
        }
    }
}


//--------------------------------------
impl ComponentEval for Lineshape1D {

    //--------------------------------------
    fn nd(&self) -> usize {
        2
    }

    //--------------------------------------
    fn n(&self) -> usize {
        self.x.len()
    }

    //--------------------------------------
    fn eval(&mut self, p: &[f64], mut y: YSlice, mut dy: DySlice) {

        // Initialize y to zero
        y.fill(0.0);

        let tol = self.tol;

        // Loop through each peak
        for (i, j) in (0 .. p.len()).step_by(4).enumerate() {

            // Slicing work arrays
            let group = &self.groups[i];
            let range = &self.ranges[*group];
            
            let x = self.x.slice(s![range.clone()]);
            let y_temp = self.y.slice_mut(s![range.clone(), ..]);

            let (gf, gh, gw, gp) = dy.multi_slice_mut((
                s![j + 3, range.clone(), ..],
                s![j + 2, range.clone(), ..],
                s![j + 1, range.clone(), ..],
                s![j, range.clone(), ..],
            ));

            // Generate lineshape
            f_peak_gradients(&p[j..(j+4)], x, y_temp, (gp, gw, gh, gf), tol);

            // And then update cumulative y
            let mut y_slice = y.slice_mut(s![range.clone(), ..]);
            let y_temp = self.y.slice_mut(s![range.clone(), ..]);

            y_slice += &y_temp;
        }
    }   
}

/*
//=============================================================================
// Container for 1D lineshape that maps repeated values
#[allow(dead_code)]
pub struct Lineshape1DMap {

    // Reference to x indices for convenience
    x_map: ArcArray1<usize>,

    // Lineshape object for unique values
    lineshape: Lineshape1D,

    // Work arrays passed to Lineshape1D
    y: Array2<f64>,
    dydp: Array3<f64>,
}


//--------------------------------------
#[allow(dead_code)]
impl Lineshape1DMap {


    //--------------------------------------
    fn nd(&self) -> usize {
        2
    }

    //--------------------------------------
    fn n(&self) -> usize {
        self.x_map.len()
    }

    //--------------------------------------
    // Generate a map of unique values
    fn map(x: ArcArray1<f64>) -> (ArcArray1<f64>, ArcArray1<usize>) {

        let mut map = HashMap::new();
        let mut unique: Vec<f64> = Vec::new();
        let mut index: Vec<usize> = Vec::new();
        let mut i = 0;

        for x in x.into_iter() {
        
            // x-values are assumed to be normalized 0-1
            // and imprecision of 1e-8 is acceptable
            let key = (x * 1e8) as i64;

            match map.get(&key) {
                Some(j) => {
                    index.push(*j);  
                },
                None => {
                    map.insert(key, i);
                    unique.push(x);
                    index.push(i);  
                    i += 1;
                }
            }
        }

        let unique = Array::from_shape_vec((unique.len(),), unique)
            .unwrap()
            .into_shared();

        let index = Array::from_shape_vec((index.len(),), index)
            .unwrap()
            .into_shared();
        
        ( unique, index )
    }


    //--------------------------------------
    pub fn new(x: ArcArray1<f64>, x_map: ArcArray1<usize>, np: usize, tol: f64)
        -> Lineshape1DMap {

        let n = x.len();

        Lineshape1DMap {
            x_map: x_map,

            lineshape: Lineshape1D::new(x, tol),
            
            y: Array::zeros((n, 2)),
            dydp: Array::zeros((np, n, 2)),
        }
    }

    //--------------------------------------
    pub fn eval(&mut self, p: &[f64], mut y: YSlice, mut dy: DySlice) {       

        // First, evaluate unique values
        let y_source = self.y.slice_mut(s![.., ..]);
        let dy_source = self.dydp.slice_mut(s![.., .., ..]);
        self.lineshape.eval(p, y_source, dy_source);

        // And then map unique values to repeats
        let np = p.len();

        // First, the y values themselves
        let iterator = self.x_map.iter()
            .zip(y.iter_mut().tuples::<(_,_)>());

        // Peaks are accumulated in Lineshape1D and simply assigned here
        for (i, (y_r, y_i)) in iterator {

            *y_r = self.y[[*i, 0]];
            *y_i = self.y[[*i, 1]];
        }

        // And then repeating the process for derivatives
        for i in 0 .. np {

            let temp = self.dydp.slice(s![i, .., ..]);
            let mut dydp = dy.slice_mut(s![i, .., ..]);

            let iterator = self.x_map.iter()
                .zip(dydp.iter_mut().tuples::<(_,_)>());

            for (j, (dydp_r, dydp_i)) in iterator {

                *dydp_r = temp[[*j, 0]];
                *dydp_i = temp[[*j, 1]];
            }
        }
    }
}


//==============================================================================
// General 2D lineshape

pub struct Lineshape2D {

    // Pairs of 1D lineshapes and ranges mapping parameters to resonances
    lineshapes: Vec<(Lineshape1DMap, Lineshape1DMap)>,
    ranges: Vec<(Range<usize>, Range<usize>)>,

    // Work arrays for individual dimensions
    y_direct: Array2<f64>,
    y_indirect: Array2<f64>,

    dy_direct: Array3<f64>,
    dy_indirect: Array3<f64>,
}


//--------------------------------------
impl Lineshape2D {

    //--------------------------------------
    // parameter_map maps every peak parameter to a unique resonance/dimension,
    // i.e, (1, 1), (1, 2) would say that parameter 1 is for resonance 1
    // in the direct dimension while parameter 2 is for resonance 1 in
    // the indirect dimension.
    pub fn new(
            x_direct: ArcArray1<f64>, 
            x_indirect: ArcArray1<f64>, 
            parameter_map: Vec<(usize, usize)>,
            tol: f64
        ) -> Lineshape2D {

        // Parameter map must be sorted to ensure correct grouping
        for w in parameter_map[..].windows(2) {
            let a = w[0];
            let b = w[1];

            if (a.0 > b.0) || (a.1 > b.1) {
                panic!("Parameter map must be sorted.")
            }
        }

        let n = x_direct.len();

        // Breaking down x values into a subset of unique values and indexes
        let (x_direct, x_direct_map) = Lineshape1DMap::map(x_direct);
        let (x_indirect, x_indirect_map) = Lineshape1DMap::map(x_indirect);

        // Counting unique resonances/dimensions for conversion into ranges
        // (parameter map is currently assumed to have at least one element)
        let mut iterator = parameter_map.into_iter();
        let mut last_value = iterator.next().unwrap();

        // Total number of parameters
        let mut np = 1; 

        // Parameters split by resonance/dimension
        let mut counts: Vec<usize> = vec![1]; 

        // Current index of resonance/dimension count
        let mut i = 0;

        while let Some(value) = iterator.next() {
            if (last_value.0 != value.0) || (last_value.1 != value.1) {
                counts.push(0);
                i += 1;
                last_value = value;
            } 

            np += 1;
            counts[i] += 1;
        }

        // Chunk counts into 2 (direct/indirect) and generate lineshapes
        let mut lineshapes: Vec<(Lineshape1DMap, Lineshape1DMap)> = Vec::new();

        for (np_direct, np_indirect) in counts.iter().tuples() {

            let direct_lineshape = Lineshape1DMap::new(
                x_direct.clone(), x_direct_map.clone(), *np_direct, tol
            );

            let indirect_lineshape = Lineshape1DMap::new(
                x_indirect.clone(), x_indirect_map.clone(), *np_indirect, tol
            );

            lineshapes.push((direct_lineshape, indirect_lineshape));
        }

        // Then converting counts into contiguous ranges
        let mut ranges: Vec<Range<usize>> = Vec::new();
        let mut i = 0;

        for count in counts.iter() {
            let j = i + count;
            ranges.push(i..j);
            i = j;
        }

        // Chunk and group
        let ranges: Vec<(Range<usize>, Range<usize>)> = ranges.into_iter()
            .tuples::<(_,_)>()
            .collect();

        Lineshape2D {
            lineshapes: lineshapes,
            ranges: ranges,

            y_direct: Array2::zeros((n, 2)),
            y_indirect: Array2::zeros((n, 2)),

            dy_direct: Array3::zeros((np, n, 2)),
            dy_indirect: Array3::zeros((np, n, 2)),
        }
    }
}


impl ComponentEval for Lineshape2D {


    //--------------------------------------
    fn nd(&self) -> usize {
        4
    }

    //--------------------------------------
    fn n(&self) -> usize {
        self.y_direct.shape()[0]
    }

    //--------------------------------------
    fn eval(&mut self, p: &[f64], mut y: YSlice, mut dy: DySlice) {       

        // Loop over the resonances
        let iterator = self.lineshapes.iter_mut()
            .zip(self.ranges.iter());

        // Initialize y to zero
        y.fill(0.0);

        for ((lineshape_1, lineshape_2), (range_1, range_2)) in iterator {

            // First, perform individual evaluations
            let y_direct = self.y_direct.slice_mut(s![.., ..]);
            let dy_direct = self.dy_direct.slice_mut(
                s![range_1.clone(), .., ..]
            );
            lineshape_1.eval(&p[range_1.clone()], y_direct, dy_direct);


            let y_indirect = self.y_indirect.slice_mut(s![.., ..]);
            let dy_indirect = self.dy_indirect.slice_mut(
                s![range_2.clone(), .., ..]
            );
            lineshape_2.eval(&p[range_2.clone()], y_indirect, dy_indirect);

            // Apply cross-product on a series of chunks
            // (size 2 chunks for 1D, r/i, size 4 chunks for 2D, rr/ri/ir/ii)
            
            // First y * y
            let iterator = izip!(
                self.y_direct.iter().tuples::<(_,_)>(),
                self.y_indirect.iter().tuples::<(_,_)>(),
                y.iter_mut().tuples::<(_,_,_,_)>()
            );

            for ((d_r, d_i), (i_r, i_i), o) in iterator {
            
                let (o_rr, o_ri, o_ir, o_ii) = o;

                // y is accumulated
                *o_rr += d_r * i_r;
                *o_ri += d_r * i_i;
                *o_ir += d_i * i_r;
                *o_ii += d_i * i_i;
            }

            // Then dy * y (cycling y as it's not the same length)
            let dy_direct = self.dy_direct.slice(
                s![range_1.clone(), .., ..]
            );
            let mut output = dy.slice_mut(s![range_1.clone(), .., ..]);
            
            let iterator = izip!(
                dy_direct.iter().tuples::<(_,_)>(),
                self.y_indirect.iter().tuples::<(_,_)>().cycle(),
                output.iter_mut().tuples::<(_,_,_,_)>()
            );

            for ((d_r, d_i), (i_r, i_i), o) in iterator {
            
                let (o_rr, o_ri, o_ir, o_ii) = o;

                // dy is set
                *o_rr = d_r * i_r;
                *o_ri = d_r * i_i;
                *o_ir = d_i * i_r;
                *o_ii = d_i * i_i;
            }


            // Then y * dy (cycling y as it's not the same length)
            let dy_indirect = self.dy_indirect.slice(
                s![range_2.clone(), .., ..]
            );
            let mut output = dy.slice_mut(s![range_2.clone(), .., ..]);
            
            let iterator = izip!(
                self.y_direct.iter().tuples::<(_,_)>().cycle(),
                dy_indirect.iter().tuples::<(_,_)>(),
                output.iter_mut().tuples::<(_,_,_,_)>()
            );

            for ((d_r, d_i), (i_r, i_i), o) in iterator {
            
                let (o_rr, o_ri, o_ir, o_ii) = o;

                // dy is set
                *o_rr = d_r * i_r;
                *o_ri = d_r * i_i;
                *o_ir = d_i * i_r;
                *o_ii = d_i * i_i;
            }
        }
    }
}
*/

//=============================================================================
// Unit tests


#[cfg(test)]
mod tests {

    use ndarray::prelude::*;
    use std::iter;

    use crate::fitting::common::ComponentEval;


    /*
    //--------------------------------------
    #[test]
    fn map() {

        let n = 5;
        let x = Array::linspace(0.0, 1.0, n).into_shared();

        let xi: Vec<usize> = (0..n).map(|x| iter::repeat(x).take(2))
            .flatten()
            .collect();
        let ni = xi.len();
        let xi = Array::from_shape_vec((ni,), xi).unwrap().into_shared();

        let p = vec![
            // peak
            0.5, 0.5, 0.5, 0.5, 
            // peak
            0.3, 0.3, 0.3, 0.3
        ];
        let np = p.len();

        // Unique
        let mut unique = super::Lineshape1D::new(x.clone(), np, 1e-10);
        let mut y_unique: Array2<f64> = Array::zeros((n, 2)); 
        let mut dy_unique: Array3<f64> = Array::zeros((np, n, 2)); 

        let y = y_unique.slice_mut(s![.., ..]);
        let dy = dy_unique.slice_mut(s![.., .., ..]);

        unique.eval(&p, y, dy);

        // Mapped
        let mut mapped = super::Lineshape1DMap::new(
            x.clone(), xi.clone(), np, 1e-10
        );
        let mut y_mapped: Array2<f64> = Array::zeros((ni, 2)); 
        let mut dy_mapped: Array3<f64> = Array::zeros((np, ni, 2)); 

        let y = y_mapped.slice_mut(s![.., ..]);
        let dy = dy_mapped.slice_mut(s![.., .., ..]);

        mapped.eval(&p, y, dy);

        for (i, j) in mapped.x_map.iter().enumerate() {
            assert!(y_unique[[*j, 0]] == y_mapped[[i, 0]], "map y mismatch"); 
            assert!(y_unique[[*j, 1]] == y_mapped[[i, 1]], "map y mismatch"); 
        }

        for k in 0 .. np {
            for (i, j) in mapped.x_map.iter().enumerate() {
                assert!(
                    dy_unique[[k, *j, 0]] == dy_mapped[[k, i, 0]], 
                    "map dydp mismatch"
                ); 

                assert!(
                    dy_unique[[k, *j, 1]] == dy_mapped[[k, i, 1]], 
                    "map dydp mismatch"
                ); 
            }
        }
    }
    */


    //--------------------------------------
    #[test]
    fn lineshape_1d_gradient() {

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

        let mut lineshape = super::Lineshape1D::new(
            x, vec![0..n], vec![0, 0], 1e-10
        );

        // Summing up errors increases threshold
        lineshape.grad_test(&p, 5e-6);
    }


    /*
    //--------------------------------------
    #[test]
    fn lineshape_2d_gradient() {

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

        let p = vec![
            // peak
            0.5, 0.5, 0.5, 0.5, 
            // peak
            0.3, 0.3, 0.3, 0.3
        ];

        // Place 2nd peak in indirect dimension
        let parameter_map: Vec<(usize, usize)> = iter::repeat((1,1)).take(4)
            .chain(iter::repeat((1,2)).take(4))
            .collect();

        let mut lineshape = super::Lineshape2D::new(
            x1, x2, parameter_map, 1e-10
        );

        // Summing up errors increases threshold
        lineshape.grad_test(&p, 5e-6);
    }
    */
}
