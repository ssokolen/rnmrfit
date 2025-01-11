use anyhow::{Result, anyhow, bail};
use cached::proc_macro::cached;
use cascade::cascade;
use getset::{Getters, Setters};
use itertools::{Itertools, izip};
use itertools_num;
use ndarray::{
    prelude::*,
    concatenate,
};
use rustfft::num_complex::Complex;
use serde::Deserialize;
use std::{
    iter,
    convert::TryInto,
    ops::Range,
};

use crate::fitting;
use crate::processing::Bound;


//=============================================================================
// Template

#[derive(Clone, Debug, Default, Deserialize, Getters)]
#[serde(deny_unknown_fields)]
#[getset(get = "pub")]
pub struct Baseline1DTemplate {

    #[serde(default)]
    span: Option<f64>,
    #[serde(default)]
    degree: Option<usize>,
    #[serde(default)]
    bounds: Bound,

}


impl Baseline1DTemplate {

    pub fn update(&mut self, template: &Baseline1DTemplate) {

        if self.span.is_none() {
            self.span = *template.span();
        }

        if self.degree.is_none() {
            self.degree = *template.degree();
        }

        self.bounds.update(template.bounds())
    }
}


// Sensible global defaults
pub const DEFAULT_BASELINE: Baseline1DTemplate = 
    Baseline1DTemplate {
    
        span: Some(f64::NAN),
        degree: Some(3),
        bounds: Bound { lower: -1.0, upper: 1.0 }

};


//=============================================================================
// Baseline1D

#[derive(Clone, Debug, Getters, Setters)]
pub struct Baseline1D {

    knots: Vec<f64>,
    degree: usize,

    terms: Vec<f64>,
    
    #[getset(get = "pub", set = "pub")]
    bounds: Bound,

    #[getset(get = "pub")]
    np: usize,
}


//=============================================================================
// Constructors

impl Baseline1D {

    //-------------------------------------------------------------------------
    /// Baseline1D from a span length and x data
    pub fn from_span(
            x: &[f64], span: f64, degree: usize
        ) -> Result<Baseline1D> { 

        // Rounding span based on evenly spread knots
        let first_x = x.first().expect("Error generating baseline.");
        let last_x = x.last().expect("Error generating baseline.");

        assert!(first_x > last_x, "Incorrect x sort order.");

        // Default n is 100
        let n = if span.is_nan() {
            20.0
        } else {
            ( (first_x - last_x)/span ).ceil()
        };

        let span = (first_x - last_x)/n;

        let knots: Vec<f64> = (0..)
            .map(|i| last_x + (i as f64) * span)
            .take_while(|x| x <= first_x)
            .collect();

        Baseline1D::from_knots(knots, degree)
    }


    //-------------------------------------------------------------------------
    /// Baseline1D from a complete set of knots (boundary and internal)
    pub fn from_knots(
            knots: Vec<f64>, degree: usize
        ) -> Result<Baseline1D> { 

        let first_knot = *knots.first()
            .ok_or(anyhow!("Cannot process knot values"))?;
        let last_knot = *knots.last()
            .ok_or(anyhow!("Cannot process knot values"))?;

        // Repeating boundary knots
        let knots: Vec<f64> = iter::repeat(first_knot).take(degree)
            .chain(knots.into_iter())
            .chain(iter::repeat(last_knot).take(degree))
            .collect();

        let np = knots.len() - (degree + 2) + 1;

        let baseline = Baseline1D {
            knots: knots,
            degree: degree,
            terms: vec![0.0; np*2],
            bounds: Bound { lower: -1.0, upper: 1.0 },
            np: np,
        };

        Ok( baseline )
    }
}


//=============================================================================
// Methods

impl Baseline1D {

    //-------------------------------------------------------------------------
    pub fn terms(&self, ppm: &[ [f64; 2] ]) -> Vec<f64> {

        self.term_ranges(ppm).into_iter()
            .map(|range| self.terms[range].to_vec())
            .flatten()
            .collect()
    }


    //-------------------------------------------------------------------------
    pub fn set_terms(&mut self, terms: &[f64], ppm: &[ [f64; 2] ]) {

        let indexes: Vec<usize> = self.term_ranges(ppm)
            .into_iter()
            .map(|range| range.into_iter())
            .flatten()
            .collect();

        if terms.len() != indexes.len() {
            panic!("Invalid number of terms.")
        }

        for (i, value) in indexes.into_iter().zip(terms.into_iter()) {
            self.terms[i] = *value;
        }
    }


    //-------------------------------------------------------------------------
    pub fn basis(&self, x: &[f64], ppm: &[ [f64; 2] ]) -> Array2<f64> {

        // Generate one basis per range
        let baseis: Vec<_> = self.knot_ranges(ppm)
            .into_iter()
            .map(|range| {
                let knots = self.knots[range].to_vec();
                generate_basis(x, &knots, self.degree)
            })
            .collect();

        // Concatenation requires views
        let views: Vec::<_> = baseis.iter()
            .map(|basis| basis.slice(s![..,..]))
            .collect();

        // Then combine them column-wise
        concatenate(Axis(1), &views).expect("Error generating basis.")
    }


    //-------------------------------------------------------------------------
    pub fn baseline(&self, x: &[f64]) -> Vec<Complex<f64>> {

        // Real and imaginary components share the same basis
        let basis = generate_basis(
            x, &self.knots, self.degree
        );

        let n = x.len();

        // But values are calculated separately
        let mut components: Vec<Vec<f64>> = vec![];

        for offset in 0..2 {

            let terms: Vec<_> = self.terms.iter()
                .skip(offset)
                .step_by(2)
                .copied()
                .collect();

            let terms = Array::from(terms);
            let y = basis.dot(&terms);
            
            components.push(y.into_iter().collect());
        }

        // Then, combine components
        components[0].iter()
            .zip(components[1].iter())
            .map(|(r, i)| Complex::new(*r, *i))
            .collect()
    }


    //-------------------------------------------------------------------------
    // Match term indexes to ppm
    fn knot_range(&self, ppm: &[f64; 2]) -> Range<usize> {

        // Find first knot bigger than ppm cutoff
        // (this should avoid issues with repeating lower bounds)
        let mut i = self.knots.iter()
            .enumerate()
            .skip_while(|(_, x)| **x <= ppm[0])
            .next()
            .map(|(i, _)| i)
            .expect("Invalid x domain for baseline.");

        // And subtract one to get lower bound
        i -= 1;
        i = i.max(0);

        // Then factor in degrees
        i -= self.degree;

        // Repeat for upper bounds
        let mut j = self.knots.iter()
            .enumerate()
            .skip_while(|(_, x)| **x < ppm[1])
            .next()
            .map(|(i, _)| i)
            .expect("Invalid x domain for baseline.");

        // Factor in degrees
        j += self.degree;

        // Making sure to actually include j in the range...
        i..(j+1)
    }

    //-------------------------------------------------------------------------
    // Combine terms from potentially overlapping ranges
    fn knot_ranges(&self, ppm: &[ [f64; 2] ]) -> Vec<Range<usize>> {

        let mut ranges = ppm.iter()
            .map(|range| self.knot_range(range))
            .sorted_by(|a, b| (a.start).cmp(&b.start));

        let first = ranges.next().expect("Error segmenting baseline.");

        // Combining
        ranges.into_iter()
            .fold(vec![first], |mut a, r| {
                let n = a.len();

                if r.start <= a[n-1].end {
                    a[n-1].end = a[n-1].end.max(r.end);
                } else {
                    a.push(r)
                }

                a
            })
    }


    //-------------------------------------------------------------------------
    // Match term indexes to ppm
    fn term_ranges(&self, ppm: &[ [f64; 2] ]) -> Vec<Range<usize>> {

        // Start with the knot ranges
        let knot_ranges = self.knot_ranges(ppm);

        let mut term_ranges = knot_ranges.into_iter()
            .map(|range| {
                // The start point maps directly to term index.
                // All we need to do is multipy by 2 to account for 
                // real/imaginary.
                let i = range.start * 2;

                // To get the end point, we calculate how many more knots
                // there are on the right, and use the same offset for terms
                // (multiplied by 2 to account for real/imaginary).
                let offset = self.knots.len() - range.end;
                let j = self.terms.len() - offset*2;

                // Since j starts at length and not length - 1, it already
                // accounts for the fact that the actual index j is not 
                // included in the range.
                i..j
            });

        let first = term_ranges.next().expect("Error segmenting baseline.");

        // Combining
        term_ranges.into_iter()
            .fold(vec![first], |mut a, r| {
                let n = a.len();

                if r.start <= a[n-1].end {
                    a[n-1].end = a[n-1].end.max(r.end);
                } else {
                    a.push(r)
                }

                a
            })
    }

    /*
    //-------------------------------------------------------------------------
    // Combine terms from potentially overlapping ranges
    fn term_ranges(&self, ppm: &[ [f64; 2] ]) -> Vec<Range<usize>> {

        let mut ranges = ppm.iter()
            .map(|range| self.term_range(range))
            .sorted_by(|a, b| (a.start).cmp(&b.start));

        let first = ranges.next().expect("Error segmenting baseline.");

        // Combining
        ranges.into_iter()
            .fold(vec![first], |mut a, r| {
                let n = a.len();

                if r.start <= a[n-1].end {
                    a[n-1].end = a[n-1].end.max(r.end);
                } else {
                    a.push(r)
                }

                a
            })
    }
    */
}


//=============================================================================
// Basis function

// Based on https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/
// bspline-ex-1.html
// Knots have to be sorted smallest to biggest and repeated if necessary
fn generate_basis(x: &[f64], knots: &[f64], degree: usize) -> Array2<f64> {

    let first_knot = knots.first().expect("Error generating basis");
    let last_knot = knots.last().expect("Error generating basis");

    assert!(first_knot < last_knot, "Incorrect knot sort order");

    // Each basis column/parameter is determined by degree + 2 knots.
    // For reference, boundary knots need to be repeated degree + 1 times,
    // i.e., adding an extra degree repeats to an existing boundary.
    let np = knots.len() - (degree + 2) + 1;
    assert!(np > 0, "Insufficient knots to generate basis");
    
    let n = x.len();

    let mut matrix = Array::zeros((n, np));

    for j in 0 .. np {
        let mut column = matrix.slice_mut(s![..,j]);
        for (i, value) in column.iter_mut().enumerate() {
            *value = basis(x[i], j, degree, &knots);
        }
    }

    // Note, the very last matrix parameter is always 0 since the recursuve
    // function is not inclusive of the last knot. However, this also 
    // technically means that the last x value becomes undefined and the
    // overall matrix becomes singular... For the sake of convenience, 
    // we can set the last parameter to 1, therefore "extrapolating"
    // the final spline to account for the last x value.
    // Note, where the last value is depends on x order.
    if x[0] > x[n-1] {

        matrix[[0, np-1]] = 1.0;
    
    } else {
        
        matrix[[n-1, np-1]] = 1.0;

    }

    matrix 
}


fn weight(x: f64, i: usize, j: usize, knots: &[f64]) -> f64 {
    if knots[i+j] != knots[i] {
        (x - knots[i])/(knots[i+j] - knots[i])
    } else {
        0.0
    }
}


fn basis(x: f64, i: usize, j: usize, knots: &[f64]) -> f64 {
    if j == 0 {
        if (x >= knots[i]) && (x < knots[i+1]) {
            1.0
        } else {
            0.0
        }
    } else {
        weight(x, i, j, knots)*basis(x, i, j-1, knots) +
        ( 1.0 - weight(x, i+1, j, knots) )*basis(x, i+1, j-1, knots)
    }
}


//=============================================================================
// Unit tests

#[cfg(test)]
mod tests {

    use rustfft::num_complex::Complex;

    //--------------------------------------
    #[test]
    fn baseline_basis() {

        // Comparing values to output of R bSpline (from splines2)

        let x: Vec<f64> = vec![0.0, 0.25, 0.5, 0.75, 1.0];
        let knots: Vec<f64> = vec![0.0, 0.5, 1.0];
        let basis = super::generate_basis(&x, &knots, 0);

        assert!( (basis[[1, 0]] - 1.0).abs() < 1e-8 ); 
        assert!( (basis[[2, 0]] - 0.0).abs() < 1e-8 ); 

        assert!( (basis[[1, 1]] - 0.0).abs() < 1e-8 ); 
        assert!( (basis[[3, 1]] - 1.0).abs() < 1e-8 ); 


        let x: Vec<f64> = vec![2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 0.8, 0.6, 0.4];
        let knots: Vec<f64> = vec![
            0.4, 0.4, 0.4, 0.4, 0.4, 
            0.4, 0.5, 1.0, 1.5, 2.0,
            2.0, 2.0, 2.0, 2.0, 2.0
        ];
        let basis = super::generate_basis(&x, &knots, 5);

        assert!( (basis[[7, 1]] - 0.158024691).abs() < 1e-8 ); 
        assert!( (basis[[8, 1]] - 0.0).abs() < 1e-8 ); 

        assert!( (basis[[0, 3]] - 0.0).abs() < 1e-8 ); 
        assert!( (basis[[1, 3]] - 0.0001666667).abs() < 1e-8 ); 

        assert!( (basis[[0, 8]] - 1.0).abs() < 1e-8 ); 
        assert!( (basis[[2, 8]] - 0.00032).abs() < 1e-8 ); 
        assert!( (basis[[3, 8]] - 0.0).abs() < 1e-8 ); 
    }    


    //--------------------------------------
    #[test]
    fn baseline_term_extraction() {

        let x: Vec<f64> = vec![0.0, 0.25, 0.5, 0.75, 1.0].into_iter()
            .rev()
            .collect();
        let knots: Vec<f64> = vec![0.0, 0.5, 1.0];

        let mut baseline = super::Baseline1D::from_knots(knots.clone(), 1)
            .unwrap();

        baseline.terms = (0..baseline.terms.len())
            .map(|x| x as f64)
            .collect();

        assert!(baseline.terms(&[[0.6, 0.9]]) == baseline.terms[2..].to_vec());
    }


    //--------------------------------------
    #[test]
    fn baseline_basis_extraction() {

        let x: Vec<f64> = vec![0.0, 0.25, 0.5, 0.75, 1.0].into_iter()
            .rev()
            .collect();
        let knots: Vec<f64> = vec![0.0, 0.5, 1.0];

        let mut baseline = super::Baseline1D::from_knots(knots.clone(), 1)
            .unwrap();

        let actual = baseline.basis(&x, &[[0.5, 1.0]]);

        // Generating basis from scratch
        let knots: Vec<f64> = vec![0.0, 0.5, 1.0, 1.0];

        let expected = super::generate_basis(&x, &knots, 1);

        assert!(expected == actual);
    }


    //--------------------------------------
    #[test]
    fn baseline_baseline_values() {

        let x: Vec<f64> = vec![0.0, 0.25, 0.5, 0.75, 1.0].into_iter()
            .rev()
            .collect();
        let knots: Vec<f64> = vec![0.0, 0.5, 1.0];

        let mut baseline = super::Baseline1D::from_knots(knots.clone(), 1)
            .unwrap();

        baseline.terms = (0..baseline.terms.len())
            .map(|x| x as f64)
            .collect();

        let x: Vec<f64> = vec![0.0, 0.25, 0.75, 1.0].into_iter()
            .rev()
            .collect();

        let expected: Vec<Complex<f64>> = vec![4.0, 3.0, 1.0, 0.0].into_iter()
            .zip(vec![5.0, 4.0, 2.0, 1.0].into_iter())
            .map(|(r, i)| Complex::new(r, i))
            .collect();

        let actual = baseline.baseline(&x);

        assert!(expected == actual);
    }
}
