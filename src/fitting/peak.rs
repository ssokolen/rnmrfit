use itertools::{
    izip,
    Itertools
};
use rustfft::num_complex::Complex;
use std::f64::consts::{PI, SQRT_2};

use crate::faddeeva;

// Some constants
const SQRT_PI: f64 = 1.772453850905515881919427556567825376987457275390625_f64;


//=============================================================================
// Generic peak functions


//--------------------------------------
pub fn f_peak_area(p: &[f64], tol: f64) -> f64 {

    let f = p[3];
    
    if f < 1e-8 {
        f_lorentz_area(p)
    } else if f < (1.0 - 1e-8) {
        f_voigt_area(p, tol)
    //} else if f <= 1.0 {
    //    panic!("Gauss peaks not currently supported.")
    } else {
        f64::NAN
    }
}


//--------------------------------------
pub fn f_peak<'a, X, Y>(p: &[f64], x: X, y: Y, tol: f64)
        where X: IntoIterator<Item = &'a f64>,
              Y: IntoIterator<Item = &'a mut f64> {

    let f = p[3];
    
    if f < 1e-8 {
        f_lorentz(p, x, y)
    } else if f < (1.0 - 1e-8) {
        f_voigt(p, x, y, tol)
    //} else if f <= 1.0 {
    //    panic!("Gauss peaks not currently supported.")
    } else {
        panic!("Invalid fraction Gauss.")
    }
}


//--------------------------------------
pub fn f_peak_gradients<'a, X, Y>(
        p: &[f64], x: X, y: Y, g: (Y, Y, Y, Y), tol: f64
    )
        where X: IntoIterator<Item = &'a f64>,
              Y: IntoIterator<Item = &'a mut f64> {

    let f = p[3];
    
    if f < 1e-8 {
        f_lorentz_gradients(p, x, y, (g.0, g.1, g.2))
    } else if f < (1.0 - 1e-8) {
        f_voigt_gradients(p, x, y, g, tol)
    } else {
        panic!("Gauss peaks not currently supported.")
    }
}


//=============================================================================
// Lorentz


//--------------------------------------
fn f_lorentz_area(p: &[f64]) -> f64 {

    let h = p[2];
    let w = p[1];
    //let p = p[0];

    PI * w * h
}


//--------------------------------------
fn f_lorentz<'a, X, Y>(p: &[f64], x: X, y: Y)
        where X: IntoIterator<Item = &'a f64>,
              Y: IntoIterator<Item = &'a mut f64> {

    let h = p[2];
    let w = p[1];
    let p = p[0];

    // y is broken down into re/im pairs
    let iterator = x.into_iter()
        .zip(y.into_iter().tuples::<(_, _)>());

    for (x, (y_r, y_i)) in iterator {

        let z = (p - x) / w;
        let y = h * Complex::new(1.0, z)  / (z*z + 1.0);

        *y_r = y.re;
        *y_i = y.im;
    }
}


//--------------------------------------
fn f_lorentz_gradients<'a, X, Y>(p: &[f64], x: X, y: Y, g: (Y, Y, Y))
        where X: IntoIterator<Item = &'a f64>,
              Y: IntoIterator<Item = &'a mut f64> {

    let h = p[2];
    let w = p[1];
    let p = p[0];

    let (gp, gw, gh) = g;

    // Intermediate terms
    let dzdp = 1.0 / w;

    // y and all the gradients are broken down into re/im pairs
    let iterator = izip!(
        x.into_iter(),
        y.into_iter().tuples::<(_,_)>(),
        gp.into_iter().tuples::<(_,_)>(),
        gw.into_iter().tuples::<(_,_)>(),
        gh.into_iter().tuples::<(_,_)>(),
    );

    for (x, (y_r, y_i), gp, gw, gh) in iterator {

        let (dydp_r, dydp_i) = gp;
        let (dydw_r, dydw_i) = gw;
        let (dydh_r, dydh_i) = gh;

        let z = (p - x) / w;
        let z_sq = z * z;
        let yo = Complex::new(1.0, z)  / (z_sq + 1.0);

        // Derivative of y
        let dydz = Complex::new(-2.0 * z, 1.0 - z_sq) * 
                   h / ( (z_sq + 1.0) * (z_sq + 1.0) );

        // Position gradient
        let dydp = dydz * dzdp;

        *dydp_r = dydp.re;
        *dydp_i = dydp.im;

        // Lorentz width gradient
        let dzdw = -z / w;
        let dydw = dydz * dzdw;

        *dydw_r = dydw.re;
        *dydw_i = dydw.im;

        // Height gradient
        *dydh_r = yo.re;
        *dydh_i = yo.im;

        // Value
        let y = h * yo;

        *y_r = y.re;
        *y_i = y.im; 
    }
}


//=============================================================================
// Voigt


//--------------------------------------
fn f_voigt_area(p: &[f64], tol: f64) -> f64 {

    let f = p[3];
    let h = p[2];
    let w = p[1];
    //let p = p[0];

    let wg = w * f / (1.0 - f);

    let zn = Complex::new( 0.0, w / (SQRT_2 * wg) );
    let yn = faddeeva::w( zn, tol );

    let area = SQRT_2 * SQRT_PI * wg * h / yn;
    area.re
}


//--------------------------------------
fn f_voigt<'a, X, Y>(p: &[f64], x: X, y: Y, tol: f64)
        where X: IntoIterator<Item = &'a f64>,
              Y: IntoIterator<Item = &'a mut f64> {

    let f = p[3];
    let h = p[2];
    let w = p[1];
    let p = p[0];

    let wg = w * f / (1.0 - f);

    // Intermediate terms
    let zn = Complex::new( 0.0, w / (SQRT_2 * wg) );
    let yn = faddeeva::w( zn, tol );

    // y is broken down into re/im pairs
    let y_chunks = y.into_iter().chunks(2);

    let iterator = x.into_iter()
        .zip(y_chunks.into_iter());

    for (x, mut y) in iterator {

        let y_re = y.next().unwrap();
        let y_im = y.next().unwrap();

        let z = Complex::new( p - x, w )  / ( SQRT_2 * wg);
        let y = h / yn * faddeeva::w( z, tol );

        *y_re = y.re;
        *y_im = y.im;
    }
}


//--------------------------------------
fn f_voigt_gradients<'a, X, Y>(
        p: &[f64], x: X, y: Y, g: (Y, Y, Y, Y), tol: f64
    )
        where X: IntoIterator<Item = &'a f64>,
              Y: IntoIterator<Item = &'a mut f64> {

    let f = p[3];
    let h = p[2];
    let w = p[1];
    let p = p[0];

    let wg = w * f / (1.0 - f);

    let (gp, gw, gh, gf) = g;

    // Intermediate terms
    let zn = Complex::new( 0.0, w / (SQRT_2 * wg) );
    let yn = faddeeva::w( zn, tol );

    let a = h / yn;
    let b = wg / (w * f * f);

    let dyndzn = -2.0 * zn * yn + Complex::new( 0.0, 2.0/SQRT_PI );  
    let dzndf = -zn * b;

    let dzdp = 1.0/( SQRT_2 * wg );

    // y is broken down into re/im pairs along with the gradients
    let iterator = izip!(
        x.into_iter(),
        y.into_iter().tuples::<(_,_)>(),
        gp.into_iter().tuples::<(_,_)>(),
        gw.into_iter().tuples::<(_,_)>(),
        gh.into_iter().tuples::<(_,_)>(),
        gf.into_iter().tuples::<(_,_)>(),
    );

    for (x, (y_r, y_i), gp, gw, gh, gf) in iterator {

        let (dydp_r, dydp_i) = gp;
        let (dydw_r, dydw_i) = gw;
        let (dydh_r, dydh_i) = gh;
        let (dydf_r, dydf_i) = gf;

        let z =  Complex::new( p - x, w )  / ( SQRT_2 * wg);
        let yo = faddeeva::w( z, tol );

        // Derivative of y
        let dyodz = -2.0 * z * yo + Complex::new( 0.0, 2.0/SQRT_PI );  

        // Position gradient
        let dydp = a * dyodz * dzdp;

        *dydp_r = dydp.re;
        *dydp_i = dydp.im;

        // Lorentz width gradient
        let dzdw = (x - p)/(SQRT_2 * w * wg);
        let dydw = a * dyodz * dzdw;

        *dydw_r = dydw.re;
        *dydw_i = dydw.im;

        // Height gradient
        let dydh = yo / yn;

        *dydh_r = dydh.re;
        *dydh_i = dydh.im;

        // Fraction gradient
        let dzdf = -z * b;
        let dydf = a * (dyodz * dzdf * yn - 
                        dyndzn * dzndf * yo) / yn;

        *dydf_r = dydf.re;
        *dydf_i = dydf.im;

        // Value
        let y = a * yo;

        *y_r = y.re;
        *y_i = y.im; 
    }
}


//=============================================================================


#[cfg(test)]
mod tests {
    
    use itertools::{Itertools, izip};
    use nlopt;

    use super::{f_peak, f_peak_gradients};

    fn check_gradient(p: &[f64]) {

        // x must be hard coded due to nlopt
        let x = vec![0.45];

        // Analytical gradients

        let mut y = vec![0.0; 2];

        let mut gp = vec![0.0; 2];
        let mut gw = vec![0.0; 2];
        let mut gh = vec![0.0; 2];
        let mut gf = vec![0.0; 2];

        let g = (&mut gp[..], &mut gw[..], &mut gh[..], &mut gf[..]);

        f_peak_gradients(p, &x[..], &mut y[..], g, 1e-10);
        println!("{:?}", vec![gp.clone(), gw.clone(), gh.clone(), gf.clone()]);

        // Numerical gradients

        fn f_real(p: &[f64]) -> f64 {
            let x = vec![0.45];
            let mut y = vec![0.0, 0.0];
            f_peak(p, &x[..], &mut y[..], 1e-10);
            y[0]
        }

        let mut grad_real = vec![0.0; 4];
        nlopt::approximate_gradient(&p, f_real, &mut grad_real);
        println!("{:?}", grad_real);

        fn f_imag(p: &[f64]) -> f64 {
            let x = vec![0.45];
            let mut y = vec![0.0, 0.0];
            f_peak(p, &x[..], &mut y[..], 1e-10);
            y[1]
        }

        let mut grad_imag = vec![0.0; 4];
        nlopt::approximate_gradient(&p, f_imag, &mut grad_imag);
        println!("{:?}", grad_imag);

        // Comparison

        let domains = vec!["real", "imaginary"];
        let parameters = vec!["position", "width", "height", "fraction"];

        let mut gradients = vec![gp, gw, gh];

        if p[3] > 1e-6 {
            gradients.push(gf);
        }

        let iterator = izip!(
            gradients.into_iter().flatten(),
            grad_real.into_iter()
                .interleave(grad_imag.into_iter()),
            domains.into_iter().cycle(),
            parameters.into_iter().cycle()
        );

        for (dy_true, dy_estimate, domain, parameter) in iterator {

            assert!((dy_true - dy_estimate).abs() < 1e-6, 
                    "Gradient error -- {} {} ({} vs {} true)", 
                    domain, parameter, dy_estimate, dy_true);

        }
    }


    #[test]
    fn lorentz_gradient() {

        let p = vec![0.5, 0.5, 0.5, 0.0];
        check_gradient(&p[..]);
    }


    #[test]
    fn voigt_gradient() {

        let p = vec![0.5, 0.5, 0.5, 0.5];
        check_gradient(&p[..]);
    }


    #[test]
    fn voigt_repeatability() {

        let p = vec![0.5, 0.5, 0.5, 0.5];
        let x = vec![0.45];

        
        let diffs: f64 = (0 .. 1_000)
            .map(|_| {

                let mut y = vec![0.0, 0.0];
                f_peak(&p, &x[..], &mut y[..], 1e-6);
                let y1 = y[0];

                let mut y = vec![0.0, 0.0];
                f_peak(&p, &x[..], &mut y[..], 1e-6);
                let y2 = y[0];

                (y2 - y1).abs()
            })
            .sum();

        println!("{}", diffs);

        assert!(diffs == 0.0);
    }
}
            
