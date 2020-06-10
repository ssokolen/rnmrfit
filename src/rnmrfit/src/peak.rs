use num::complex::Complex;
use std::f64::consts::SQRT_2;

use faddeeva;

// Some constants
const SQRT_PI: f64 = 1.772453850905515881919427556567825376987457275390625_f64;

//==============================================================================
// Generic peak function traits

pub trait PeakFunctions {
    fn value(&mut self, x: f64) -> Complex<f64>;
    fn gradients(&mut self, x: f64, grad: &mut[Complex<f64>]) -> Complex<f64>;
}


//==============================================================================
// Generic peak structure

pub enum Peak {
    Lorentz(Lorentz),
    Voigt(Voigt),
}

impl Peak {

    //------------------------------------------------
    pub fn new(p: f64, w: f64, h: f64, f: f64) -> Peak {

        if f < 1e-8 {
            Peak::Lorentz( Lorentz::new(p, w, h) )
        } else if f < (1.0 - 1e-8) {
            Peak::Voigt( Voigt::new(p, w, h, f) )
        } else {
            panic!("Gauss peaks not currently supported.");
        }

    }

    //--------------------------------------
    pub fn from_vec(p: &[f64]) -> Peak {

        Peak::new(p[0], p[1], p[2], p[3])

    }

}

impl PeakFunctions for Peak {

    //--------------------------------------
    fn value(&mut self, x: f64) -> Complex<f64> {

        match self {
            Peak::Lorentz( lorentz ) => lorentz.value(x),
            Peak::Voigt( voigt ) => voigt.value(x),
        }
    }

    //--------------------------------------
    fn gradients(&mut self, x: f64, grad: &mut[Complex<f64>]) -> Complex<f64> {

        match self {
            Peak::Lorentz( lorentz ) => lorentz.gradients(x, grad),
            Peak::Voigt( voigt ) => voigt.gradients(x, grad),
        }  
    }

}

//==============================================================================
// Lorentz

pub struct Lorentz {

    // Main parameters
    p: f64,
    w: f64,
    h: f64,
    
    // Intermediate terms independent of x
    dzdp: f64,

    // Intermediate terms dependent on x
    z: f64,
    z_2: f64,
    yo: Complex<f64>

}

//------------------------------------------------------------------------------
// Constructors

impl Lorentz {

    //--------------------------------------
    pub fn new(p:f64, w:f64, h: f64) -> Lorentz {

        Lorentz {
            p: p,
            w: w,
            h: h,

            dzdp: 1.0 / w,

            z: 0.0,
            z_2: 0.01,
            yo: Complex::new(0.0, 0.0),
        }
    }

    //--------------------------------------
    pub fn from_vec(p: &[f64]) -> Lorentz {

        Lorentz::new(p[0], p[1], p[2])

    }

    //--------------------------------------
    fn intermediates(&mut self, x: f64) {

        self.z = (self.p - x) / self.w;
        self.z_2 = self.z * self.z;

        self.yo = Complex::new(1.0, self.z)  / (self.z_2 + 1.0);

    }
}

//------------------------------------------------------------------------------
// Calculations

impl PeakFunctions for Lorentz {

    //--------------------------------------
    fn value(&mut self, x: f64) -> Complex<f64> {

        // Generate common intermediates
        self.intermediates(x);

        // Final value
        self.yo * self.h
    }

    //--------------------------------------
    fn gradients(&mut self, x: f64, grad: &mut[Complex<f64>]) -> Complex<f64> {

        // Generate common intermediates 
        self.intermediates(x);

        // Derivative of y (for chain rule)
        let dydz = Complex::new(-2.0 * self.z, 1.0 - self.z_2) * 
                   self.h / ( (self.z_2 + 1.0) * (self.z_2 + 1.0) );

        // Derivatives of z (for chain rule)
        let dzdw = -self.z / self.w;

        // Position gradient
        grad[0] = dydz * self.dzdp;

        // Lorentz width gradient
        grad[1] = dydz * dzdw;

        // Height gradient (consistent at yo)
        grad[2] = self.yo;

        // Final value
        self.yo * self.h     
    }

}

//==============================================================================
// Voigt

pub struct Voigt {

    // Main parameters
    p: f64,
    w: f64,
    _h: f64,
    _f: f64,
    wg: f64,
    
    // Intermediate terms independent of x
    yn: Complex<f64>,

    a: Complex<f64>,
    b: f64,
    
    dyndzn: Complex<f64>,
    dzndf: Complex<f64>,

    dzdp: f64,

    // Intermediate terms dependent on x
    z: Complex<f64>,
    yo: Complex<f64>,
}

//------------------------------------------------------------------------------
// Constructors

impl Voigt {

    //--------------------------------------
    pub fn new(p:f64, w:f64, h: f64, f: f64) -> Voigt {

        let wg = w * f / (1.0 - f);
        let zn = Complex::new( 0.0, w / (SQRT_2 * wg) );
        let yn = faddeeva::w( zn, 1e-6 );

        let a = h / yn;
        let b = wg / (w * f * f);

        let dyndzn = -2.0 * zn * yn + Complex::new( 0.0, 2.0/SQRT_PI );  
        let dzndf = -zn * b;

        let dzdp = 1.0/( SQRT_2 * wg );

        Voigt {
            p: p,
            w: w,
            _h: h,
            _f: f,
            wg: wg,

            yn: yn,

            a: a,
            b: b,
            
            dyndzn: dyndzn,
            dzndf: dzndf,

            dzdp: dzdp,

            z: Complex::new(0.0, 0.0),
            yo: Complex::new(0.0, 0.0),
        }
    }

    //--------------------------------------
    pub fn from_vec(p: &[f64]) -> Voigt {

        Voigt::new(p[0], p[1], p[2], p[3])

    }

    //--------------------------------------
    fn intermediates(&mut self, x: f64) {

        self.z =  Complex::new( self.p - x, self.w )  / ( SQRT_2 * self.wg);
        self.yo = faddeeva::w( self.z, 1e-6 );

    }
}

//------------------------------------------------------------------------------
// Calculations

impl PeakFunctions for Voigt {

    //--------------------------------------
    fn value(&mut self, x: f64) -> Complex<f64> {

        // Generate common intermediates
        self.intermediates(x);

        // Final value
        self.yo * self.a
    }

    //--------------------------------------
    fn gradients(&mut self, x: f64, grad: &mut[Complex<f64>]) -> Complex<f64> {

        // Generate common intermediates 
        self.intermediates(x);

        // Derivative of y (for chain rule)
        let dyodz = -2.0 * self.z * self.yo + Complex::new( 0.0, 2.0/SQRT_PI );  

        // Derivatives of z (for chain rule)
        let dzdw = -(self.p - x)/(SQRT_2 * self.w * self.wg);
        let dzdf = -self.z * self.b;

        // Position gradient
        grad[0] = self.a * dyodz * self.dzdp;

        // Lorentz width gradient
        grad[1] = self.a * dyodz * dzdw;

        // Height gradient
        grad[2] = self.yo / self.yn;

        // Fraction gradient
        grad[3] = self.a * (dyodz * dzdf * self.yn - 
                            self.dyndzn * self.dzndf * self.yo) / self.yn;

        // Final value
        self.yo * self.a     
    }

}

#[cfg(test)]
mod tests {
    
    use num::complex::Complex;
    use nlopt;

    use super::PeakFunctions;
    use super::Peak;

    fn check_gradient(p: Vec<f64>) {

        // x must be hard coded due to nlopt
        let x: f64 = 0.45;

        // Analytical gradients
        let mut peak = Peak::from_vec(&p);
        let mut grad = vec![Complex::new(0.0, 0.0); 4];
        peak.gradients(x, &mut grad);

        // Comparing to real numerical gradients
        fn f_real(p: &[f64]) -> f64 {
            let x = 0.45;
            let mut peak = Peak::from_vec(p);
            peak.value(x).re
        };

        let mut grad_real = vec![0.0; 4];
        nlopt::approximate_gradient(&p, f_real, &mut grad_real);

        assert!((grad[0].re - grad_real[0]).abs() < 1e-5, 
                "Position gradient error -- real domain");
        assert!((grad[1].re - grad_real[1]).abs() < 1e-5, 
                "Width gradient error -- real domain");
        assert!((grad[2].re - grad_real[2]).abs() < 1e-5, 
                "Height gradient error -- real domain");
        assert!((grad[3].re - grad_real[3]).abs() < 1e-5, 
                "Fraction gradient error -- real domain");

        // Comparing to maginary numerical gradients
        fn f_imag(p: &[f64]) -> f64 {
            let x = 0.45;
            let mut peak = Peak::from_vec(p);
            peak.value(x).im
        };

        let mut grad_imag = vec![0.0; 4];
        nlopt::approximate_gradient(&p, f_imag, &mut grad_imag);

        assert!((grad[0].im - grad_imag[0]).abs() < 1e-5, 
                "Position gradient error -- imaginary domain");
        assert!((grad[1].im - grad_imag[1]).abs() < 1e-5, 
                "Width gradient error -- imaginary domain");
        assert!((grad[2].im - grad_imag[2]).abs() < 1e-5, 
                "Height gradient error -- imaginary domain");
        assert!((grad[3].re - grad_real[3]).abs() < 1e-5, 
                "Fraction gradient error -- imaginary domain");

    }

    #[test]
    fn lorentz_gradient() {

        let p = vec![0.5, 0.5, 0.5, 0.0];
        check_gradient(p);

    }

    #[test]
    fn voigt_gradient() {

        let p = vec![0.5, 0.5, 0.5, 0.5];
        check_gradient(p);

    }
}
            
