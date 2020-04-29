

//==============================================================================
// C defintions

//------------------------------------------------------------------------------
// Ensuring C-compliant complex numbers

use num::complex::Complex;

#[derive(PartialEq, Copy, Clone, Hash, Debug, Default)]
#[repr(C)]
struct CComplex<T> {
    pub re: T,
    pub im: T,
}

impl<T> CComplex<T> {
    
    pub fn to_complex(self) -> Complex<T> {
        Complex::new(self.re, self.im)
    }

    pub fn from_complex(z: Complex<T>) -> CComplex<T> {
        CComplex { 
            re: z.re, 
            im: z.im,
        }
    }

}

//------------------------------------------------------------------------------
// Function definitions

extern "C" {
    fn Faddeeva_w(z: CComplex<f64>, relerr: f64) -> CComplex<f64>;
    //pub fn Faddeeva_w_im(x: f64) -> f64;
    //pub fn Faddeeva_erfcx(z: CComplex<f64>, relerr: f64) -> CComplex<f64>;
    //pub fn Faddeeva_erfcx_re(x: f64) -> f64;
    //pub fn Faddeeva_erf(z: CComplex<f64>, relerr: f64) -> CComplex<f64>;
    //pub fn Faddeeva_erf_re(x: f64) -> f64;
    //pub fn Faddeeva_erfi(z: CComplex<f64>, relerr: f64) -> CComplex<f64>;
    //pub fn Faddeeva_erfi_re(x: f64) -> f64;
    //pub fn Faddeeva_erfc(z: CComplex<f64>, relerr: f64) -> CComplex<f64>;
    //pub fn Faddeeva_erfc_re(x: f64) -> f64;
    //pub fn Faddeeva_Dawson(z: CComplex<f64>, relerr: f64) -> CComplex<f64>;
    //pub fn Faddeeva_Dawson_re(x: f64) -> f64;
}


//==============================================================================
// 

pub fn w(z: Complex<f64>, relerr: f64) -> Complex<f64> {
    let out = unsafe { 
        Faddeeva_w(CComplex::from_complex(z), relerr)
    };
    out.to_complex()
}

//==============================================================================
// Testing

#[cfg(test)]
mod tests {

    use assert_approx_eq::assert_approx_eq;
    use num::complex::Complex;

    use super::w;

    #[test]
    fn faddeeva_w() {

        let z: Complex<f64> = Complex::new(1.0, 1.0);
        let out: Complex<f64> = w(z, 1e-8);

        assert_approx_eq!(out.re, 0.3047442, 1e-7);
        assert_approx_eq!(out.im, 0.2082189, 1e-7);
    }
}

