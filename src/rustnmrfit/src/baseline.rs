use ndarray::{prelude::*, ArcArray1, stack};

//==============================================================================
// General 1D baseline

pub struct Baseline1D {

    // Number of parameters per domain
    np: usize,

    // Basis
    basis: Array2<f64>,

    // Baseline function
    f: fn (&mut Baseline1D, &Array<f64, Ix1>),

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
            basis = Baseline1D::gen_basis(x, knots, np);

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
    pub fn eval(&mut self, p: &Array<f64, Ix1>) {       

        (self.f)(self, p);

    }

    //--------------------------------------
    pub fn eval_0(&mut self, _p: &Array<f64, Ix1>) {       

        // Nothing to do -- if there are no baseline terms, then y is unchanged

    }

    //--------------------------------------
    pub fn eval_n(&mut self, p: &Array<f64, Ix1>) {       

        // Assume that p is divided in half -- real parameters then imaginary ones,
        // so loop over the two halves
        for i in 0 .. 2 {
            let p_slice = p.slice(s![(self.np*i) .. (self.np*(i+1))]);
            let mut y_slice = self.y.slice_mut(s![i, ..]);
            y_slice.assign( &(self.basis.dot(&p_slice)) );
        }

    }

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

}
