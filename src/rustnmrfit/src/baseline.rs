use ndarray::{prelude::*, stack};

//==============================================================================

pub fn gen_basis(x: &Array1<f64>, knots: &Array1<f64>, np: usize) -> Array2<f64> { 

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
        for j in 0 .. n_basis {
            basis[[i,j]] = beta(x[i], j, degree, &knots);
        }
    }

    // Fudging last value by assuming boundary knots correspond to end points
    basis[[x.len() - 1, n_basis - 1]] = 1.0;

    basis
}

//==============================================================================
// General 1D baseline

pub struct Baseline1D {

    // Basis
    basis: Array2<f64>,

    // Baseline function
    f: fn (&mut Baseline1D, &Array<f64, Ix1>, &Array<f64, Ix1>),

    // Final output
    pub y: Array2<f64>,
    pub dydp: Array3<f64>,

}

//--------------------------------------
impl Baseline1D {

    //--------------------------------------
    pub fn new(x: &Array1<f64>, knots: &Array1<f64>, np: usize) -> Baseline1D {

        let n: usize = x.len();

        let basis: Array2<f64> = Array::zeros((n, np));
        let dydp: Array3<f64> = Array::zeros((2, np, n));
        let f = Baseline1D::eval_0;

        if np > 0 {
            basis = gen_basis(x, knots, np);

            let mut dydp_r = dydp.slice_mut(s![0, .., ..]);
            dydp_r.assign(&basis.t());

            let mut dydp_i = dydp.slice_mut(s![1, .., ..]);
            dydp_i.assign(&basis.t());

            f = Baseline1D::eval_n;
        }

        Baseline1D {
            basis: basis,

            f: f,

            y: Array::zeros((2, n)),
            dydp: dydp,
        }
    }

    //--------------------------------------
    pub fn eval(&mut self, x: &Array<f64, Ix1>, p: &Array<f64, Ix1>) {       

        (self.f)(self, x, p);

    }

    //--------------------------------------
    pub fn eval_0(&mut self, _x: &Array<f64, Ix1>, _p: &Array<f64, Ix1>) {       

        // Nothing to do -- if there are no baseline terms, then y is unchanged

    }

    //--------------------------------------
    pub fn eval_n(&mut self, _x: &Array<f64, Ix1>, p: &Array<f64, Ix1>) {       

        //... TODO: Figure out whether input should be 1 vs 2 parameters (and how 2D will compare)
        // math is just a dot product of basis and p terms


    }

}
