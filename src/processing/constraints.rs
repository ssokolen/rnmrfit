use ndarray::prelude::*;
use serde::{Deserialize};
use std::iter;

//==============================================================================
#[derive(Clone, Debug, Deserialize)]
#[serde(default)]
pub struct ConstraintLeeways {

    #[serde(default)]
    pub position: f64,
    #[serde(default)]
    pub width: f64,
    #[serde(default)]
    pub height: f64,
    #[serde(default)]
    pub fraction: f64,
    #[serde(default)]
    pub area: f64,

}


impl Default for ConstraintLeeways {

    fn default() -> Self {
        ConstraintLeeways {
            position: 0.0,
            width: 0.0,
            height: 0.0,
            fraction: 0.0,
            area: 0.0,
        }
    }
}


//==============================================================================
#[derive(Clone, Debug)]
pub struct PositionConstraint {
    indexes: (usize, usize),
    value: f64,
    leeway: f64
}

#[derive(Clone, Debug)]
pub struct WidthConstraint {
    indexes: (usize, usize),
    leeway: f64
}

#[derive(Clone, Debug)]
pub struct HeightConstraint {
    indexes: (usize, usize),
    value: f64,
    leeway: f64
}

#[derive(Clone, Debug)]
pub struct FractionConstraint {
    indexes: (usize, usize),
    leeway: f64
}

#[derive(Clone, Debug)]
pub struct AreaConstraint {
    _indexes: (Vec<usize>, Vec<usize>),
    _value: f64,
    leeway: f64
}


//==============================================================================
#[derive(Clone, Debug)]
pub enum Constraint {
    Position(PositionConstraint),
    Width(WidthConstraint),
    Height(HeightConstraint),
    Fraction(FractionConstraint),
    Area(AreaConstraint)
}


impl Constraint {

    //--------------------------------------------------------------------------
    pub fn position(indexes: (usize, usize), value: f64, leeway: f64) -> Self {
        Constraint::Position( PositionConstraint {
            indexes: (indexes.0 * 4, indexes.1 * 4),
            value: value,
            leeway: leeway,
        })
    }

    pub fn width(indexes: (usize, usize), leeway: f64) -> Self {
        Constraint::Width( WidthConstraint {
            indexes: (indexes.0 * 4 + 1, indexes.1 * 4 + 1),
            leeway: leeway,
        })
    }

    pub fn height(indexes: (usize, usize), value: f64, leeway: f64) -> Self {
        Constraint::Height( HeightConstraint {
            indexes: (indexes.0 * 4 + 2, indexes.1 * 4 + 2),
            value: value, 
            leeway: leeway,
        })
    }

    pub fn fraction(indexes: (usize, usize), leeway: f64) -> Self {
        Constraint::Fraction( FractionConstraint {
            indexes: (indexes.0 * 4 + 3, indexes.1 * 4 + 3),
            leeway: leeway,
        })
    }

    pub fn area(
        indexes: (Vec<usize>, Vec<usize>), 
        value: f64, 
        leeway: f64
    ) -> Self {
        Constraint::Area( AreaConstraint {
            _indexes: indexes, 
            _value: value, 
            leeway: leeway,
        })
    }


    //--------------------------------------------------------------------------
    pub fn update(&mut self, leeways: &ConstraintLeeways) {
        
        let (old_leeway, new_leeway) = match self {
            Constraint::Position(object) => {
                (&mut object.leeway, leeways.position)
            },
            Constraint::Height(object) => {
                (&mut object.leeway, leeways.height)
            },
            Constraint::Width(object) => {
                (&mut object.leeway, leeways.width)
            },
            Constraint::Fraction(object) => {
                (&mut object.leeway, leeways.fraction)
            },
            Constraint::Area(object) => {
                (&mut object.leeway, leeways.area)
            }
        };

        // Update leeways only if bigger
        if new_leeway > *old_leeway {
            *old_leeway = new_leeway;
        }
    }
}


//==============================================================================
#[derive(Clone, Debug)]
pub struct PeakEquality {
    mat: Array2<f64>,
    vec: Array1<f64>,
    x: Array1<f64>,
    gradient: Vec<f64>
}


impl PeakEquality {

    //--------------------------------------------------------------------------
    pub fn new(
        constraints: &Vec<Constraint>, 
        n_parameters: usize,
        xcoeff: f64
    ) -> PeakEquality {

        // Extract equality constraints
        let constraints: Vec<_> = constraints.iter()
            .filter(|x| {
                let leeway = match x {
                    Constraint::Position(obj) => obj.leeway,
                    Constraint::Width(obj) => obj.leeway,
                    Constraint::Height(obj) => obj.leeway,
                    Constraint::Fraction(obj) => obj.leeway,
                    _ => f64::NAN
                };

                leeway == 0.0
            })
            .map(|x| x.clone())
            .collect();

        // Generate parameter matrix
        let n_row = constraints.len();

        let mut mat = Array2::<f64>::zeros((n_row, n_parameters));
        let mut vec = Array1::<f64>::zeros((n_row,));

        for (i, constraint) in constraints.into_iter().enumerate() {
            match constraint {
                Constraint::Position(obj) => {
                    mat[[i, obj.indexes.1]] = 1.0;
                    mat[[i, obj.indexes.0]] = -1.0;

                    // Scale position difference to match x
                    vec[[i]] = obj.value/xcoeff;
                }

                Constraint::Width(obj) => {
                    mat[[i, obj.indexes.1]] = 1.0;
                    mat[[i, obj.indexes.0]] = -1.0;
                    
                    // vec remains 0
                }

                Constraint::Height(obj) => {
                    mat[[i, obj.indexes.1]] = 1.0;
                    mat[[i, obj.indexes.0]] = -obj.value;

                    // vec remains 0
                }

                Constraint::Fraction(obj) => {
                    mat[[i, obj.indexes.1]] = 1.0;
                    mat[[i, obj.indexes.0]] = -1.0;
                    
                    // vec remains 0
                }

                _ => {},
            }
        }

        PeakEquality {
            mat: mat.clone(),
            vec: vec,
            x: Array1::<f64>::zeros((n_parameters,)),
            gradient: mat.into_iter().collect()
        }
    }

    //--------------------------------------------------------------------------
    pub fn m(&self) -> usize {
        self.mat.nrows()
    }


    //--------------------------------------------------------------------------
    pub fn mobj(
        result: &mut [f64], 
        x: &[f64], 
        gradient: Option<&mut [f64]>, 
        user_data: &mut PeakEquality
    ) {
        // Fill in Array objs
        for (x, y) in user_data.x.iter_mut().zip(x.iter()) {
            *x = *y;
        }
        
        let out = user_data.mat.dot(&user_data.x) - &user_data.vec;

        for (x, y) in result.iter_mut().zip(out.iter()) {
            *x = *y;        
        }

        // For parameters, the gradient is just matrix A
        if let Some(gradient) = gradient {

            for (x, y) in gradient.iter_mut().zip(user_data.gradient.iter()) {
                *x = *y;
            }
        }
    }
}



//==============================================================================
#[derive(Clone, Debug)]
pub struct PeakInequality {
    mat: Array2<f64>,
    vec: Array1<f64>,
    x: Array1<f64>,
    gradient: Vec<f64>
}


impl PeakInequality {

    //--------------------------------------------------------------------------
    pub fn new(
        constraints: &Vec<Constraint>, 
        n_parameters: usize,
        xcoeff: f64
    ) -> PeakInequality {

        // Extract equality constraints
        let constraints: Vec<_> = constraints.iter()
            .filter(|x| {
                let leeway = match x {
                    Constraint::Position(obj) => obj.leeway,
                    Constraint::Width(obj) => obj.leeway,
                    Constraint::Height(obj) => obj.leeway,
                    Constraint::Fraction(obj) => obj.leeway,
                    _ => f64::NAN
                };

                leeway != 0.0
            })
            .map(|x| x.clone())
            .collect();

        // Generate parameter matrix
        let n_con = constraints.len();
        let n_row = n_con*2;

        let mut mat = Array2::<f64>::zeros((n_row, n_parameters));
        let mut vec = Array1::<f64>::zeros((n_row,));

        for (i, constraint) in constraints.into_iter().enumerate() {
            match constraint {
                Constraint::Position(obj) => {
                    mat[[i, obj.indexes.1]] = 1.0;
                    mat[[i, obj.indexes.0]] = -1.0;

                    mat[[i+n_con, obj.indexes.1]] = -1.0;
                    mat[[i+n_con, obj.indexes.0]] = 1.0;

                    // Scale position difference to match x
                    vec[[i]] = (obj.value + obj.leeway)/xcoeff;
                    vec[[i+n_con]] = (-obj.value + obj.leeway)/xcoeff;
                }

                Constraint::Width(obj) => {
                    mat[[i, obj.indexes.1]] = 1.0;
                    mat[[i, obj.indexes.0]] = -(1.0 + obj.leeway);
                    
                    mat[[i+n_con, obj.indexes.1]] = -1.0;
                    mat[[i+n_con, obj.indexes.0]] = 1.0 - obj.leeway;
                    
                    // vec remains 0
                }

                Constraint::Height(obj) => {
                    mat[[i, obj.indexes.1]] = 1.0;
                    mat[[i, obj.indexes.0]] = -(obj.value + obj.leeway);

                    mat[[i+n_con, obj.indexes.1]] = -1.0;
                    mat[[i+n_con, obj.indexes.0]] = obj.value - obj.leeway;

                    // vec remains 0
                }

                Constraint::Fraction(obj) => {
                    mat[[i, obj.indexes.1]] = 1.0;
                    mat[[i, obj.indexes.0]] = -1.0;

                    mat[[i, obj.indexes.1]] = -1.0;
                    mat[[i, obj.indexes.0]] = 1.0;
                    
                    vec[[i]] = obj.leeway;
                    vec[[i+n_con]] = obj.leeway;
                }

                _ => {},
            }
        }

        PeakInequality {
            mat: mat.clone(),
            vec: vec,
            x: Array1::<f64>::zeros((n_parameters,)),
            gradient: mat.into_iter().collect()
        }
    }

    //--------------------------------------------------------------------------
    pub fn m(&self) -> usize {
        self.mat.nrows()
    }


    //--------------------------------------------------------------------------
    pub fn mobj(
        result: &mut [f64], 
        x: &[f64], 
        gradient: Option<&mut [f64]>, 
        user_data: &mut PeakInequality 
    ) {
        // Fill in Array objs
        for (x, y) in user_data.x.iter_mut().zip(x.iter()) {
            *x = *y;
        }
        
        let out = user_data.mat.dot(&user_data.x) - &user_data.vec;

        for (x, y) in result.iter_mut().zip(out.iter()) {
            *x = *y;        
        }

        // For parameters, the gradient is just matrix A
        if let Some(gradient) = gradient {

            for (x, y) in gradient.iter_mut().zip(user_data.gradient.iter()) {
                *x = *y;
            }
        }
    }
}

//==============================================================================
// Only applicable for phase order of 1
#[derive(Clone, Debug)]
pub struct PhaseInequality {
    index: usize,
    constraint: f64,
    gradient: Vec<f64>
}


impl PhaseInequality {

    //--------------------------------------------------------------------------
    pub fn new(constraint: f64, n_parameters: usize) -> PhaseInequality {

        // The two inequalities are theta_0 +/- theta_1 < constraint
        // (so the gradients are 0, 0, ..., 1, 1, 0, 0, ..., 1, -1)
        // As well as theta_0 +/- theta_1 > -constraint

        PhaseInequality {
            index: n_parameters - 2,
            constraint: constraint,
            gradient: iter::repeat(0.0)
                .take(n_parameters - 2)
                .chain(vec![1.0, 1.0])
                .chain(iter::repeat(0.0).take(n_parameters - 2))
                .chain(vec![1.0, -1.0])
                .chain(iter::repeat(0.0).take(n_parameters - 2))
                .chain(vec![-1.0, -1.0])
                .chain(iter::repeat(0.0).take(n_parameters - 2))
                .chain(vec![-1.0, 1.0])
                .collect()
        }
    }

    //--------------------------------------------------------------------------
    pub fn m(&self) -> usize {
        4
    }


    //--------------------------------------------------------------------------
    pub fn mobj(
        result: &mut [f64], 
        x: &[f64], 
        gradient: Option<&mut [f64]>, 
        user_data: &mut PhaseInequality 
    ) {
        // Calculating the two limits
        let i = user_data.index;

        result[0] =  x[i] + x[i+1] - user_data.constraint; 
        result[1] =  x[i] - x[i+1] - user_data.constraint; 
        result[2] = -x[i] - x[i+1] - user_data.constraint; 
        result[3] = -x[i] + x[i+1] - user_data.constraint; 

        // The gradient is constant
        if let Some(gradient) = gradient {

            for (x, y) in gradient.iter_mut().zip(user_data.gradient.iter()) {
                *x = *y;
            }
        }
    }
}
