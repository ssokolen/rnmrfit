pub mod peak;
pub mod common;
pub mod constraint;
pub mod lineshape;
pub mod baseline;
pub mod phase;
pub mod fit;

pub use fit::Fit1D;
//pub use baseline::gen_basis;
pub use peak::{f_peak, f_peak_area};
