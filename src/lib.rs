pub mod mat;
pub mod mesh;
pub mod vct;

pub use mat::Mat;
pub use mesh::Mesh;
pub use vct::Vct;
pub type Flt = f64;

pub const EPS: Flt = 1e-5;
pub const PI: Flt = std::f64::consts::PI as Flt;
