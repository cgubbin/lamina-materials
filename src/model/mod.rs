//! Generic analytical optical material models.
//!
//! Types in this module represent mathematical constitutive models. They do not
//! carry database identity, literature references, or spectral validity ranges.
//! Named parameter sets belong in [`crate::database`].

mod cauchy;
mod delegate;
mod drude;
mod drude_lorentz;
mod lorentz;
mod rational;
mod sellmeir;
mod tabulated;
mod tensor;

pub use cauchy::Cauchy;
pub use drude::Drude;
pub use drude_lorentz::DrudeLorentz;
pub use lorentz::{Lorentz, LorentzOscillator};
pub use sellmeir::{Sellmeier, SellmeierParameterRates, SellmeierTerm};
pub use tabulated::{Extrapolation, Tabulated};
pub use tensor::{
    ConstantTensor, CrystalAxis, PrincipalAxes, PrincipalModel, SellmeierAxes, TensorParameter,
    TensorRate,
};
