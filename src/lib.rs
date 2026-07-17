//! Optical material models for the `tmm` ecosystem.
//!
//! Analytical models implement the material capability traits from `tmm-core`.
//! Tabulated optical data implements only real-axis evaluation because a
//! real-axis interpolation does not define a physically meaningful complex
//! continuation.
//!
//! The spectral coordinate is vacuum angular wavenumber,
//! `k₀ = ω / c = 2π / λ₀`. Model parameters must use consistent units.

pub mod database;
pub mod error;
pub mod metadata;
pub mod model;
pub mod units;

pub use error::{MaterialModelError, TabulatedMaterialError};
pub use metadata::{MaterialRecord, SpectralDomain};
pub use model::{
    Cauchy, Drude, DrudeLorentz, Extrapolation, Lorentz, LorentzOscillator, Sellmeier,
    SellmeierTerm, Tabulated,
};
pub use units::{LengthUnit, Wavelength};
