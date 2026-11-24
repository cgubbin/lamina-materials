//! Optical material models for the `lamina` ecosystem.
//!
//! Analytical models implement the material capability traits from `lamina-core`.
//! Tabulated optical data implements only real-axis evaluation because a
//! real-axis interpolation does not define a physically meaningful complex
//! continuation.
//!
//! # Spectral convention
//!
//! Material evaluation uses vacuum angular wavenumber
//!
//! `k₀ = ω / c = 2π / λ₀`
//!
//! in inverse centimetres (`cm⁻¹`) as its canonical spectral coordinate.
//! Public material constructors accept semantic spectral quantities from
//! `lamina-units` and convert them to this canonical representation before
//! validation and storage.

pub mod database;
pub mod error;
pub mod metadata;
pub mod model;

pub use error::{MaterialModelError, TabulatedMaterialError};
pub use metadata::{MaterialRecord, SpectralDomain};
pub use model::{
    Cauchy, Drude, DrudeLorentz, Extrapolation, Lorentz, LorentzOscillator, Sellmeier,
    SellmeierParameterRates, SellmeierTerm, Tabulated,
};
