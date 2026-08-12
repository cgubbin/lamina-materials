//! Optical material models for the `lamina` ecosystem.
//!
//! Analytical models implement the material capability traits from `lamina-core`.
//! Tabulated optical data implements only real-axis evaluation because a
//! real-axis interpolation does not define a physically meaningful complex
//! continuation.
//!
//! # Spectral convention
//!
//! The canonical spectral coordinate is the vacuum angular wavenumber
//!
//! `k₀ = ω / c = 2π / λ₀`,
//!
//! expressed in inverse centimetres (`cm⁻¹`). Spectral material parameters are
//! expressed in units consistent with this coordinate unless explicitly stated
//! otherwise. Constructors accepting wavelength-domain parameters convert them
//! to the canonical inverse-centimetre representation.

pub mod database;
pub mod error;
pub mod metadata;
pub mod model;

pub use error::{MaterialModelError, TabulatedMaterialError};
pub use metadata::{MaterialRecord, SpectralDomain};
pub use model::{
    Cauchy, Drude, DrudeLorentz, Extrapolation, Lorentz, LorentzOscillator, Sellmeier,
    SellmeierTerm, Tabulated,
};
