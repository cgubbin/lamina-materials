//! Drude free-carrier model.

use crate::{
    MaterialModelError,
    model::{DrudeLorentz, delegate::delegate_analytical_material},
};
use num_traits::Float;
use std::fmt::Debug;

/// Drude free-carrier relative-permittivity model.
///
/// The model is
///
/// `ε(k₀) = ε∞ - Ωᴅ² / (k₀² + i Γᴅ k₀)`,
///
/// where `k₀`, the plasma angular wavenumber `Ωᴅ`, and the damping angular
/// wavenumber `Γᴅ` are expressed in inverse centimetres (`cm⁻¹`).
///
/// `ε∞` is dimensionless.
#[derive(Clone, Debug, PartialEq)]
pub struct Drude<R> {
    pub(crate) inner: DrudeLorentz<R>,
}

impl<R> Drude<R>
where
    R: Float + Debug,
{
    /// Construct a Drude model.
    ///
    /// `epsilon_infinity` is the dimensionless high-frequency relative
    /// permittivity. `plasma_wavenumber` and `damping` are expressed in inverse
    /// centimetres (`cm⁻¹`).
    pub fn new(
        epsilon_infinity: R,
        plasma_wavenumber: R,
        damping: R,
    ) -> Result<Self, MaterialModelError<R>> {
        Ok(Self {
            inner: DrudeLorentz::from_plasma_wavenumber(
                epsilon_infinity,
                plasma_wavenumber,
                damping,
                Vec::new(),
            )?,
        })
    }

    /// Return the dimensionless high-frequency relative permittivity.
    pub fn epsilon_infinity(&self) -> R {
        self.inner.epsilon_infinity()
    }

    /// Return the plasma angular wavenumber in inverse centimetres (`cm⁻¹`).
    pub fn plasma_wavenumber(&self) -> R {
        self.inner.plasma_wavenumber()
    }

    /// Return the damping angular wavenumber in inverse centimetres (`cm⁻¹`).
    pub fn damping(&self) -> R {
        self.inner.drude_damping()
    }
}

delegate_analytical_material!(Drude);
