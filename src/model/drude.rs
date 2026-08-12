//! Drude free-carrier model.

use crate::{
    MaterialModelError,
    model::{DrudeLorentz, delegate::delegate_analytical_material},
};

use lamina_units::SpectralQuantity;
use num_traits::{Float, FromPrimitive};
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
    /// Construct a Drude model from semantic spectral parameters.
    ///
    /// `epsilon_infinity` is the dimensionless high-frequency relative
    /// permittivity. `plasma_wavenumber` and `damping` may use any spectral
    /// representation implementing [`SpectralQuantity`]; both are converted
    /// internally to vacuum angular wavenumber in inverse centimetres (`cm⁻¹`).
    pub fn new<P, D>(
        epsilon_infinity: R,
        plasma_wavenumber: P,
        damping: D,
    ) -> Result<Self, MaterialModelError<R>>
    where
        R: FromPrimitive,
        P: SpectralQuantity<R>,
        D: SpectralQuantity<R>,
    {
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

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use lamina_core::{DerivativeOrder, DifferentiableMeromorphicMaterial, MeromorphicMaterial};
    use lamina_units::{Frequency, FrequencyUnit, VacuumAngularWavenumber, VacuumWavenumber};
    use ndarray::arr0;
    use num_complex::Complex64;

    const TOLERANCE: f64 = 1.0e-12;

    fn assert_complex_close(actual: Complex64, expected: Complex64) {
        let error = (actual - expected).norm();
        let scale = expected.norm().max(1.0);

        assert!(
            error <= TOLERANCE * scale,
            "actual = {actual:?}, expected = {expected:?}, relative error = {}",
            error / scale
        );
    }

    #[test]
    fn constructor_canonicalises_angular_wavenumbers() {
        let material = Drude::new(
            2.0,
            VacuumAngularWavenumber::inverse_centimetres(12_000.0),
            VacuumAngularWavenumber::inverse_centimetres(150.0),
        )
        .unwrap();

        assert_eq!(material.epsilon_infinity(), 2.0);
        assert_eq!(material.plasma_wavenumber(), 12_000.0);
        assert_eq!(material.damping(), 150.0);
    }

    #[test]
    fn equivalent_spectral_representations_construct_same_model() {
        let plasma_wavenumber = 12_000.0;
        let damping_wavenumber = 150.0;

        let canonical = Drude::new(
            2.0,
            VacuumAngularWavenumber::inverse_centimetres(plasma_wavenumber),
            VacuumAngularWavenumber::inverse_centimetres(damping_wavenumber),
        )
        .unwrap();

        let spectroscopic = Drude::new(
            2.0,
            VacuumWavenumber::inverse_centimetres(plasma_wavenumber / std::f64::consts::TAU),
            VacuumWavenumber::inverse_centimetres(damping_wavenumber / std::f64::consts::TAU),
        )
        .unwrap();

        assert_relative_eq!(
            canonical.plasma_wavenumber(),
            spectroscopic.plasma_wavenumber(),
            epsilon = 1.0e-12,
        );

        assert_relative_eq!(
            canonical.damping(),
            spectroscopic.damping(),
            epsilon = 1.0e-12,
        );
    }

    #[test]
    fn frequency_input_is_canonicalised() {
        let frequency_hz = 2.0e14;

        let material = Drude::new(
            1.0,
            Frequency::new(frequency_hz, FrequencyUnit::Hertz),
            VacuumAngularWavenumber::inverse_centimetres(100.0),
        )
        .unwrap();

        let expected = std::f64::consts::TAU * frequency_hz / 2.997_924_58e10;

        assert_relative_eq!(material.plasma_wavenumber(), expected, epsilon = 1.0e-10,);
    }

    #[test]
    fn negative_canonicalised_plasma_wavenumber_is_rejected() {
        let result = Drude::new(
            1.0,
            VacuumAngularWavenumber::inverse_centimetres(-1.0),
            VacuumAngularWavenumber::inverse_centimetres(0.1),
        );

        assert!(matches!(
            result,
            Err(MaterialModelError::NegativeParameter {
                name: "plasma_wavenumber",
                ..
            })
        ));
    }

    #[test]
    fn negative_canonicalised_damping_is_rejected() {
        let result = Drude::new(
            1.0,
            VacuumAngularWavenumber::inverse_centimetres(1.0),
            VacuumAngularWavenumber::inverse_centimetres(-0.1),
        );

        assert!(matches!(
            result,
            Err(MaterialModelError::NegativeParameter {
                name: "drude_damping",
                ..
            })
        ));
    }

    #[test]
    fn value_matches_closed_form() {
        let epsilon_infinity = 1.4;
        let plasma = 2.1;
        let damping = 0.4;

        let material = Drude::new(
            epsilon_infinity,
            VacuumAngularWavenumber::inverse_centimetres(plasma),
            VacuumAngularWavenumber::inverse_centimetres(damping),
        )
        .unwrap();

        let k0 = Complex64::new(1.3, 0.1);

        let denominator = k0 * k0 + Complex64::i() * damping * k0;

        let expected = epsilon_infinity - plasma * plasma / denominator;

        let actual = material.relative_permittivity_complex::<_, Complex64>(arr0(k0))[()];

        assert_complex_close(actual, expected);
    }

    #[test]
    fn derivatives_match_closed_form() {
        let plasma = 2.1;
        let damping = 0.4;

        let material = Drude::new(
            1.0,
            VacuumAngularWavenumber::inverse_centimetres(plasma),
            VacuumAngularWavenumber::inverse_centimetres(damping),
        )
        .unwrap();

        let k0 = Complex64::new(1.3, 0.1);

        let strength = plasma * plasma;
        let igamma = Complex64::i() * damping;
        let d = k0 * k0 + igamma * k0;
        let d1 = 2.0 * k0 + igamma;
        let d2 = Complex64::new(2.0, 0.0);

        let expected_first = strength * d1 / (d * d);

        let expected_second = strength * (d2 * d - 2.0 * d1 * d1) / (d * d * d);

        let expected_third = 6.0 * strength * d1 * (d1 * d1 - d2 * d) / (d * d * d * d);

        for (order, expected) in [
            (DerivativeOrder::First, expected_first),
            (DerivativeOrder::Second, expected_second),
            (DerivativeOrder::Third, expected_third),
        ] {
            let actual = material
                .relative_permittivity_complex_derivative::<_, Complex64>(arr0(k0), order)[()];

            assert_complex_close(actual, expected);
        }
    }
}
