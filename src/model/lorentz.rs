//! Lorentz oscillator model.

use crate::{
    MaterialModelError,
    model::{DrudeLorentz, delegate::delegate_analytical_material},
};
use lamina_core::{ComplexScalar, DerivativeOrder};
use lamina_units::{SpectralQuantity, VacuumAngularWavenumberSquared};
use num_traits::{Float, FromPrimitive};
use std::fmt::Debug;

/// One Lorentz oscillator contribution to relative permittivity.
///
/// The contribution is
///
/// `Δε(k₀) = S / (Ω² - k₀² - i Γ k₀)`,
///
/// where `k₀`, the resonance angular wavenumber `Ω`, and damping angular
/// wavenumber `Γ` are expressed in inverse centimetres (`cm⁻¹`).
///
/// The primitive oscillator strength `S` therefore has units of `cm⁻²`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LorentzOscillator<R> {
    strength: R,
    resonance: R,
    damping: R,
}

impl<R> LorentzOscillator<R>
where
    R: Float + Debug,
{
    /// Construct from primitive oscillator strength.
    ///
    /// `strength` has units of `cm⁻²`. `resonance` and `damping` are angular
    /// wavenumbers in inverse centimetres (`cm⁻¹`).
    pub(crate) fn new(
        strength: R,
        resonance: R,
        damping: R,
    ) -> Result<Self, MaterialModelError<R>> {
        validate_nonnegative("strength", strength)?;
        validate_positive("resonance", resonance)?;
        validate_nonnegative("damping", damping)?;
        Ok(Self {
            strength,
            resonance,
            damping,
        })
    }

    /// Construct from primitive oscillator strength.
    ///
    /// `strength` must be expressed as [`VacuumAngularWavenumberSquared`] through the lamina-units
    /// typing system. `resonance` and `damping` may use any
    /// spectral representation implementing [`SpectralQuantity`].
    ///
    /// All three parameters are converted immediately to the canonical
    /// inverse-centimetre representation before validation and storage.
    pub fn from_strength<Q, D>(
        strength: VacuumAngularWavenumberSquared<R>,
        resonance: Q,
        damping: D,
    ) -> Result<Self, MaterialModelError<R>>
    where
        R: FromPrimitive,
        Q: SpectralQuantity<R>,
        D: SpectralQuantity<R>,
    {
        Self::new(
            strength.as_inverse_centimetres_squared(),
            resonance
                .vacuum_angular_wavenumber()
                .as_inverse_centimetres(),
            damping.vacuum_angular_wavenumber().as_inverse_centimetres(),
        )
    }

    /// Construct from the dimensionless zero-frequency contribution `Δε`.
    ///
    /// The primitive strength is obtained as `S = Δε Ω²`.
    /// `resonance` and `damping` may use any spectral representation implementing
    /// [`SpectralQuantity`] and are converted immediately to canonical vacuum
    /// angular wavenumber in inverse centimetres (`cm⁻¹`).
    pub fn from_delta_epsilon<Q, D>(
        delta_epsilon: R,
        resonance: Q,
        damping: D,
    ) -> Result<Self, MaterialModelError<R>>
    where
        R: FromPrimitive,
        Q: SpectralQuantity<R>,
        D: SpectralQuantity<R>,
    {
        validate_nonnegative("delta_epsilon", delta_epsilon)?;

        let resonance = resonance
            .vacuum_angular_wavenumber()
            .as_inverse_centimetres();

        validate_nonnegative("resonance", resonance)?;

        Self::new(
            delta_epsilon * resonance * resonance,
            resonance,
            damping.vacuum_angular_wavenumber().as_inverse_centimetres(),
        )
    }

    /// Return the primitive oscillator strength in `cm⁻²`.
    pub fn strength(&self) -> R {
        self.strength
    }

    /// Return the dimensionless zero-frequency contribution `S / Ω²`.
    pub fn delta_epsilon(&self) -> R {
        self.strength / (self.resonance * self.resonance)
    }

    /// Return the resonance angular wavenumber in inverse centimetres (`cm⁻¹`).
    pub fn resonance(&self) -> R {
        self.resonance
    }

    /// Return the damping angular wavenumber in inverse centimetres (`cm⁻¹`).
    pub fn damping(&self) -> R {
        self.damping
    }

    pub(crate) fn value_at<C>(&self, k0: C) -> C
    where
        C: ComplexScalar<RealField = R> + Copy,
    {
        let denominator = C::from_real(self.resonance * self.resonance)
            - k0 * k0
            - C::i() * C::from_real(self.damping) * k0;

        C::from_real(self.strength) / denominator
    }

    pub(super) fn derivative_at<C>(&self, k0: C, order: DerivativeOrder) -> C
    where
        C: ComplexScalar<RealField = R> + Copy,
    {
        let two = C::from_real(R::one() + R::one());
        let igamma = C::i() * C::from_real(self.damping);

        let d = C::from_real(self.resonance * self.resonance) - k0 * k0 - igamma * k0;

        let d1 = -two * k0 - igamma;
        let strength = C::from_real(self.strength);

        let d2 = d * d;

        match order {
            DerivativeOrder::First => -strength * d1 / d2,

            DerivativeOrder::Second => {
                let d3 = d2 * d;

                strength * (two * d1 * d1 + two * d) / d3
            }

            DerivativeOrder::Third => {
                let six = two + two + two;
                let d3 = d2 * d;
                let d4 = d3 * d;

                six * strength * d1 * (-two * d - d1 * d1) / d4
            }
        }
    }
}

/// Pure Lorentz relative-permittivity model.
///
/// The model consists of a dimensionless high-frequency permittivity plus one
/// or more [`LorentzOscillator`] contributions. All spectral oscillator
/// parameters use the canonical inverse-centimetre convention.
#[derive(Clone, Debug, PartialEq)]
pub struct Lorentz<R> {
    pub(crate) inner: DrudeLorentz<R>,
}

impl<R> Lorentz<R>
where
    R: Float + Debug,
{
    /// Construct a Lorentz model.
    pub fn new(
        epsilon_infinity: R,
        oscillators: Vec<LorentzOscillator<R>>,
    ) -> Result<Self, MaterialModelError<R>> {
        Ok(Self {
            inner: DrudeLorentz::lorentz_only(epsilon_infinity, oscillators)?,
        })
    }

    /// Return the high-frequency permittivity.
    pub fn epsilon_infinity(&self) -> R {
        self.inner.epsilon_infinity()
    }

    /// Return oscillator terms.
    pub fn oscillators(&self) -> &[LorentzOscillator<R>] {
        self.inner.oscillators()
    }
}

delegate_analytical_material!(Lorentz);

fn validate_finite<R>(name: &'static str, value: R) -> Result<(), MaterialModelError<R>>
where
    R: Float + Debug,
{
    if !value.is_finite() {
        return Err(MaterialModelError::NonFiniteParameter { name, value });
    }
    Ok(())
}

fn validate_nonnegative<R>(name: &'static str, value: R) -> Result<(), MaterialModelError<R>>
where
    R: Float + Debug,
{
    validate_finite(name, value)?;
    if value < R::zero() {
        return Err(MaterialModelError::NegativeParameter { name, value });
    }
    Ok(())
}

fn validate_positive<R>(name: &'static str, value: R) -> Result<(), MaterialModelError<R>>
where
    R: Float + Debug,
{
    validate_finite(name, value)?;
    if value <= R::zero() {
        return Err(MaterialModelError::NonPositiveParameter { name, value });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    use approx::assert_relative_eq;
    use lamina_units::{VacuumAngularWavenumber, VacuumAngularWavenumberSquared, VacuumWavenumber};
    use num_complex::Complex64;

    const TOLERANCE: f64 = 1e-12;

    fn assert_complex_close(actual: Complex64, expected: Complex64) {
        let error = (actual - expected).norm();
        let scale = expected.norm().max(1.0);

        assert!(
            error <= TOLERANCE * scale,
            "actual = {actual:?}, expected = {expected:?}, relative error = {}",
            error / scale
        );
    }

    fn expected_value(strength: f64, resonance: f64, damping: f64, k0: Complex64) -> Complex64 {
        let denominator = resonance * resonance - k0 * k0 - Complex64::i() * damping * k0;

        strength / denominator
    }

    fn expected_first(strength: f64, resonance: f64, damping: f64, k0: Complex64) -> Complex64 {
        let igamma = Complex64::i() * damping;
        let d = resonance * resonance - k0 * k0 - igamma * k0;
        let d1 = -2.0 * k0 - igamma;

        -strength * d1 / (d * d)
    }

    fn expected_second(strength: f64, resonance: f64, damping: f64, k0: Complex64) -> Complex64 {
        let igamma = Complex64::i() * damping;
        let d = resonance * resonance - k0 * k0 - igamma * k0;
        let d1 = -2.0 * k0 - igamma;
        let d2 = -2.0;

        strength * (2.0 * d1 * d1 - d2 * d) / (d * d * d)
    }

    fn expected_third(strength: f64, resonance: f64, damping: f64, k0: Complex64) -> Complex64 {
        let igamma = Complex64::i() * damping;
        let d = resonance * resonance - k0 * k0 - igamma * k0;
        let d1 = -2.0 * k0 - igamma;
        let d2 = -2.0;

        6.0 * strength * d1 * (d2 * d - d1 * d1) / (d * d * d * d)
    }

    #[test]
    fn from_strength_canonicalises_inputs() {
        let oscillator = LorentzOscillator::from_strength(
            VacuumAngularWavenumberSquared::inverse_centimetres_squared(400.0),
            VacuumAngularWavenumber::inverse_centimetres(20.0),
            VacuumAngularWavenumber::inverse_centimetres(2.0),
        )
        .unwrap();

        assert_eq!(oscillator.strength(), 400.0);
        assert_eq!(oscillator.resonance(), 20.0);
        assert_eq!(oscillator.damping(), 2.0);
    }

    #[test]
    fn from_strength_accepts_equivalent_spectral_representations() {
        let canonical = LorentzOscillator::from_strength(
            VacuumAngularWavenumberSquared::inverse_centimetres_squared(400.0),
            VacuumAngularWavenumber::inverse_centimetres(20.0),
            VacuumAngularWavenumber::inverse_centimetres(2.0),
        )
        .unwrap();

        let spectroscopic = LorentzOscillator::from_strength(
            VacuumAngularWavenumberSquared::inverse_centimetres_squared(400.0),
            VacuumWavenumber::inverse_centimetres(20.0 / std::f64::consts::TAU),
            VacuumWavenumber::inverse_centimetres(2.0 / std::f64::consts::TAU),
        )
        .unwrap();

        assert_relative_eq!(
            canonical.resonance(),
            spectroscopic.resonance(),
            epsilon = 1.0e-12,
        );

        assert_relative_eq!(
            canonical.damping(),
            spectroscopic.damping(),
            epsilon = 1.0e-12,
        );

        assert_relative_eq!(
            canonical.strength(),
            spectroscopic.strength(),
            epsilon = 1.0e-12,
        );
    }

    #[test]
    fn from_delta_epsilon_constructs_expected_strength() {
        let oscillator = LorentzOscillator::from_delta_epsilon(
            2.0,
            VacuumAngularWavenumber::inverse_centimetres(3.0),
            VacuumAngularWavenumber::inverse_centimetres(0.1),
        )
        .unwrap();

        assert_eq!(oscillator.strength(), 18.0);
        assert_eq!(oscillator.delta_epsilon(), 2.0);
        assert_eq!(oscillator.resonance(), 3.0);
        assert_eq!(oscillator.damping(), 0.1);
    }

    #[test]
    fn from_delta_epsilon_canonicalises_spectral_inputs_before_forming_strength() {
        let oscillator = LorentzOscillator::from_delta_epsilon(
            2.0,
            VacuumWavenumber::inverse_centimetres(3.0 / std::f64::consts::TAU),
            VacuumWavenumber::inverse_centimetres(0.1 / std::f64::consts::TAU),
        )
        .unwrap();

        assert_relative_eq!(oscillator.resonance(), 3.0, epsilon = 1.0e-12,);

        assert_relative_eq!(oscillator.damping(), 0.1, epsilon = 1.0e-12,);

        assert_relative_eq!(oscillator.strength(), 18.0, epsilon = 1.0e-11,);
    }

    #[test]
    fn negative_strength_is_rejected_after_canonicalisation() {
        let result = LorentzOscillator::from_strength(
            VacuumAngularWavenumberSquared::inverse_centimetres_squared(-1.0),
            VacuumAngularWavenumber::inverse_centimetres(2.0),
            VacuumAngularWavenumber::inverse_centimetres(0.1),
        );

        assert!(matches!(
            result,
            Err(MaterialModelError::NegativeParameter {
                name: "strength",
                ..
            })
        ));
    }

    #[test]
    fn zero_resonance_is_rejected_after_canonicalisation() {
        let result = LorentzOscillator::from_delta_epsilon(
            1.0,
            VacuumAngularWavenumber::inverse_centimetres(0.0),
            VacuumAngularWavenumber::inverse_centimetres(0.1),
        );

        assert!(matches!(
            result,
            Err(MaterialModelError::NonPositiveParameter {
                name: "resonance",
                ..
            })
        ));
    }

    #[test]
    fn negative_damping_is_rejected_after_canonicalisation() {
        let result = LorentzOscillator::from_delta_epsilon(
            1.0,
            VacuumAngularWavenumber::inverse_centimetres(2.0),
            VacuumAngularWavenumber::inverse_centimetres(-0.1),
        );

        assert!(matches!(
            result,
            Err(MaterialModelError::NegativeParameter {
                name: "damping",
                ..
            })
        ));
    }

    #[test]
    fn negative_delta_epsilon_is_rejected() {
        let result = LorentzOscillator::from_delta_epsilon(
            -1.0,
            VacuumAngularWavenumber::inverse_centimetres(2.0),
            VacuumAngularWavenumber::inverse_centimetres(0.1),
        );

        assert!(matches!(
            result,
            Err(MaterialModelError::NegativeParameter {
                name: "delta_epsilon",
                ..
            })
        ));
    }

    #[test]
    fn value_matches_closed_form() {
        let oscillator = LorentzOscillator::new(1.7, 2.6, 0.35).unwrap();

        let k0 = Complex64::new(0.9, 0.15);

        assert_complex_close(oscillator.value_at(k0), expected_value(1.7, 2.6, 0.35, k0));
    }

    #[test]
    fn first_derivative_matches_closed_form() {
        let oscillator = LorentzOscillator::new(1.7, 2.6, 0.35).unwrap();

        let k0 = Complex64::new(0.9, 0.15);

        assert_complex_close(
            oscillator.derivative_at(k0, DerivativeOrder::First),
            expected_first(1.7, 2.6, 0.35, k0),
        );
    }

    #[test]
    fn second_derivative_matches_closed_form() {
        let oscillator = LorentzOscillator::new(1.7, 2.6, 0.35).unwrap();

        let k0 = Complex64::new(0.9, 0.15);

        assert_complex_close(
            oscillator.derivative_at(k0, DerivativeOrder::Second),
            expected_second(1.7, 2.6, 0.35, k0),
        );
    }

    #[test]
    fn third_derivative_matches_closed_form() {
        let oscillator = LorentzOscillator::new(1.7, 2.6, 0.35).unwrap();

        let k0 = Complex64::new(0.9, 0.15);

        assert_complex_close(
            oscillator.derivative_at(k0, DerivativeOrder::Third),
            expected_third(1.7, 2.6, 0.35, k0),
        );
    }
}
