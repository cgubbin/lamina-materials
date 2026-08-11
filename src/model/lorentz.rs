//! Lorentz oscillator model.

use crate::{
    MaterialModelError,
    model::{DrudeLorentz, delegate::delegate_analytical_material},
};
use lamina_core::{ComplexScalar, DerivativeOrder};
use num_traits::Float;
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
    pub fn new(strength: R, resonance: R, damping: R) -> Result<Self, MaterialModelError<R>> {
        Self::from_strength(strength, resonance, damping)
    }

    /// Construct from primitive oscillator strength.
    ///
    /// `strength` has units of `cm⁻²`. `resonance` and `damping` are angular
    /// wavenumbers in inverse centimetres (`cm⁻¹`).
    pub fn from_strength(
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

    /// Construct from the dimensionless zero-frequency contribution `Δε`.
    ///
    /// The primitive strength is obtained as `S = Δε Ω²`.
    /// `resonance` and `damping` are expressed in inverse centimetres (`cm⁻¹`).
    pub fn from_delta_epsilon(
        delta_epsilon: R,
        resonance: R,
        damping: R,
    ) -> Result<Self, MaterialModelError<R>> {
        validate_nonnegative("delta_epsilon", delta_epsilon)?;
        validate_positive("resonance", resonance)?;

        Self::from_strength(delta_epsilon * resonance * resonance, resonance, damping)
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
