//! Combined Drude–Lorentz model.

use crate::{MaterialModelError, model::LorentzOscillator};
use lamina_core::{
    ComplexScalar, DerivativeOrder, DifferentiableMaterial, DifferentiableMeromorphicMaterial,
    Material, MeromorphicMaterial, Sampled,
};
use lamina_units::{SpectralQuantity, VacuumAngularWavenumberSquared};
use num_traits::{Float, FromPrimitive};
use std::fmt::Debug;

/// Combined Drude–Lorentz relative-permittivity model.
///
/// The model is
///
/// `ε(k₀) = ε∞ - Ωᴅ²/(k₀² + iΓᴅk₀)
///          + Σ Sⱼ/(Ωⱼ² - k₀² - iΓⱼk₀)`.
///
/// All angular wavenumbers are expressed in inverse centimetres (`cm⁻¹`).
/// Consequently, the primitive Drude strength `Ωᴅ²` and Lorentz oscillator
/// strengths `Sⱼ` have units of `cm⁻²`. `ε∞` is dimensionless.
#[derive(Clone, Debug, PartialEq)]
pub struct DrudeLorentz<R> {
    epsilon_infinity: R,
    drude_strength: R,
    drude_damping: R,
    oscillators: Vec<LorentzOscillator<R>>,
}

impl<R> DrudeLorentz<R>
where
    R: Float + Debug,
{
    /// Construct from primitive Drude strength.
    ///
    /// `epsilon_infinity` is dimensionless. `drude_strength` has units of `cm⁻²`,
    /// `drude_damping` has units of `cm⁻¹`, and each oscillator must use the same
    /// canonical inverse-centimetre convention.
    pub(crate) fn new(
        epsilon_infinity: R,
        drude_strength: R,
        drude_damping: R,
        oscillators: Vec<LorentzOscillator<R>>,
    ) -> Result<Self, MaterialModelError<R>> {
        validate_finite("epsilon_infinity", epsilon_infinity)?;
        validate_nonnegative("drude_strength", drude_strength)?;
        validate_nonnegative("drude_damping", drude_damping)?;

        Ok(Self {
            epsilon_infinity,
            drude_strength,
            drude_damping,
            oscillators,
        })
    }

    /// Construct from unit-aware primitive Drude parameters.
    ///
    /// `drude_strength` is converted to canonical `cm⁻²`, while
    /// `drude_damping` is converted to canonical vacuum angular wavenumber
    /// in `cm⁻¹`.
    pub fn from_strength<D>(
        epsilon_infinity: R,
        drude_strength: VacuumAngularWavenumberSquared<R>,
        drude_damping: D,
        oscillators: Vec<LorentzOscillator<R>>,
    ) -> Result<Self, MaterialModelError<R>>
    where
        R: FromPrimitive,
        D: SpectralQuantity<R>,
    {
        Self::new(
            epsilon_infinity,
            drude_strength.as_inverse_centimetres_squared(),
            drude_damping
                .vacuum_angular_wavenumber()
                .as_inverse_centimetres(),
            oscillators,
        )
    }

    /// Construct from the Drude plasma angular wavenumber.
    ///
    /// `epsilon_infinity` is dimensionless. `plasma_wavenumber` and
    /// `drude_damping` may use any spectral representation implementing
    /// [`SpectralQuantity`]. Both are converted immediately to canonical
    /// vacuum angular wavenumber in inverse centimetres (`cm⁻¹`).
    ///
    /// The primitive Drude strength is formed as `Ωᴅ²` after canonicalisation.
    pub fn from_plasma_wavenumber<P, D>(
        epsilon_infinity: R,
        plasma_wavenumber: P,
        drude_damping: D,
        oscillators: Vec<LorentzOscillator<R>>,
    ) -> Result<Self, MaterialModelError<R>>
    where
        R: FromPrimitive,
        P: SpectralQuantity<R>,
        D: SpectralQuantity<R>,
    {
        let plasma_wavenumber = plasma_wavenumber
            .vacuum_angular_wavenumber()
            .as_inverse_centimetres();

        validate_nonnegative("plasma_wavenumber", plasma_wavenumber)?;

        Self::new(
            epsilon_infinity,
            plasma_wavenumber * plasma_wavenumber,
            drude_damping
                .vacuum_angular_wavenumber()
                .as_inverse_centimetres(),
            oscillators,
        )
    }

    /// Construct a Lorentz-only model.
    ///
    /// `epsilon_infinity` is dimensionless. Oscillator spectral parameters use
    /// the canonical inverse-centimetre convention.
    pub fn lorentz_only(
        epsilon_infinity: R,
        oscillators: Vec<LorentzOscillator<R>>,
    ) -> Result<Self, MaterialModelError<R>> {
        Self::new(epsilon_infinity, R::zero(), R::zero(), oscillators)
    }

    /// Return the dimensionless high-frequency relative permittivity.
    pub fn epsilon_infinity(&self) -> R {
        self.epsilon_infinity
    }

    /// Return the primitive Drude strength `Ωᴅ²` in `cm⁻²`.
    pub fn drude_strength(&self) -> R {
        self.drude_strength
    }

    /// Return the plasma angular wavenumber in inverse centimetres (`cm⁻¹`).
    pub fn plasma_wavenumber(&self) -> R {
        self.drude_strength.sqrt()
    }

    /// Return the Drude damping angular wavenumber in inverse centimetres (`cm⁻¹`).
    pub fn drude_damping(&self) -> R {
        self.drude_damping
    }

    /// Return the Lorentz oscillator terms.
    pub fn oscillators(&self) -> &[LorentzOscillator<R>] {
        &self.oscillators
    }

    fn relative_permittivity_at<C>(&self, k0: C) -> C
    where
        C: ComplexScalar<RealField = R> + Copy,
    {
        let mut epsilon = C::from_real(self.epsilon_infinity);

        if self.drude_strength > R::zero() {
            epsilon += self.drude_value(k0);
        }

        for oscillator in &self.oscillators {
            epsilon += oscillator.value_at(k0);
        }

        epsilon
    }

    fn derivative_at<C>(&self, k0: C, order: DerivativeOrder) -> C
    where
        C: ComplexScalar<RealField = R> + Copy,
    {
        let mut derivative = C::zero();

        if self.drude_strength > R::zero() {
            derivative += self.drude_derivative_at(k0, order);
        }

        for oscillator in &self.oscillators {
            derivative += oscillator.derivative_at(k0, order);
        }

        derivative
    }

    fn drude_value<C>(&self, k0: C) -> C
    where
        C: ComplexScalar<RealField = R> + Copy,
    {
        let denominator = k0 * k0 + C::i() * C::from_real(self.drude_damping) * k0;

        -C::from_real(self.drude_strength) / denominator
    }

    fn drude_derivative_at<C>(&self, k0: C, order: DerivativeOrder) -> C
    where
        C: ComplexScalar<RealField = R> + Copy,
    {
        let two = C::from_real(R::one() + R::one());
        let igamma = C::i() * C::from_real(self.drude_damping);
        let d = k0 * k0 + igamma * k0;
        let d1 = two * k0 + igamma;
        let strength = C::from_real(self.drude_strength);

        match order {
            DerivativeOrder::First => strength * d1 / (d * d),
            DerivativeOrder::Second => {
                let d2 = two;
                strength * (d2 * d - two * d1 * d1) / (d * d * d)
            }
            DerivativeOrder::Third => {
                let d2 = two;
                let six = two + two + two;
                six * strength * d1 * (d1 * d1 - d2 * d) / (d * d * d * d)
            }
        }
    }
}

impl<R> Material for DrudeLorentz<R>
where
    R: Float + Debug,
{
    type Real = R;

    fn relative_permeability<I, C>(&self, vacuum_angular_wavenumber: I) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = R> + Copy,
        I: Sampled<Elem = R>,
    {
        vacuum_angular_wavenumber.map(|_| C::one())
    }

    fn relative_permittivity<I, C>(&self, vacuum_angular_wavenumber: I) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = R> + Copy,
        I: Sampled<Elem = R>,
    {
        vacuum_angular_wavenumber.map(|k0| self.relative_permittivity_at::<C>(C::from_real(k0)))
    }
}

impl<R> DifferentiableMaterial for DrudeLorentz<R>
where
    R: Float + Debug,
{
    fn relative_permittivity_derivative<I, C>(
        &self,
        vacuum_angular_wavenumber: I,
        order: DerivativeOrder,
    ) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = R> + Copy,
        I: Sampled<Elem = R>,
    {
        vacuum_angular_wavenumber.map(|k0| self.derivative_at::<C>(C::from_real(k0), order))
    }
}

impl<R> MeromorphicMaterial for DrudeLorentz<R>
where
    R: Float + Debug,
{
    fn relative_permeability_complex<I, C>(&self, vacuum_angular_wavenumber: I) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = R> + Copy,
        I: Sampled<Elem = C>,
    {
        vacuum_angular_wavenumber.map(|_| C::one())
    }

    fn relative_permittivity_complex<I, C>(&self, vacuum_angular_wavenumber: I) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = R> + Copy,
        I: Sampled<Elem = C>,
    {
        vacuum_angular_wavenumber.map(|k0| self.relative_permittivity_at(k0))
    }
}

impl<R> DifferentiableMeromorphicMaterial for DrudeLorentz<R>
where
    R: Float + Debug,
{
    fn relative_permittivity_complex_derivative<I, C>(
        &self,
        vacuum_angular_wavenumber: I,
        order: DerivativeOrder,
    ) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = R> + Copy,
        I: Sampled<Elem = C>,
    {
        vacuum_angular_wavenumber.map(|k0| self.derivative_at(k0, order))
    }
}

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

    #[test]
    fn from_strength_canonicalises_parameters() {
        let model = DrudeLorentz::from_strength(
            1.5,
            VacuumAngularWavenumberSquared::inverse_centimetres_squared(400.0),
            VacuumAngularWavenumber::inverse_centimetres(2.0),
            Vec::new(),
        )
        .unwrap();

        assert_eq!(model.epsilon_infinity(), 1.5);
        assert_eq!(model.drude_strength(), 400.0);
        assert_eq!(model.plasma_wavenumber(), 20.0);
        assert_eq!(model.drude_damping(), 2.0);
    }

    #[test]
    fn from_plasma_wavenumber_canonicalises_parameters() {
        let model = DrudeLorentz::from_plasma_wavenumber(
            1.5,
            VacuumAngularWavenumber::inverse_centimetres(20.0),
            VacuumAngularWavenumber::inverse_centimetres(2.0),
            Vec::new(),
        )
        .unwrap();

        assert_eq!(model.epsilon_infinity(), 1.5);
        assert_eq!(model.drude_strength(), 400.0);
        assert_eq!(model.plasma_wavenumber(), 20.0);
        assert_eq!(model.drude_damping(), 2.0);
    }

    #[test]
    fn equivalent_spectral_representations_construct_same_model() {
        let canonical = DrudeLorentz::from_plasma_wavenumber(
            1.5,
            VacuumAngularWavenumber::inverse_centimetres(20.0),
            VacuumAngularWavenumber::inverse_centimetres(2.0),
            Vec::new(),
        )
        .unwrap();

        let spectroscopic = DrudeLorentz::from_plasma_wavenumber(
            1.5,
            VacuumWavenumber::inverse_centimetres(20.0 / std::f64::consts::TAU),
            VacuumWavenumber::inverse_centimetres(2.0 / std::f64::consts::TAU),
            Vec::new(),
        )
        .unwrap();

        assert_relative_eq!(
            canonical.drude_strength(),
            spectroscopic.drude_strength(),
            epsilon = 1.0e-12,
        );

        assert_relative_eq!(
            canonical.plasma_wavenumber(),
            spectroscopic.plasma_wavenumber(),
            epsilon = 1.0e-12,
        );

        assert_relative_eq!(
            canonical.drude_damping(),
            spectroscopic.drude_damping(),
            epsilon = 1.0e-12,
        );
    }

    #[test]
    fn negative_plasma_wavenumber_is_rejected_before_squaring() {
        let result = DrudeLorentz::from_plasma_wavenumber(
            1.0,
            VacuumAngularWavenumber::inverse_centimetres(-2.0),
            VacuumAngularWavenumber::inverse_centimetres(0.1),
            Vec::new(),
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
    fn negative_primitive_strength_is_rejected() {
        let result = DrudeLorentz::from_strength(
            1.0,
            VacuumAngularWavenumberSquared::inverse_centimetres_squared(-1.0),
            VacuumAngularWavenumber::inverse_centimetres(0.1),
            Vec::new(),
        );

        assert!(matches!(
            result,
            Err(MaterialModelError::NegativeParameter {
                name: "drude_strength",
                ..
            })
        ));
    }

    #[test]
    fn negative_damping_is_rejected_after_canonicalisation() {
        let result = DrudeLorentz::from_plasma_wavenumber(
            1.0,
            VacuumAngularWavenumber::inverse_centimetres(2.0),
            VacuumAngularWavenumber::inverse_centimetres(-0.1),
            Vec::new(),
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
    fn non_finite_epsilon_infinity_is_rejected() {
        let result = DrudeLorentz::from_plasma_wavenumber(
            f64::NAN,
            VacuumAngularWavenumber::inverse_centimetres(2.0),
            VacuumAngularWavenumber::inverse_centimetres(0.1),
            Vec::new(),
        );

        assert!(matches!(
            result,
            Err(MaterialModelError::NonFiniteParameter {
                name: "epsilon_infinity",
                ..
            })
        ));
    }

    #[test]
    fn lorentz_only_has_zero_drude_contribution() {
        let oscillator = LorentzOscillator::new(0.8, 2.7, 0.25).unwrap();

        let model = DrudeLorentz::lorentz_only(1.4, vec![oscillator]).unwrap();

        assert_eq!(model.drude_strength(), 0.0);
        assert_eq!(model.drude_damping(), 0.0);
        assert_eq!(model.plasma_wavenumber(), 0.0);
        assert_eq!(model.oscillators(), &[oscillator]);
    }

    #[test]
    fn drude_value_matches_closed_form() {
        let strength = 2.1;
        let damping = 0.4;

        let model = DrudeLorentz::new(1.0, strength, damping, Vec::new()).unwrap();

        let k0 = Complex64::new(1.3, 0.1);

        let d = k0 * k0 + Complex64::i() * damping * k0;

        let expected = -strength / d;

        assert_complex_close(model.drude_value(k0), expected);
    }

    #[test]
    fn drude_first_derivative_matches_closed_form() {
        let strength = 2.1;
        let damping = 0.4;

        let model = DrudeLorentz::new(1.0, strength, damping, Vec::new()).unwrap();

        let k0 = Complex64::new(1.3, 0.1);

        let igamma = Complex64::i() * damping;
        let d = k0 * k0 + igamma * k0;
        let d1 = 2.0 * k0 + igamma;

        let expected = strength * d1 / (d * d);

        assert_complex_close(
            model.drude_derivative_at(k0, DerivativeOrder::First),
            expected,
        );
    }

    #[test]
    fn drude_second_derivative_matches_closed_form() {
        let strength = 2.1;
        let damping = 0.4;

        let model = DrudeLorentz::new(1.0, strength, damping, Vec::new()).unwrap();

        let k0 = Complex64::new(1.3, 0.1);

        let igamma = Complex64::i() * damping;
        let d = k0 * k0 + igamma * k0;
        let d1 = 2.0 * k0 + igamma;
        let d2 = Complex64::new(2.0, 0.0);

        let expected = strength * (d2 * d - 2.0 * d1 * d1) / (d * d * d);

        assert_complex_close(
            model.drude_derivative_at(k0, DerivativeOrder::Second),
            expected,
        );
    }

    #[test]
    fn drude_third_derivative_matches_closed_form() {
        let strength = 2.1;
        let damping = 0.4;

        let model = DrudeLorentz::new(1.0, strength, damping, Vec::new()).unwrap();

        let k0 = Complex64::new(1.3, 0.1);

        let igamma = Complex64::i() * damping;
        let d = k0 * k0 + igamma * k0;
        let d1 = 2.0 * k0 + igamma;
        let d2 = Complex64::new(2.0, 0.0);

        let expected = 6.0 * strength * d1 * (d1 * d1 - d2 * d) / (d * d * d * d);

        assert_complex_close(
            model.drude_derivative_at(k0, DerivativeOrder::Third),
            expected,
        );
    }

    #[test]
    fn value_is_sum_of_constant_drude_and_lorentz_contributions() {
        let oscillator = LorentzOscillator::new(0.8, 2.7, 0.25).unwrap();

        let model = DrudeLorentz::new(1.4, 1.9, 0.3, vec![oscillator]).unwrap();

        let k0 = Complex64::new(1.1, 0.2);

        let expected = Complex64::new(1.4, 0.0) + model.drude_value(k0) + oscillator.value_at(k0);

        assert_complex_close(model.relative_permittivity_at(k0), expected);
    }

    #[test]
    fn derivative_is_sum_of_requested_contributions() {
        let oscillator = LorentzOscillator::new(0.8, 2.7, 0.25).unwrap();

        let model = DrudeLorentz::new(1.4, 1.9, 0.3, vec![oscillator]).unwrap();

        let k0 = Complex64::new(1.1, 0.2);

        for order in [
            DerivativeOrder::First,
            DerivativeOrder::Second,
            DerivativeOrder::Third,
        ] {
            let expected =
                model.drude_derivative_at(k0, order) + oscillator.derivative_at(k0, order);

            assert_complex_close(model.derivative_at(k0, order), expected);
        }
    }
}
