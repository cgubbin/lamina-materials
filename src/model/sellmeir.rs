//! Sellmeier model in canonical wavenumber form.

use crate::MaterialModelError;
use lamina_core::{
    ComplexScalar, DerivativeOrder, DifferentiableMaterial, DifferentiableMeromorphicMaterial,
    Material, MeromorphicMaterial, Sampled,
};
use lamina_units::{SpectralQuantity, Wavelength};
use num_traits::{Float, FloatConst, FromPrimitive};
use std::fmt::Debug;

/// One Sellmeier contribution in canonical angular-wavenumber form.
///
/// The contribution is
///
/// `B q² / (q² - k₀²)`,
///
/// where `B` is dimensionless and the pole angular wavenumber `q` is stored
/// canonically in inverse centimetres (`cm⁻¹`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SellmeierTerm<R> {
    coefficient: R,
    pole_wavenumber: R,
}

impl<R> SellmeierTerm<R>
where
    R: Float + Debug,
{
    /// Construct from dimensionless coefficient `B` and a spectral pole.
    ///
    /// `pole` may be supplied in any representation implementing
    /// [`SpectralQuantity`]. It is converted immediately to the canonical
    /// vacuum angular wavenumber in inverse centimetres (`cm⁻¹`).
    pub fn new<P>(coefficient: R, pole: P) -> Result<Self, MaterialModelError<R>>
    where
        R: FromPrimitive,
        P: SpectralQuantity<R>,
    {
        validate_finite("coefficient", coefficient)?;

        let pole_wavenumber = pole.vacuum_angular_wavenumber().as_inverse_centimetres();
        validate_positive("pole_wavenumber", pole_wavenumber)?;

        Ok(Self {
            coefficient,
            pole_wavenumber,
        })
    }

    /// Construct from a wavelength-domain Sellmeier pole.
    ///
    /// Literature Sellmeier relations are commonly written as
    ///
    /// ```text
    /// B λ² / (λ² - C).
    /// ```
    ///
    /// Pass `pole_wavelength = sqrt(C)`. The wavelength is converted to the
    /// canonical vacuum angular wavenumber `q = 2π / pole_wavelength` in
    /// inverse centimetres (`cm⁻¹`).
    pub fn from_wavelength_pole(
        coefficient: R,
        pole_wavelength: Wavelength<R>,
    ) -> Result<Self, MaterialModelError<R>>
    where
        R: FloatConst + FromPrimitive,
    {
        Self::new(coefficient, pole_wavelength.vacuum_angular_wavenumber())
    }

    /// Return the dimensionless coefficient `B`.
    pub fn coefficient(&self) -> R {
        self.coefficient
    }

    /// Return the canonical pole angular wavenumber in inverse centimetres.
    pub fn pole_wavenumber(&self) -> R {
        self.pole_wavenumber
    }

    fn value_at<C>(&self, k0: C) -> C
    where
        C: ComplexScalar<RealField = R> + Copy,
    {
        let q2 = self.pole_wavenumber * self.pole_wavenumber;
        C::from_real(self.coefficient * q2) / (C::from_real(q2) - k0 * k0)
    }

    fn derivative_at<C>(&self, k0: C, order: DerivativeOrder) -> C
    where
        C: ComplexScalar<RealField = R> + Copy,
    {
        let two = C::from_real(R::one() + R::one());

        let q2 = self.pole_wavenumber * self.pole_wavenumber;
        let numerator = C::from_real(self.coefficient * q2);
        let denominator = C::from_real(q2) - k0 * k0;
        let denominator2 = denominator * denominator;

        match order {
            DerivativeOrder::First => numerator * two * k0 / denominator2,
            DerivativeOrder::Second => {
                let four = two + two;
                let eight = four + four;
                let k02 = k0 * k0;
                let denominator3 = denominator2 * denominator;
                numerator * (two / denominator2 + eight * k02 / denominator3)
            }
            DerivativeOrder::Third => {
                let four = two + two;
                let eight = four + four;
                let twenty_four = eight + eight + eight;
                let forty_eight = twenty_four + twenty_four;
                let k02 = k0 * k0;
                let denominator3 = denominator2 * denominator;
                let denominator4 = denominator3 * denominator;
                numerator
                    * (twenty_four * k0 / denominator3 + forty_eight * k0 * k02 / denominator4)
            }
        }
    }
}

/// Sellmeier relative-permittivity model.
///
/// The canonical form is
///
/// `ε(k₀) = A + Σ Bⱼ qⱼ² / (qⱼ² - k₀²)`,
///
/// where `A` and `Bⱼ` are dimensionless and every pole angular wavenumber
/// `qⱼ` is expressed in inverse centimetres (`cm⁻¹`).
#[derive(Clone, Debug, PartialEq)]
pub struct Sellmeier<R> {
    constant: R,
    terms: Vec<SellmeierTerm<R>>,
}

impl<R> Sellmeier<R>
where
    R: Float + Debug,
{
    /// Construct with an explicit constant term.
    pub fn new(constant: R, terms: Vec<SellmeierTerm<R>>) -> Result<Self, MaterialModelError<R>> {
        validate_finite("constant", constant)?;
        Ok(Self { constant, terms })
    }

    /// Construct the common unit-constant form.
    pub fn with_unit_constant(terms: Vec<SellmeierTerm<R>>) -> Result<Self, MaterialModelError<R>> {
        Self::new(R::one(), terms)
    }

    /// Return the constant term.
    pub fn constant(&self) -> R {
        self.constant
    }

    /// Return the Sellmeier terms.
    pub fn terms(&self) -> &[SellmeierTerm<R>] {
        &self.terms
    }

    fn relative_permittivity_at<C>(&self, k0: C) -> C
    where
        C: ComplexScalar<RealField = R> + Copy,
    {
        self.terms
            .iter()
            .fold(C::from_real(self.constant), |epsilon, term| {
                epsilon + term.value_at(k0)
            })
    }

    fn derivative_at<C>(&self, k0: C, order: DerivativeOrder) -> C
    where
        C: ComplexScalar<RealField = R> + Copy,
    {
        self.terms.iter().fold(C::zero(), |derivative, term| {
            derivative + term.derivative_at(k0, order)
        })
    }
}

impl<R> Material for Sellmeier<R>
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

impl<R> DifferentiableMaterial for Sellmeier<R>
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

impl<R> MeromorphicMaterial for Sellmeier<R>
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

impl<R> DifferentiableMeromorphicMaterial for Sellmeier<R>
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
    use lamina_units::{LengthUnit, VacuumAngularWavenumber, VacuumWavenumber};
    use ndarray::{arr0, arr1};
    use num_complex::Complex64;

    const TOLERANCE: f64 = 1e-12;

    fn canonical_pole(value: f64) -> VacuumAngularWavenumber<f64> {
        VacuumAngularWavenumber::inverse_centimetres(value)
    }

    fn assert_complex_close(actual: Complex64, expected: Complex64) {
        let error = (actual - expected).norm();
        let scale = expected.norm().max(1.0);

        assert!(
            error <= TOLERANCE * scale,
            "actual = {actual:?}, expected = {expected:?}, relative error = {}",
            error / scale
        );
    }

    fn expected_term_value(coefficient: f64, pole: f64, k0: Complex64) -> Complex64 {
        let q2 = pole * pole;
        let numerator = coefficient * q2;
        let denominator = Complex64::new(q2, 0.0) - k0 * k0;

        numerator / denominator
    }

    fn expected_term_first(coefficient: f64, pole: f64, k0: Complex64) -> Complex64 {
        let q2 = pole * pole;
        let numerator = coefficient * q2;
        let denominator = Complex64::new(q2, 0.0) - k0 * k0;

        2.0 * numerator * k0 / denominator.powu(2)
    }

    fn expected_term_second(coefficient: f64, pole: f64, k0: Complex64) -> Complex64 {
        let q2 = pole * pole;
        let numerator = coefficient * q2;
        let denominator = Complex64::new(q2, 0.0) - k0 * k0;

        numerator * (2.0 / denominator.powu(2) + 8.0 * k0 * k0 / denominator.powu(3))
    }

    fn expected_term_third(coefficient: f64, pole: f64, k0: Complex64) -> Complex64 {
        let q2 = pole * pole;
        let numerator = coefficient * q2;
        let denominator = Complex64::new(q2, 0.0) - k0 * k0;
        let k02 = k0 * k0;

        numerator * (24.0 * k0 / denominator.powu(3) + 48.0 * k0 * k02 / denominator.powu(4))
    }

    #[test]
    fn term_constructor_canonicalises_pole() {
        let term =
            SellmeierTerm::new(0.7, VacuumAngularWavenumber::inverse_centimetres(2.8)).unwrap();

        assert_eq!(term.coefficient(), 0.7);
        assert_eq!(term.pole_wavenumber(), 2.8);
    }

    #[test]
    fn equivalent_spectral_poles_construct_same_term() {
        let angular =
            SellmeierTerm::new(0.7, VacuumAngularWavenumber::inverse_centimetres(2.8)).unwrap();

        let ordinary = SellmeierTerm::new(
            0.7,
            VacuumWavenumber::inverse_centimetres(2.8 / std::f64::consts::TAU),
        )
        .unwrap();

        assert_relative_eq!(
            angular.pole_wavenumber(),
            ordinary.pole_wavenumber(),
            epsilon = TOLERANCE,
        );

        assert_eq!(angular.coefficient(), ordinary.coefficient());
    }

    #[test]
    fn wavelength_pole_is_converted_to_angular_wavenumber() {
        let wavelength = Wavelength::new(2.0, LengthUnit::Micrometre);

        let term = SellmeierTerm::from_wavelength_pole(0.7, wavelength).unwrap();

        let wavelength_cm = 2.0e-4;
        let expected = std::f64::consts::TAU / wavelength_cm;

        assert_relative_eq!(term.pole_wavenumber(), expected, epsilon = 1e-10,);
    }

    #[test]
    fn wavelength_and_angular_wavenumber_poles_are_equivalent() {
        let wavelength = Wavelength::new(2.0, LengthUnit::Micrometre);

        let from_wavelength = SellmeierTerm::from_wavelength_pole(0.7, wavelength).unwrap();

        let from_wavenumber = SellmeierTerm::new(
            0.7,
            VacuumAngularWavenumber::inverse_centimetres(std::f64::consts::TAU / 2.0e-4),
        )
        .unwrap();

        assert_relative_eq!(
            from_wavelength.pole_wavenumber(),
            from_wavenumber.pole_wavenumber(),
            epsilon = 1e-10,
        );
    }

    #[test]
    fn zero_pole_is_rejected_after_canonicalisation() {
        let result = SellmeierTerm::new(0.7, VacuumAngularWavenumber::inverse_centimetres(0.0));

        assert!(matches!(
            result,
            Err(MaterialModelError::NonPositiveParameter {
                name: "pole_wavenumber",
                ..
            })
        ));
    }

    #[test]
    fn negative_pole_is_rejected_after_canonicalisation() {
        let result = SellmeierTerm::new(0.7, VacuumAngularWavenumber::inverse_centimetres(-1.0));

        assert!(matches!(
            result,
            Err(MaterialModelError::NonPositiveParameter {
                name: "pole_wavenumber",
                ..
            })
        ));
    }

    #[test]
    fn non_finite_coefficient_is_rejected() {
        let result = SellmeierTerm::new(f64::NAN, canonical_pole(2.8));

        assert!(matches!(
            result,
            Err(MaterialModelError::NonFiniteParameter {
                name: "coefficient",
                ..
            })
        ));
    }

    #[test]
    fn non_finite_constant_is_rejected() {
        let result = Sellmeier::new(f64::NAN, Vec::new());

        assert!(matches!(
            result,
            Err(MaterialModelError::NonFiniteParameter {
                name: "constant",
                ..
            })
        ));
    }

    #[test]
    fn with_unit_constant_sets_constant_to_one() {
        let material = Sellmeier::with_unit_constant(vec![
            SellmeierTerm::new(0.7, canonical_pole(2.8)).unwrap(),
        ])
        .unwrap();

        assert_eq!(material.constant(), 1.0);
        assert_eq!(material.terms().len(), 1);
    }

    #[test]
    fn term_value_matches_closed_form() {
        let coefficient = 0.7;
        let pole = 2.8;
        let k0 = Complex64::new(0.9, 0.15);

        let term = SellmeierTerm::new(coefficient, canonical_pole(pole)).unwrap();

        assert_complex_close(
            term.value_at(k0),
            expected_term_value(coefficient, pole, k0),
        );
    }

    #[test]
    fn term_first_derivative_matches_closed_form() {
        let coefficient = 0.7;
        let pole = 2.8;
        let k0 = Complex64::new(0.9, 0.15);

        let term = SellmeierTerm::new(coefficient, canonical_pole(pole)).unwrap();

        assert_complex_close(
            term.derivative_at(k0, DerivativeOrder::First),
            expected_term_first(coefficient, pole, k0),
        );
    }

    #[test]
    fn term_second_derivative_matches_closed_form() {
        let coefficient = 0.7;
        let pole = 2.8;
        let k0 = Complex64::new(0.9, 0.15);

        let term = SellmeierTerm::new(coefficient, canonical_pole(pole)).unwrap();

        assert_complex_close(
            term.derivative_at(k0, DerivativeOrder::Second),
            expected_term_second(coefficient, pole, k0),
        );
    }

    #[test]
    fn term_third_derivative_matches_closed_form() {
        let coefficient = 0.7;
        let pole = 2.8;
        let k0 = Complex64::new(0.9, 0.15);

        let term = SellmeierTerm::new(coefficient, canonical_pole(pole)).unwrap();

        assert_complex_close(
            term.derivative_at(k0, DerivativeOrder::Third),
            expected_term_third(coefficient, pole, k0),
        );
    }

    #[test]
    fn material_value_is_constant_plus_sum_of_terms() {
        let term_a = SellmeierTerm::new(0.7, canonical_pole(2.8)).unwrap();

        let term_b = SellmeierTerm::new(0.3, canonical_pole(4.1)).unwrap();

        let material = Sellmeier::new(1.2, vec![term_a, term_b]).unwrap();

        let k0 = Complex64::new(0.9, 0.12);

        let expected = Complex64::new(1.2, 0.0) + term_a.value_at(k0) + term_b.value_at(k0);

        assert_complex_close(material.relative_permittivity_at(k0), expected);
    }

    #[test]
    fn material_derivatives_are_sum_of_term_derivatives() {
        let term_a = SellmeierTerm::new(0.7, canonical_pole(2.8)).unwrap();

        let term_b = SellmeierTerm::new(0.3, canonical_pole(4.1)).unwrap();

        let material = Sellmeier::new(1.2, vec![term_a, term_b]).unwrap();

        let k0 = Complex64::new(0.9, 0.12);

        for order in [
            DerivativeOrder::First,
            DerivativeOrder::Second,
            DerivativeOrder::Third,
        ] {
            let expected = term_a.derivative_at(k0, order) + term_b.derivative_at(k0, order);

            assert_complex_close(material.derivative_at(k0, order), expected);
        }
    }

    #[test]
    fn real_axis_and_complex_axis_derivatives_agree() {
        let material = Sellmeier::new(
            1.0,
            vec![
                SellmeierTerm::new(0.7, canonical_pole(2.8)).unwrap(),
                SellmeierTerm::new(0.3, canonical_pole(4.1)).unwrap(),
            ],
        )
        .unwrap();

        let k0 = 0.9;

        for order in [
            DerivativeOrder::First,
            DerivativeOrder::Second,
            DerivativeOrder::Third,
        ] {
            let real = material.relative_permittivity_derivative::<_, Complex64>(arr0(k0), order);

            let complex = material.relative_permittivity_complex_derivative::<_, Complex64>(
                arr0(Complex64::new(k0, 0.0)),
                order,
            );

            assert_complex_close(real[()], complex[()]);
        }
    }

    #[test]
    fn scalar_and_sampled_derivatives_agree() {
        let material = Sellmeier::new(
            1.0,
            vec![
                SellmeierTerm::new(0.7, canonical_pole(2.8)).unwrap(),
                SellmeierTerm::new(0.3, canonical_pole(4.1)).unwrap(),
            ],
        )
        .unwrap();

        let k0 = 0.9;

        for order in [
            DerivativeOrder::First,
            DerivativeOrder::Second,
            DerivativeOrder::Third,
        ] {
            let scalar = material.relative_permittivity_derivative::<_, Complex64>(arr0(k0), order);

            let sampled =
                material.relative_permittivity_derivative::<_, Complex64>(arr1(&[k0]), order);

            assert_complex_close(sampled[0], scalar[()]);
        }
    }
}
