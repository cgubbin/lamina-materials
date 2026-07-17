//! Sellmeier model in canonical wavenumber form.

use crate::{MaterialModelError, Wavelength, model::drude_lorentz::convert_derivative_variable};
use num_traits::{Float, FloatConst};
use std::fmt::Debug;
use tmm_core::{
    ComplexScalar, DerivativeOrder, DifferentiableMaterial, DifferentiableMeromorphicMaterial,
    Material, MeromorphicMaterial, Sampled, SpectralVariable,
};

/// One canonical Sellmeier term.
///
/// The contribution is `B q² / (q² - k₀²)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SellmeierTerm<R> {
    coefficient: R,
    pole_wavenumber: R,
}

impl<R> SellmeierTerm<R>
where
    R: Float + Debug,
{
    /// Construct a term from dimensionless coefficient `B` and pole wavenumber `q`.
    pub fn new(coefficient: R, pole_wavenumber: R) -> Result<Self, MaterialModelError<R>> {
        validate_finite("coefficient", coefficient)?;
        validate_positive("pole_wavenumber", pole_wavenumber)?;

        Ok(Self {
            coefficient,
            pole_wavenumber,
        })
    }

    // Construct from the wavelength-domain Sellmeier pole.
    ///
    /// For the usual term
    ///
    /// ```text
    /// B λ² / (λ² - C),
    /// ```
    ///
    /// pass `pole_wavelength = sqrt(C)` in the units used by the source.
    pub fn from_wavelength_pole(
        coefficient: R,
        pole_wavelength: Wavelength<R>,
    ) -> Result<Self, MaterialModelError<R>>
    where
        R: FloatConst,
    {
        Self::new(coefficient, pole_wavelength.vacuum_angular_wavenumber())
    }

    /// Return coefficient `B`.
    pub fn coefficient(&self) -> R {
        self.coefficient
    }

    /// Return pole wavenumber `q`.
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

    fn derivatives_at<C>(&self, k0: C) -> (C, C)
    where
        C: ComplexScalar<RealField = R> + Copy,
    {
        let two = C::from_real(R::one() + R::one());
        let eight = two + two + two + two;
        let q2 = self.pole_wavenumber * self.pole_wavenumber;
        let numerator = C::from_real(self.coefficient * q2);
        let denominator = C::from_real(q2) - k0 * k0;

        let first = numerator * two * k0 / (denominator * denominator);

        let second = numerator
            * (two / (denominator * denominator)
                + eight * k0 * k0 / (denominator * denominator * denominator));

        (first, second)
    }
}

/// Sellmeier relative-permittivity model.
///
/// `ε(k₀) = A + Σ Bⱼ qⱼ²/(qⱼ²-k₀²)`.
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

    fn derivative_at<C>(&self, k0: C, order: DerivativeOrder, variable: SpectralVariable) -> C
    where
        C: ComplexScalar<RealField = R> + Copy,
    {
        let (first, second) =
            self.terms
                .iter()
                .fold((C::zero(), C::zero()), |(first, second), term| {
                    let (a, b) = term.derivatives_at(k0);
                    (first + a, second + b)
                });

        convert_derivative_variable(k0, first, second, order, variable)
    }
}

impl<R> Material for Sellmeier<R>
where
    R: Float + Debug,
{
    type Real = R;

    fn relative_permeability<I, C>(&self, vacuum_wavenumber: I) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = R> + Copy,
        I: Sampled<Elem = R>,
    {
        vacuum_wavenumber.map(|_| C::one())
    }

    fn relative_permittivity<I, C>(&self, vacuum_wavenumber: I) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = R> + Copy,
        I: Sampled<Elem = R>,
    {
        vacuum_wavenumber.map(|k0| self.relative_permittivity_at::<C>(C::from_real(k0)))
    }
}

impl<R> DifferentiableMaterial for Sellmeier<R>
where
    R: Float + Debug,
{
    fn relative_permittivity_derivative<I, C>(
        &self,
        vacuum_wavenumber: I,
        order: DerivativeOrder,
        variable: SpectralVariable,
    ) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = R> + Copy,
        I: Sampled<Elem = R>,
    {
        vacuum_wavenumber.map(|k0| self.derivative_at::<C>(C::from_real(k0), order, variable))
    }
}

impl<R> MeromorphicMaterial for Sellmeier<R>
where
    R: Float + Debug,
{
    fn relative_permeability_complex<I, C>(&self, vacuum_wavenumber: I) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = R> + Copy,
        I: Sampled<Elem = C>,
    {
        vacuum_wavenumber.map(|_| C::one())
    }

    fn relative_permittivity_complex<I, C>(&self, vacuum_wavenumber: I) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = R> + Copy,
        I: Sampled<Elem = C>,
    {
        vacuum_wavenumber.map(|k0| self.relative_permittivity_at(k0))
    }
}

impl<R> DifferentiableMeromorphicMaterial for Sellmeier<R>
where
    R: Float + Debug,
{
    fn relative_permittivity_complex_derivative<I, C>(
        &self,
        vacuum_wavenumber: I,
        order: DerivativeOrder,
        variable: SpectralVariable,
    ) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = R> + Copy,
        I: Sampled<Elem = C>,
    {
        vacuum_wavenumber.map(|k0| self.derivative_at(k0, order, variable))
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
