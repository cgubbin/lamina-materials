//! Cauchy refractive-index model in canonical wavenumber form.

use crate::{MaterialModelError, model::drude_lorentz::convert_derivative_variable};
use num_traits::Float;
use std::fmt::Debug;
use tmm_core::{
    ComplexScalar, DerivativeOrder, DifferentiableMaterial, DifferentiableMeromorphicMaterial,
    Material, MeromorphicMaterial, Sampled, SpectralVariable,
};

/// Cauchy refractive-index model.
///
/// Coefficients define `n(k₀) = a₀ + a₁k₀² + a₂k₀⁴ + ...`; relative
/// permittivity is `ε = n²`.
#[derive(Clone, Debug, PartialEq)]
pub struct Cauchy<R> {
    coefficients: Vec<R>,
}

impl<R> Cauchy<R>
where
    R: Float + Debug,
{
    /// Construct from coefficients ordered by ascending powers of `k₀²`.
    pub fn new(coefficients: Vec<R>) -> Result<Self, MaterialModelError<R>> {
        if coefficients.is_empty() {
            return Err(MaterialModelError::EmptyCoefficients {
                name: "coefficients",
            });
        }

        for &coefficient in &coefficients {
            if !coefficient.is_finite() {
                return Err(MaterialModelError::NonFiniteParameter {
                    name: "coefficient",
                    value: coefficient,
                });
            }
        }

        Ok(Self { coefficients })
    }

    /// Return coefficients ordered by ascending powers of `k₀²`.
    pub fn coefficients(&self) -> &[R] {
        &self.coefficients
    }

    fn refractive_index_and_derivatives_at<C>(&self, k0: C) -> (C, C, C)
    where
        C: ComplexScalar<RealField = R> + Copy,
    {
        let two_r = R::one() + R::one();
        let x = k0 * k0;

        let mut n = C::zero();
        let mut dn_dx = C::zero();
        let mut d2n_dx2 = C::zero();

        for (index, &coefficient) in self.coefficients.iter().enumerate() {
            n = n + C::from_real(coefficient) * powi(x, index);

            if index >= 1 {
                let index_r = R::from(index).expect("usize converts to float");
                dn_dx = dn_dx + C::from_real(coefficient * index_r) * powi(x, index - 1);
            }

            if index >= 2 {
                let index_r = R::from(index).expect("usize converts to float");
                d2n_dx2 = d2n_dx2
                    + C::from_real(coefficient * index_r * (index_r - R::one()))
                        * powi(x, index - 2);
            }
        }

        let dn_dk = C::from_real(two_r) * k0 * dn_dx;
        let d2n_dk2 = C::from_real(two_r) * dn_dx + C::from_real(two_r * two_r) * k0 * k0 * d2n_dx2;

        (n, dn_dk, d2n_dk2)
    }

    fn relative_permittivity_and_derivatives_at<C>(&self, k0: C) -> (C, C, C)
    where
        C: ComplexScalar<RealField = R> + Copy,
    {
        let two = C::from_real(R::one() + R::one());
        let (n, n_first, n_second) = self.refractive_index_and_derivatives_at(k0);

        (
            n * n,
            two * n * n_first,
            two * (n_first * n_first + n * n_second),
        )
    }

    fn derivative_at<C>(&self, k0: C, order: DerivativeOrder, variable: SpectralVariable) -> C
    where
        C: ComplexScalar<RealField = R> + Copy,
    {
        let (_, first, second) = self.relative_permittivity_and_derivatives_at(k0);

        convert_derivative_variable(k0, first, second, order, variable)
    }
}

impl<R> Material for Cauchy<R>
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
        vacuum_wavenumber.map(|k0| {
            self.relative_permittivity_and_derivatives_at::<C>(C::from_real(k0))
                .0
        })
    }
}

impl<R> DifferentiableMaterial for Cauchy<R>
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

impl<R> MeromorphicMaterial for Cauchy<R>
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
        vacuum_wavenumber.map(|k0| self.relative_permittivity_and_derivatives_at(k0).0)
    }
}

impl<R> DifferentiableMeromorphicMaterial for Cauchy<R>
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

fn powi<C>(value: C, exponent: usize) -> C
where
    C: ComplexScalar + Copy,
{
    (0..exponent).fold(C::one(), |result, _| result * value)
}
