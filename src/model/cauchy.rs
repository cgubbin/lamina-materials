//! Cauchy refractive-index model in canonical wavenumber form.

use crate::MaterialModelError;
use lamina_core::{
    ComplexScalar, DerivativeOrder, DifferentiableMaterial, DifferentiableMeromorphicMaterial,
    Material, MeromorphicMaterial, Sampled,
};
use num_traits::Float;
use std::fmt::Debug;

/// Cauchy refractive-index model in canonical angular-wavenumber form.
///
/// Coefficients define
///
/// `n(k₀) = a₀ + a₁k₀² + a₂k₀⁴ + ...`,
///
/// with relative permittivity `ε = n²`.
///
/// The vacuum angular wavenumber `k₀` is expressed in inverse centimetres
/// (`cm⁻¹`). To keep `n` dimensionless, coefficient `aᵢ` has units
/// `cm^(2i)`. In particular, `a₀` is dimensionless.
#[derive(Clone, Debug, PartialEq)]
pub struct Cauchy<R> {
    coefficients: Vec<R>,
}

impl<R> Cauchy<R>
where
    R: Float + Debug,
{
    /// Construct from coefficients ordered by ascending powers of `k₀²`.
    ///
    /// For coefficient index `i`, `aᵢ` has units `cm^(2i)` because `k₀` is
    /// expressed in `cm⁻¹`.
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
    ///
    /// Coefficient `aᵢ` has units `cm^(2i)`.
    pub fn coefficients(&self) -> &[R] {
        &self.coefficients
    }

    fn refractive_index_at<C>(&self, k0: C) -> C
    where
        C: ComplexScalar<RealField = R> + Copy,
    {
        let x = k0 * k0;

        self.coefficients
            .iter()
            .rev()
            .fold(C::zero(), |value, &coefficient| {
                value * x + C::from_real(coefficient)
            })
    }

    fn relative_permittivity_at<C>(&self, k0: C) -> C
    where
        C: ComplexScalar<RealField = R> + Copy,
    {
        let n = self.refractive_index_at(k0);
        n * n
    }

    fn refractive_index_derivative_at<C>(&self, k0: C, order: DerivativeOrder) -> C
    where
        C: ComplexScalar<RealField = R> + Copy,
    {
        let one = R::one();
        let two_r = one + one;

        let two = C::from_real(two_r);
        let four = two + two;
        let eight = four + four;
        let twelve = eight + four;

        let x = k0 * k0;

        match order {
            DerivativeOrder::First => {
                let mut dn_dx = C::zero();

                for (index, &coefficient) in self.coefficients.iter().enumerate().skip(1) {
                    let index_r = R::from(index).expect("usize converts to float");

                    dn_dx = dn_dx + C::from_real(coefficient * index_r) * powi(x, index - 1);
                }

                two * k0 * dn_dx
            }

            DerivativeOrder::Second => {
                let mut dn_dx = C::zero();
                let mut d2n_dx2 = C::zero();

                for (index, &coefficient) in self.coefficients.iter().enumerate().skip(1) {
                    let index_r = R::from(index).expect("usize converts to float");

                    dn_dx = dn_dx + C::from_real(coefficient * index_r) * powi(x, index - 1);

                    if index >= 2 {
                        d2n_dx2 = d2n_dx2
                            + C::from_real(coefficient * index_r * (index_r - one))
                                * powi(x, index - 2);
                    }
                }

                two * dn_dx + four * x * d2n_dx2
            }

            DerivativeOrder::Third => {
                let mut d2n_dx2 = C::zero();
                let mut d3n_dx3 = C::zero();

                for (index, &coefficient) in self.coefficients.iter().enumerate().skip(2) {
                    let index_r = R::from(index).expect("usize converts to float");

                    d2n_dx2 = d2n_dx2
                        + C::from_real(coefficient * index_r * (index_r - one))
                            * powi(x, index - 2);

                    if index >= 3 {
                        d3n_dx3 = d3n_dx3
                            + C::from_real(
                                coefficient * index_r * (index_r - one) * (index_r - two_r),
                            ) * powi(x, index - 3);
                    }
                }

                twelve * k0 * d2n_dx2 + eight * k0 * x * d3n_dx3
            }
        }
    }

    fn derivative_at<C>(&self, k0: C, order: DerivativeOrder) -> C
    where
        C: ComplexScalar<RealField = R> + Copy,
    {
        let two = C::from_real(R::one() + R::one());

        match order {
            DerivativeOrder::First => {
                let n = self.refractive_index_at(k0);
                let n_first = self.refractive_index_derivative_at(k0, DerivativeOrder::First);

                two * n * n_first
            }

            DerivativeOrder::Second => {
                let n = self.refractive_index_at(k0);

                let n_first = self.refractive_index_derivative_at(k0, DerivativeOrder::First);

                let n_second = self.refractive_index_derivative_at(k0, DerivativeOrder::Second);

                two * (n_first * n_first + n * n_second)
            }

            DerivativeOrder::Third => {
                let n = self.refractive_index_at(k0);

                let n_first = self.refractive_index_derivative_at(k0, DerivativeOrder::First);

                let n_second = self.refractive_index_derivative_at(k0, DerivativeOrder::Second);

                let n_third = self.refractive_index_derivative_at(k0, DerivativeOrder::Third);

                let three = C::from_real(R::one() + R::one() + R::one());

                two * (three * n_first * n_second + n * n_third)
            }
        }
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
        vacuum_wavenumber.map(|k0| self.relative_permittivity_at::<C>(C::from_real(k0)))
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
    ) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = R> + Copy,
        I: Sampled<Elem = R>,
    {
        vacuum_wavenumber.map(|k0| self.derivative_at::<C>(C::from_real(k0), order))
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
        vacuum_wavenumber.map(|k0| self.relative_permittivity_at(k0))
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
    ) -> I::Mapped<C>
    where
        C: ComplexScalar<RealField = R> + Copy,
        I: Sampled<Elem = C>,
    {
        vacuum_wavenumber.map(|k0| self.derivative_at(k0, order))
    }
}

fn powi<C>(value: C, exponent: usize) -> C
where
    C: ComplexScalar + Copy,
{
    (0..exponent).fold(C::one(), |result, _| result * value)
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

    #[test]
    fn refractive_index_value_matches_polynomial() {
        let coefficients = [1.4, 0.3, -0.07, 0.012];

        let model = Cauchy::new(coefficients.to_vec()).unwrap();

        let k = Complex64::new(0.8, 0.15);
        let k2 = k * k;
        let k4 = k2 * k2;
        let k6 = k4 * k2;

        let expected = Complex64::new(coefficients[0], 0.0)
            + coefficients[1] * k2
            + coefficients[2] * k4
            + coefficients[3] * k6;

        assert_complex_close(model.refractive_index_at(k), expected);
    }

    #[test]
    fn refractive_index_derivatives_match_polynomial() {
        let a0 = 1.4;
        let a1 = 0.3;
        let a2 = -0.07;
        let a3 = 0.012;

        let model = Cauchy::new(vec![a0, a1, a2, a3]).unwrap();

        let k = Complex64::new(0.8, 0.15);

        let k2 = k * k;
        let k3 = k2 * k;
        let k4 = k2 * k2;
        let k5 = k4 * k;

        let expected_first = 2.0 * a1 * k + 4.0 * a2 * k3 + 6.0 * a3 * k5;

        let expected_second = Complex64::new(2.0 * a1, 0.0) + 12.0 * a2 * k2 + 30.0 * a3 * k4;

        let expected_third = 24.0 * a2 * k + 120.0 * a3 * k3;

        assert_complex_close(
            model.refractive_index_derivative_at(k, DerivativeOrder::First),
            expected_first,
        );

        assert_complex_close(
            model.refractive_index_derivative_at(k, DerivativeOrder::Second),
            expected_second,
        );

        assert_complex_close(
            model.refractive_index_derivative_at(k, DerivativeOrder::Third),
            expected_third,
        );
    }

    #[test]
    fn relative_permittivity_is_square_of_refractive_index() {
        let model = Cauchy::new(vec![1.4, 0.3, -0.07, 0.012]).unwrap();

        let k = Complex64::new(0.8, 0.15);

        let n = model.refractive_index_at(k);

        assert_complex_close(model.relative_permittivity_at(k), n * n);
    }

    #[test]
    fn permittivity_derivatives_follow_square_rule() {
        let model = Cauchy::new(vec![1.4, 0.3, -0.07, 0.012]).unwrap();

        let k = Complex64::new(0.8, 0.15);

        let n = model.refractive_index_at(k);

        let n1 = model.refractive_index_derivative_at(k, DerivativeOrder::First);

        let n2 = model.refractive_index_derivative_at(k, DerivativeOrder::Second);

        let n3 = model.refractive_index_derivative_at(k, DerivativeOrder::Third);

        assert_complex_close(model.derivative_at(k, DerivativeOrder::First), 2.0 * n * n1);

        assert_complex_close(
            model.derivative_at(k, DerivativeOrder::Second),
            2.0 * (n1 * n1 + n * n2),
        );

        assert_complex_close(
            model.derivative_at(k, DerivativeOrder::Third),
            6.0 * n1 * n2 + 2.0 * n * n3,
        );
    }
}
