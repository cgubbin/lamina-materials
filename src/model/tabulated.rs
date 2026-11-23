//! Tabulated real-axis complex-permittivity model.

use crate::TabulatedMaterialError;
use lamina_core::{ComplexScalar, Material, Sampled};
use ndarray::Array1;
use num_complex::Complex;
use num_traits::Float;
use std::fmt::Debug;

/// Extrapolation policy for tabulated constitutive data.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Extrapolation {
    /// Clamp to the nearest endpoint value.
    #[default]
    Clamp,

    /// Continue the nearest linear segment.
    Linear,
}

/// Piecewise-linear tabulated relative-permittivity model.
///
/// Tabulated coordinates are vacuum angular wavenumbers in inverse centimetres
/// (`cm⁻¹`), and values are dimensionless complex relative permittivities.
///
/// This type intentionally implements only [`Material`]. Real-axis
/// interpolation does not define a physically meaningful meromorphic
/// continuation or intrinsic dispersive derivative contract.
#[derive(Clone, Debug, PartialEq)]
pub struct Tabulated<R> {
    vacuum_angular_wavenumber: Array1<R>,
    relative_permittivity: Array1<Complex<R>>,
    extrapolation: Extrapolation,
}

impl<R> Tabulated<R>
where
    R: Float + Debug,
{
    /// Construct validated tabulated constitutive data.
    ///
    /// `vacuum_angular_wavenumber` contains strictly increasing canonical vacuum angular
    /// wavenumbers in inverse centimetres (`cm⁻¹`). `relative_permittivity` contains
    /// the corresponding dimensionless complex relative permittivities.
    pub fn new(
        vacuum_angular_wavenumber: Array1<R>,
        relative_permittivity: Array1<Complex<R>>,
        extrapolation: Extrapolation,
    ) -> Result<Self, TabulatedMaterialError<R>> {
        if vacuum_angular_wavenumber.len() != relative_permittivity.len() {
            return Err(TabulatedMaterialError::LengthMismatch {
                coordinate_len: vacuum_angular_wavenumber.len(),
                value_len: relative_permittivity.len(),
            });
        }

        if vacuum_angular_wavenumber.len() < 2 {
            return Err(TabulatedMaterialError::TooFewSamples {
                len: vacuum_angular_wavenumber.len(),
            });
        }

        for (index, &coordinate) in vacuum_angular_wavenumber.iter().enumerate() {
            if !coordinate.is_finite() {
                return Err(TabulatedMaterialError::NonFiniteCoordinate {
                    index,
                    value: coordinate,
                });
            }
        }

        for (index, value) in relative_permittivity.iter().enumerate() {
            if !value.re.is_finite() || !value.im.is_finite() {
                return Err(TabulatedMaterialError::NonFiniteValue { index });
            }
        }

        for index in 1..vacuum_angular_wavenumber.len() {
            let left = vacuum_angular_wavenumber[index - 1];
            let right = vacuum_angular_wavenumber[index];

            if right <= left {
                return Err(TabulatedMaterialError::NonIncreasingCoordinates {
                    left,
                    right,
                    right_index: index,
                });
            }
        }

        Ok(Self {
            vacuum_angular_wavenumber,
            relative_permittivity,
            extrapolation,
        })
    }

    /// Return tabulated coordinates.
    pub fn vacuum_angular_wavenumber(&self) -> &Array1<R> {
        &self.vacuum_angular_wavenumber
    }

    /// Return tabulated complex relative permittivity.
    pub fn values(&self) -> &Array1<Complex<R>> {
        &self.relative_permittivity
    }

    /// Return the extrapolation policy.
    pub fn extrapolation(&self) -> Extrapolation {
        self.extrapolation
    }

    fn interpolate(&self, x: R) -> Complex<R> {
        let last = self.vacuum_angular_wavenumber.len() - 1;

        if x <= self.vacuum_angular_wavenumber[0] {
            return match self.extrapolation {
                Extrapolation::Clamp => self.relative_permittivity[0],
                Extrapolation::Linear => self.interpolate_segment(0, x),
            };
        }

        if x >= self.vacuum_angular_wavenumber[last] {
            return match self.extrapolation {
                Extrapolation::Clamp => self.relative_permittivity[last],
                Extrapolation::Linear => self.interpolate_segment(last - 1, x),
            };
        }

        let coordinates = self
            .vacuum_angular_wavenumber
            .as_slice()
            .expect("owned Array1 is contiguous");

        let right = coordinates.partition_point(|value| *value < x);
        self.interpolate_segment(right - 1, x)
    }

    fn interpolate_segment(&self, left: usize, x: R) -> Complex<R> {
        let right = left + 1;
        let x0 = self.vacuum_angular_wavenumber[left];
        let x1 = self.vacuum_angular_wavenumber[right];
        let fraction = (x - x0) / (x1 - x0);
        let y0 = self.relative_permittivity[left];
        let y1 = self.relative_permittivity[right];

        y0 + (y1 - y0) * fraction
    }
}

impl<R> Material for Tabulated<R>
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
        vacuum_angular_wavenumber.map(|k0| {
            let epsilon = self.interpolate(k0);
            C::from_parts(epsilon.re, epsilon.im)
        })
    }
}
