use lamina_units::{LengthUnit, UnitLabel};
use num_traits::{Float, FloatConst, FromPrimitive, Zero};
use std::fmt;

use crate::{ToCanonicalSpectral, VacuumAngularWavenumber, VacuumWavenumber};

/// A vacuum wavelength with an explicit physical length unit.
///
/// `Wavelength` is primarily used to construct material parameters from
/// wavelength-domain literature values. Conversion to the spectral coordinate
/// used by `lamina-core` is performed by [`Self::vacuum_angular_wavenumber`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wavelength<F> {
    value: F,
    unit: LengthUnit,
}

impl<F> Wavelength<F> {
    pub(crate) fn new(value: F, unit: LengthUnit) -> Self {
        Self { value, unit }
    }

    pub fn centimetres(value: F) -> Self {
        Self::new(value, LengthUnit::Centimetre)
    }

    pub fn millimetres(value: F) -> Self {
        Self::new(value, LengthUnit::Millimetre)
    }

    pub fn micrometres(value: F) -> Self {
        Self::new(value, LengthUnit::Micrometre)
    }

    pub fn nanometres(value: F) -> Self {
        Self::new(value, LengthUnit::Nanometre)
    }

    pub fn vacuum_wavenumber(&self) -> VacuumWavenumber<F>
    where
        F: Float + FromPrimitive,
    {
        VacuumWavenumber::inverse_centimetres(F::one() / self.as_cm())
    }

    pub fn vacuum_angular_wavenumber(&self) -> VacuumAngularWavenumber<F>
    where
        F: Float + FloatConst + FromPrimitive,
    {
        let two_pi = F::PI() + F::PI();

        VacuumAngularWavenumber::inverse_centimetres(two_pi / self.as_cm())
    }

    pub fn value(&self) -> F
    where
        F: Copy,
    {
        self.value
    }

    pub fn unit(&self) -> LengthUnit {
        self.unit
    }

    pub(crate) fn into_parts(self) -> (F, LengthUnit) {
        (self.value, self.unit)
    }

    pub(crate) fn as_cm(&self) -> F
    where
        F: Float + FromPrimitive,
    {
        self.unit.scale_to_centimetres::<F>() * self.value
    }
}

impl<F> ToCanonicalSpectral<F> for Wavelength<F>
where
    F: Float + FromPrimitive + FloatConst,
{
    fn to_canonical_spectral(self) -> F {
        self.vacuum_angular_wavenumber().canonical()
    }
}

impl<F> fmt::Display for Wavelength<F>
where
    F: Float + std::fmt::Debug + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.value, self.unit.symbol())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn micrometre_wavelength_converts_to_inverse_centimetres() {
        let wavelength = Wavelength::micrometres(2.0);

        assert_relative_eq!(
            wavelength.vacuum_angular_wavenumber(),
            std::f64::consts::TAU / 2.0e-4,
            epsilon = 1.0e-10,
        );
    }

    #[test]
    fn wavelength_units_give_same_canonical_wavenumber() {
        let expected = std::f64::consts::TAU / 1.5e-4;

        for actual in [
            Wavelength::centimetres(1.5e-4).vacuum_angular_wavenumber(),
            Wavelength::millimetres(1.5e-3).vacuum_angular_wavenumber(),
            Wavelength::micrometres(1.5).vacuum_angular_wavenumber(),
            Wavelength::nanometres(1500.0).vacuum_angular_wavenumber(),
        ] {
            assert_relative_eq!(actual, expected, epsilon = 1.0e-10);
        }
    }
}
