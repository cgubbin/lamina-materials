use crate::{Frequency, ToCanonicalSpectral, VacuumAngularWavenumber, VacuumWavenumber};
use lamina_units::{AngularFrequencyUnit, FrequencyUnit, UnitLabel};
use num_traits::{Float, FloatConst, FromPrimitive};
use std::fmt;

use super::SPEED_OF_LIGHT_CM_PER_SECOND;

/// Angular frequency `ω`.
///
/// Angular frequency is related to ordinary frequency by `ω = 2πf`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AngularFrequency<F> {
    value: F,
    unit: AngularFrequencyUnit,
}

impl<F> AngularFrequency<F> {
    pub const fn new(value: F, unit: AngularFrequencyUnit) -> Self {
        Self { value, unit }
    }

    pub fn radians_per_second(value: F) -> Self {
        Self::new(value, AngularFrequencyUnit::RadiansPerSecond)
    }

    pub fn value(&self) -> F
    where
        F: Copy,
    {
        self.value
    }

    pub const fn unit(&self) -> AngularFrequencyUnit {
        self.unit
    }

    pub fn frequency(&self) -> Frequency<F>
    where
        F: Float + FloatConst + FromPrimitive,
    {
        let omega = self.as_radians_per_second();
        let two_pi = F::PI() + F::PI();

        Frequency::new(omega / two_pi, FrequencyUnit::Hertz)
    }

    pub fn vacuum_wavenumber(&self) -> VacuumWavenumber<F>
    where
        F: Float + FloatConst + FromPrimitive,
    {
        let omega = self.as_radians_per_second();
        let two_pi = F::PI() + F::PI();
        let c = F::from_f64(SPEED_OF_LIGHT_CM_PER_SECOND).unwrap();

        VacuumWavenumber::inverse_centimetres(omega / (two_pi * c))
    }

    pub fn vacuum_angular_wavenumber(&self) -> VacuumAngularWavenumber<F>
    where
        F: Float + FromPrimitive,
    {
        let omega = self.as_radians_per_second();
        let c = F::from_f64(SPEED_OF_LIGHT_CM_PER_SECOND).unwrap();

        VacuumAngularWavenumber::inverse_centimetres(omega / c)
    }

    fn as_radians_per_second(&self) -> F
    where
        F: Float + FromPrimitive,
    {
        self.unit.scale_to_radians_per_second::<F>() * self.value
    }
}

impl<F> ToCanonicalSpectral<F> for AngularFrequency<F>
where
    F: Float + FromPrimitive,
{
    fn to_canonical_spectral(self) -> F {
        self.vacuum_angular_wavenumber().canonical()
    }
}

impl<F> fmt::Display for AngularFrequency<F>
where
    F: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.value, self.unit.symbol())
    }
}
