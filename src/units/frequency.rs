use crate::{AngularFrequency, ToCanonicalSpectral, VacuumAngularWavenumber, VacuumWavenumber};
use lamina_units::{FrequencyUnit, UnitLabel};
use num_traits::{Float, FloatConst, FromPrimitive};
use std::fmt;

use super::SPEED_OF_LIGHT_CM_PER_SECOND;

/// Ordinary frequency `f`.
///
/// A `Frequency` carries both a numerical value and its [`FrequencyUnit`].
/// It can be converted to the other vacuum spectral coordinates without
/// requiring callers to perform factors of `2π` or `c` manually.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Frequency<F> {
    value: F,
    unit: FrequencyUnit,
}

impl<F> Frequency<F> {
    pub const fn new(value: F, unit: FrequencyUnit) -> Self {
        Self { value, unit }
    }

    pub fn value(&self) -> F
    where
        F: Copy,
    {
        self.value
    }

    pub const fn unit(&self) -> FrequencyUnit {
        self.unit
    }

    pub fn angular_frequency(&self) -> AngularFrequency<F>
    where
        F: Float + FloatConst + FromPrimitive,
    {
        let frequency_hz = self.as_hertz();
        let two_pi = F::PI() + F::PI();

        AngularFrequency::radians_per_second(two_pi * frequency_hz)
    }

    pub fn vacuum_wavenumber(&self) -> VacuumWavenumber<F>
    where
        F: Float + FromPrimitive,
    {
        let frequency_hz = self.as_hertz();
        let c = F::from_f64(SPEED_OF_LIGHT_CM_PER_SECOND).unwrap();

        VacuumWavenumber::inverse_centimetres(frequency_hz / c)
    }

    pub fn vacuum_angular_wavenumber(&self) -> VacuumAngularWavenumber<F>
    where
        F: Float + FloatConst + FromPrimitive,
    {
        let frequency_hz = self.as_hertz();
        let two_pi = F::PI() + F::PI();
        let c = F::from_f64(SPEED_OF_LIGHT_CM_PER_SECOND).unwrap();

        VacuumAngularWavenumber::inverse_centimetres(two_pi * frequency_hz / c)
    }

    fn as_hertz(&self) -> F
    where
        F: Float + FromPrimitive,
    {
        self.unit.scale_to_hertz::<F>() * self.value
    }
}

impl<F> ToCanonicalSpectral<F> for Frequency<F>
where
    F: Float + FromPrimitive + FloatConst,
{
    fn to_canonical_spectral(self) -> F {
        self.vacuum_angular_wavenumber().canonical()
    }
}

impl<F> fmt::Display for Frequency<F>
where
    F: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.value, self.unit.symbol())
    }
}
