use crate::{AngularFrequency, Frequency, ToCanonicalSpectral, VacuumAngularWavenumber};
use lamina_units::{FrequencyUnit, InverseLengthUnit, UnitLabel};
use num_traits::{Float, FloatConst, FromPrimitive};
use std::fmt;

use super::SPEED_OF_LIGHT_CM_PER_SECOND;

/// Vacuum wavenumber `ν̃ = 1 / λ₀`.
///
/// This is the ordinary inverse-wavelength coordinate commonly quoted in
/// spectroscopy, for example in `cm⁻¹`.
///
/// It is distinct from [`VacuumAngularWavenumber`], which is larger by a
/// factor of `2π`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VacuumWavenumber<F> {
    value: F,
    unit: InverseLengthUnit,
}

impl<F> VacuumWavenumber<F> {
    pub const fn new(value: F, unit: InverseLengthUnit) -> Self {
        Self { value, unit }
    }

    pub fn inverse_centimetres(value: F) -> Self {
        Self::new(value, InverseLengthUnit::PerCentimetre)
    }

    pub fn value(&self) -> F
    where
        F: Copy,
    {
        self.value
    }

    pub const fn unit(&self) -> InverseLengthUnit {
        self.unit
    }

    pub fn frequency(&self) -> Frequency<F>
    where
        F: Float + FromPrimitive,
    {
        let wavenumber = self.as_inverse_centimetres();
        let c = F::from_f64(SPEED_OF_LIGHT_CM_PER_SECOND).unwrap();

        Frequency::new(c * wavenumber, FrequencyUnit::Hertz)
    }

    pub fn angular_frequency(&self) -> AngularFrequency<F>
    where
        F: Float + FloatConst + FromPrimitive,
    {
        let wavenumber = self.as_inverse_centimetres();
        let two_pi = F::PI() + F::PI();
        let c = F::from_f64(SPEED_OF_LIGHT_CM_PER_SECOND).unwrap();

        AngularFrequency::radians_per_second(two_pi * c * wavenumber)
    }

    pub fn vacuum_angular_wavenumber(&self) -> VacuumAngularWavenumber<F>
    where
        F: Float + FloatConst + FromPrimitive,
    {
        let two_pi = F::PI() + F::PI();

        VacuumAngularWavenumber::inverse_centimetres(two_pi * self.as_inverse_centimetres())
    }

    pub(crate) fn as_inverse_centimetres(&self) -> F
    where
        F: Float + FromPrimitive,
    {
        self.unit.scale_to_inverse_centimetres::<F>() * self.value
    }
}

impl<F> ToCanonicalSpectral<F> for VacuumWavenumber<F>
where
    F: Float + FromPrimitive + FloatConst,
{
    fn to_canonical_spectral(self) -> F {
        self.vacuum_angular_wavenumber().canonical()
    }
}

impl<F> fmt::Display for VacuumWavenumber<F>
where
    F: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.value, self.unit.symbol())
    }
}
