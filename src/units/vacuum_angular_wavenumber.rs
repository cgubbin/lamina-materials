use crate::{AngularFrequency, Frequency, ToCanonicalSpectral, VacuumWavenumber};
use lamina_units::{FrequencyUnit, InverseLengthUnit, UnitLabel};
use num_traits::{Float, FloatConst, FromPrimitive};
use std::fmt;

use super::SPEED_OF_LIGHT_CM_PER_SECOND;

/// Vacuum angular wavenumber `k₀ = 2π / λ₀ = ω / c`.
///
/// This is the canonical spectral coordinate used internally by `lamina`.
/// Its canonical unit is inverse centimetres (`cm⁻¹`).
///
/// Despite sharing the same physical dimension as [`VacuumWavenumber`], the
/// two quantities differ by a factor of `2π`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VacuumAngularWavenumber<F> {
    value: F,
    unit: InverseLengthUnit,
}

impl<F> VacuumAngularWavenumber<F> {
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

    /// Return the canonical numerical value in inverse centimetres.
    pub fn canonical(&self) -> F
    where
        F: Float + FromPrimitive,
    {
        self.as_inverse_centimetres()
    }

    pub fn vacuum_wavenumber(&self) -> VacuumWavenumber<F>
    where
        F: Float + FloatConst + FromPrimitive,
    {
        let two_pi = F::PI() + F::PI();

        VacuumWavenumber::inverse_centimetres(self.as_inverse_centimetres() / two_pi)
    }

    pub fn angular_frequency(&self) -> AngularFrequency<F>
    where
        F: Float + FromPrimitive,
    {
        let k0 = self.as_inverse_centimetres();
        let c = F::from_f64(SPEED_OF_LIGHT_CM_PER_SECOND).unwrap();

        AngularFrequency::radians_per_second(c * k0)
    }

    pub fn frequency(&self) -> Frequency<F>
    where
        F: Float + FloatConst + FromPrimitive,
    {
        let k0 = self.as_inverse_centimetres();
        let c = F::from_f64(SPEED_OF_LIGHT_CM_PER_SECOND).unwrap();
        let two_pi = F::PI() + F::PI();

        Frequency::new(c * k0 / two_pi, FrequencyUnit::Hertz)
    }

    pub(crate) fn as_inverse_centimetres(&self) -> F
    where
        F: Float + FromPrimitive,
    {
        self.unit.scale_to_inverse_centimetres::<F>() * self.value
    }
}

impl<F> ToCanonicalSpectral<F> for VacuumAngularWavenumber<F>
where
    F: Float + FromPrimitive + FloatConst,
{
    fn to_canonical_spectral(self) -> F {
        self.canonical()
    }
}

impl<F> fmt::Display for VacuumAngularWavenumber<F>
where
    F: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} {}", self.value, self.unit.symbol())
    }
}
