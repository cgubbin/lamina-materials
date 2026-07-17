//! Lightweight wavelength units used by literature-convention constructors.

use num_traits::{Float, FloatConst};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LengthUnit {
    Metre,
    Millimetre,
    Micrometre,
    Nanometre,
    Centimetre,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Wavelength<R> {
    value: R,
    unit: LengthUnit,
}

impl<R: Float> Wavelength<R> {
    pub fn new(value: R, unit: LengthUnit) -> Option<Self> {
        if !value.is_finite() || value <= R::zero() {
            None
        } else {
            Some(Self { value, unit })
        }
    }
    pub fn value(self) -> R {
        self.value
    }
    pub fn unit(self) -> LengthUnit {
        self.unit
    }
    pub fn in_centimetres(self) -> R {
        let scale = match self.unit {
            LengthUnit::Metre => R::from(100.0).unwrap(),
            LengthUnit::Millimetre => R::from(0.1).unwrap(),
            LengthUnit::Micrometre => R::from(1.0e-4).unwrap(),
            LengthUnit::Nanometre => R::from(1.0e-7).unwrap(),
            LengthUnit::Centimetre => R::one(),
        };
        self.value * scale
    }
    pub fn vacuum_angular_wavenumber(self) -> R
    where
        R: FloatConst,
    {
        (R::one() + R::one()) * R::PI() / self.in_centimetres()
    }
}
