//! Metadata associated with named material records.

use lamina_units::SpectralQuantity;
use num_traits::{Float, FromPrimitive};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectralDomain<R> {
    minimum: R,
    maximum: R,
}

impl<R> SpectralDomain<R>
where
    R: Float,
{
    /// Construct directly from canonical vacuum angular wavenumbers in cm⁻¹.
    pub fn new(minimum: R, maximum: R) -> Option<Self> {
        if !minimum.is_finite() || !maximum.is_finite() || minimum > maximum {
            return None;
        }

        Some(Self { minimum, maximum })
    }

    /// Construct from endpoints expressed in any supported spectral coordinate.
    ///
    /// The endpoints are converted to canonical vacuum angular wavenumber
    /// in inverse centimetres and ordered in the canonical coordinate.
    pub fn from_spectral<S, T>(first: S, second: T) -> Option<Self>
    where
        R: FromPrimitive,
        S: SpectralQuantity<R>,
        T: SpectralQuantity<R>,
    {
        let first = first.vacuum_angular_wavenumber().as_inverse_centimetres();

        let second = second.vacuum_angular_wavenumber().as_inverse_centimetres();

        if !first.is_finite() || !second.is_finite() {
            return None;
        }

        let (minimum, maximum) = if first <= second {
            (first, second)
        } else {
            (second, first)
        };

        Some(Self { minimum, maximum })
    }

    pub fn minimum(self) -> R {
        self.minimum
    }

    pub fn maximum(self) -> R {
        self.maximum
    }

    pub fn contains(self, value: R) -> bool {
        value >= self.minimum && value <= self.maximum
    }
}

pub trait MaterialRecord {
    type Real;
    fn name(&self) -> &'static str;
    fn identifier(&self) -> &'static str;
    fn source(&self) -> &'static str;
    fn spectral_domain(&self) -> Option<SpectralDomain<Self::Real>>;
    fn notes(&self) -> &'static str {
        ""
    }
}
