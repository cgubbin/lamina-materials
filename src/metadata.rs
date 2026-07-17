//! Metadata associated with named material records.

use num_traits::Float;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpectralDomain<R> {
    minimum: R,
    maximum: R,
}

impl<R: Float> SpectralDomain<R> {
    pub fn new(minimum: R, maximum: R) -> Option<Self> {
        if !minimum.is_finite() || !maximum.is_finite() || minimum > maximum {
            return None;
        }
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
