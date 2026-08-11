use crate::ToCanonicalSpectralSquared;

/// A squared inverse-length quantity.
///
/// The stored value is canonical and therefore expressed in `cm⁻²`.
///
/// This type is useful for primitive material parameters such as Drude
/// strength `Ωᴅ²` and Lorentz oscillator strength `S`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InverseLengthSquared<F> {
    canonical: F,
}

impl<F> InverseLengthSquared<F> {
    /// Construct directly from a value in `cm⁻²`.
    pub const fn inverse_centimetres_squared(value: F) -> Self {
        Self { canonical: value }
    }

    /// Return the canonical value in `cm⁻²`.
    pub fn canonical(self) -> F {
        self.canonical
    }
}

impl<F> ToCanonicalSpectralSquared<F> for InverseLengthSquared<F> {
    fn to_canonical_spectral_squared(self) -> F {
        self.canonical
    }
}
