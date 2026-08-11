use num_traits::Float;

use crate::units::InverseLengthSquared;

pub trait ToCanonicalSpectral<R> {
    fn to_canonical_spectral(self) -> R;
}

pub trait ToCanonicalSpectralSquared<R> {
    fn to_canonical_spectral_squared(self) -> R;
}

/// Convenience operations for quantities convertible to canonical spectral
/// coordinates.
pub trait SpectralExt<R>: ToCanonicalSpectral<R> + Sized
where
    R: Float,
{
    /// Convert to canonical vacuum angular wavenumber and square it.
    ///
    /// The result is stored in `cm⁻²`.
    fn squared(self) -> InverseLengthSquared<R> {
        let value = self.to_canonical_spectral();

        InverseLengthSquared::inverse_centimetres_squared(value * value)
    }
}

impl<R, T> SpectralExt<R> for T
where
    R: Float,
    T: ToCanonicalSpectral<R>,
{
}
