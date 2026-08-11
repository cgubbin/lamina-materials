mod angular_frequency;
mod frequency;
mod inverse_length_squared;
mod spectral;
mod vacuum_angular_wavenumber;
mod vacuum_wavenumber;
mod wavelength;

pub use angular_frequency::AngularFrequency;
pub use frequency::Frequency;
pub use inverse_length_squared::InverseLengthSquared;
pub use spectral::{ToCanonicalSpectral, ToCanonicalSpectralSquared};
pub use vacuum_angular_wavenumber::VacuumAngularWavenumber;
pub use vacuum_wavenumber::VacuumWavenumber;
pub use wavelength::Wavelength;

const SPEED_OF_LIGHT_CM_PER_SECOND: f64 = 2.997_924_58e10;
