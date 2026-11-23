//! Representative single-phonon polar dielectric records.

use crate::{Lorentz, LorentzOscillator, MaterialRecord, SpectralDomain};

use lamina_units::{Frequency, FrequencyUnit, VacuumWavenumber, Wavelength};

use num_traits::{Float, FloatConst, FromPrimitive};

use std::{fmt::Debug, marker::PhantomData};

fn convert<R>(value: f64) -> R
where
    R: FromPrimitive,
{
    R::from_f64(value).expect("f64 literature parameter converts to material scalar")
}

fn single_phonon_from_wavenumber<R>(
    epsilon_infinity: f64,
    transverse: f64,
    longitudinal: f64,
    damping: f64,
) -> Lorentz<R>
where
    R: Float + FloatConst + FromPrimitive + Debug,
{
    let ratio = longitudinal / transverse;

    let delta_epsilon = epsilon_infinity * (ratio * ratio - 1.0);

    let oscillator = LorentzOscillator::from_delta_epsilon(
        convert(delta_epsilon),
        VacuumWavenumber::inverse_centimetres(convert(transverse)),
        VacuumWavenumber::inverse_centimetres(convert(damping)),
    )
    .expect("single-phonon parameters are valid");

    Lorentz::new(convert(epsilon_infinity), vec![oscillator])
        .expect("single-phonon material is valid")
}

fn single_phonon_from_frequency<R>(
    epsilon_infinity: f64,
    transverse_hz: f64,
    longitudinal_hz: f64,
    damping_hz: f64,
) -> Lorentz<R>
where
    R: Float + FloatConst + FromPrimitive + Debug,
{
    let ratio = longitudinal_hz / transverse_hz;

    let delta_epsilon = epsilon_infinity * (ratio * ratio - 1.0);

    let frequency = |value: f64| Frequency::new(convert(value), FrequencyUnit::Hertz);

    let oscillator = LorentzOscillator::from_delta_epsilon(
        convert(delta_epsilon),
        frequency(transverse_hz),
        frequency(damping_hz),
    )
    .expect("single-phonon parameters are valid");

    Lorentz::new(convert(epsilon_infinity), vec![oscillator])
        .expect("single-phonon material is valid")
}

/// Silicon carbide single-phonon model used by Wang et al. (2012).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SiliconCarbideWang2012<R> {
    marker: PhantomData<R>,
}

impl<R> SiliconCarbideWang2012<R> {
    pub fn new() -> Self {
        Self {
            marker: PhantomData,
        }
    }
}

impl<R> SiliconCarbideWang2012<R>
where
    R: Float + FloatConst + FromPrimitive + Debug,
{
    pub fn model(&self) -> Lorentz<R> {
        single_phonon_from_frequency(6.7, 23.79e12, 29.07e12, 0.1428e12)
    }
}

impl<R> MaterialRecord for SiliconCarbideWang2012<R>
where
    R: Float + FloatConst + FromPrimitive,
{
    type Real = R;

    fn name(&self) -> &'static str {
        "Silicon carbide — Wang et al. (2012)"
    }

    fn identifier(&self) -> &'static str {
        "silicon-carbide-wang-2012"
    }

    fn source(&self) -> &'static str {
        "T. Wang et al., Optics Express 20, 13040–13050 (2012)"
    }

    fn spectral_domain(&self) -> Option<SpectralDomain<R>> {
        // Reststrahlen region quoted by the source.
        SpectralDomain::from_spectral(
            Wavelength::micrometres(convert(10.3)),
            Wavelength::micrometres(convert(12.6)),
        )
    }

    fn notes(&self) -> &'static str {
        "Single-phonon SiC model. The recorded domain is the quoted \
         Reststrahlen interval, not a claim of global fit validity."
    }
}

pub fn silicon_carbide_wang_2012<R>() -> SiliconCarbideWang2012<R> {
    SiliconCarbideWang2012::new()
}

/// Ordinary (`E ⟂ c`) AlN phonon response.
///
/// This is one tensor component of wurtzite AlN, represented as a scalar
/// material until anisotropic material support is available.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct AluminiumNitrideOrdinarySchoche2015<R> {
    marker: PhantomData<R>,
}

impl<R> AluminiumNitrideOrdinarySchoche2015<R> {
    pub fn new() -> Self {
        Self {
            marker: PhantomData,
        }
    }
}

impl<R> AluminiumNitrideOrdinarySchoche2015<R>
where
    R: Float + FloatConst + FromPrimitive + Debug,
{
    pub fn model(&self) -> Lorentz<R> {
        // Pure-AlN endpoint of the E1 response.
        single_phonon_from_wavenumber(4.08, 666.4, 909.9, 3.1)
    }
}

impl<R> MaterialRecord for AluminiumNitrideOrdinarySchoche2015<R>
where
    R: Float + FloatConst + FromPrimitive,
{
    type Real = R;

    fn name(&self) -> &'static str {
        "AlN ordinary — Schöche et al. (2015)"
    }

    fn identifier(&self) -> &'static str {
        "aln-ordinary-schoche-2015"
    }

    fn source(&self) -> &'static str {
        "S. Schöche et al., \
         \"Infrared dielectric function, phonon modes and free-charge carrier \
         properties of high-Al-content AlGaN alloys\", arXiv:1512.01427"
    }

    fn spectral_domain(&self) -> Option<SpectralDomain<R>> {
        SpectralDomain::from_spectral(
            VacuumWavenumber::inverse_centimetres(convert(300.0)),
            VacuumWavenumber::inverse_centimetres(convert(6000.0)),
        )
    }

    fn notes(&self) -> &'static str {
        "Ordinary E1 (E perpendicular to c) component only. \
         Wurtzite AlN is anisotropic; do not interpret this as the full \
         dielectric tensor."
    }
}

pub fn aluminium_nitride_ordinary_schoche_2015<R>() -> AluminiumNitrideOrdinarySchoche2015<R> {
    AluminiumNitrideOrdinarySchoche2015::new()
}

/// Ordinary (`E ⟂ c`) bulk GaN phonon response.
///
/// This record intentionally represents only the E1 component of wurtzite
/// GaN. It is suitable for scalar calculations where that component is the
/// physically relevant one.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct GalliumNitrideOrdinary<R> {
    marker: PhantomData<R>,
}

impl<R> GalliumNitrideOrdinary<R> {
    pub fn new() -> Self {
        Self {
            marker: PhantomData,
        }
    }
}

impl<R> GalliumNitrideOrdinary<R>
where
    R: Float + FloatConst + FromPrimitive + Debug,
{
    pub fn model(&self) -> Lorentz<R> {
        // Representative bulk ordinary E1 response:
        //
        //   ε∞  = 5.35
        //   TO  = 561 cm^-1
        //   LO  = 743 cm^-1
        //   γ   ≈ 5 cm^-1
        //
        // The TO/LO values are the bulk values used by Ratchford et al.;
        // the paper reports ~5 cm^-1 linewidths for bulk GaN films.
        single_phonon_from_wavenumber(5.35, 561.0, 743.0, 5.0)
    }
}

impl<R> MaterialRecord for GalliumNitrideOrdinary<R>
where
    R: Float + FloatConst + FromPrimitive,
{
    type Real = R;

    fn name(&self) -> &'static str {
        "GaN ordinary — bulk IR reference"
    }

    fn identifier(&self) -> &'static str {
        "gan-ordinary-bulk-ir"
    }

    fn source(&self) -> &'static str {
        "D. C. Ratchford et al., ACS Nano 13, 6730–6741 (2019); \
         A. Kasic et al., Phys. Rev. B 62, 7365–7377 (2000)"
    }

    fn spectral_domain(&self) -> Option<SpectralDomain<R>> {
        SpectralDomain::from_spectral(
            VacuumWavenumber::inverse_centimetres(convert(300.0)),
            VacuumWavenumber::inverse_centimetres(convert(1200.0)),
        )
    }

    fn notes(&self) -> &'static str {
        "Ordinary E1 (E perpendicular to c) scalar approximation only. \
         TO/LO frequencies follow the bulk values used by Ratchford et al.; \
         the high-frequency permittivity is representative of the \
         room-temperature GaN values reported by Kasic et al."
    }
}

pub fn gallium_nitride_ordinary<R>() -> GalliumNitrideOrdinary<R> {
    GalliumNitrideOrdinary::new()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polar_records_construct() {
        let _ = silicon_carbide_wang_2012::<f64>().model();

        let _ = aluminium_nitride_ordinary_schoche_2015::<f64>().model();

        let _ = gallium_nitride_ordinary::<f64>().model();
    }

    #[test]
    fn ordinary_nitride_models_have_one_oscillator() {
        assert_eq!(
            aluminium_nitride_ordinary_schoche_2015::<f64>()
                .model()
                .oscillators()
                .len(),
            1,
        );

        assert_eq!(
            gallium_nitride_ordinary::<f64>()
                .model()
                .oscillators()
                .len(),
            1,
        );
    }

    #[test]
    fn lo_to_mapping_produces_expected_static_permittivity() {
        let model = gallium_nitride_ordinary::<f64>().model();

        let expected = 5.35 * (743.0_f64 / 561.0).powi(2);

        let actual = model.epsilon_infinity() + model.oscillators()[0].delta_epsilon();

        approx::assert_relative_eq!(actual, expected, max_relative = 1e-14,);
    }
}
