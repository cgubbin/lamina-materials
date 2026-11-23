//! Room-temperature silicon and germanium parameterisations from
//! Frey, Leviton, and Madison (2006).

use crate::{MaterialRecord, Sellmeier, SellmeierTerm, SpectralDomain};

use lamina_units::Wavelength;

use num_traits::{Float, FloatConst, FromPrimitive};

use std::{fmt::Debug, marker::PhantomData};

const SOURCE: &str = "B. J. Frey, D. B. Leviton, and T. J. Madison, \
     \"Temperature-dependent refractive index of silicon and germanium\", \
     Proc. SPIE 6273, 62732J (2006), arXiv:physics/0606168";

const TEMPERATURE_K: f64 = 295.0;

fn polynomial(coefficients: [f64; 5], temperature: f64) -> f64 {
    coefficients
        .into_iter()
        .rev()
        .fold(0.0, |value, coefficient| value * temperature + coefficient)
}

fn convert<R>(value: f64) -> R
where
    R: FromPrimitive,
{
    R::from_f64(value).expect("f64 literature parameter converts to material scalar")
}

const SI_S1: [f64; 5] = [10.4907, -2.08020e-4, 4.21694e-6, -5.82298e-9, 3.44688e-12];

const SI_S2: [f64; 5] = [-1346.61, 29.1664, -0.278724, 1.05939e-3, -1.35089e-6];

const SI_S3: [f64; 5] = [4.42827e7, -1.76213e6, -7.61575e4, 678.414, 103.243];

const SI_L1: [f64; 5] = [0.299713, -1.14234e-5, 1.67134e-7, -2.51049e-10, 2.32484e-14];

const SI_L2: [f64; 5] = [-3.51710e3, 42.3892, -0.357957, 1.17504e-3, -1.13212e-6];

const SI_L3: [f64; 5] = [1.71400e6, -1.44984e5, -6.90744e3, -39.3699, 23.5770];

const GE_S1: [f64; 5] = [13.9723, 2.52809e-3, -5.02195e-6, 2.22604e-8, -4.86238e-12];

const GE_S2: [f64; 5] = [0.452096, -3.09197e-3, 2.16895e-5, -6.02290e-8, 4.12038e-11];

const GE_S3: [f64; 5] = [751.447, -14.2843, -0.238093, 2.96047e-3, -7.73454e-6];

const GE_L1: [f64; 5] = [0.386367, 2.01871e-4, -5.93448e-7, -2.27923e-10, 5.37423e-12];

const GE_L2: [f64; 5] = [1.08843, 1.16510e-3, -4.97284e-6, 1.12357e-8, 9.40201e-12];

const GE_L3: [f64; 5] = [-2893.19, -0.967948, -0.527016, 6.49364e-3, -1.95162e-5];

fn sellmeier_term<R>(
    strength_coefficients: [f64; 5],
    pole_coefficients: [f64; 5],
) -> SellmeierTerm<R>
where
    R: Float + FloatConst + FromPrimitive + Debug,
{
    let strength = polynomial(strength_coefficients, TEMPERATURE_K);

    // Only λ_i² appears in the Sellmeier relation.
    let pole_um = polynomial(pole_coefficients, TEMPERATURE_K).abs();

    SellmeierTerm::from_wavelength_pole(
        convert(strength),
        Wavelength::micrometres(convert(pole_um)),
    )
    .expect("Frey Sellmeier term is valid at 295 K")
}

/// Silicon at 295 K using the Frey–Leviton–Madison CHARMS fit.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SiliconFrey2006<R> {
    marker: PhantomData<R>,
}

impl<R> SiliconFrey2006<R> {
    pub fn new() -> Self {
        Self {
            marker: PhantomData,
        }
    }
}

impl<R> SiliconFrey2006<R>
where
    R: Float + FloatConst + FromPrimitive + Debug,
{
    pub fn model(&self) -> Sellmeier<R> {
        Sellmeier::with_unit_constant(vec![
            sellmeier_term(SI_S1, SI_L1),
            sellmeier_term(SI_S2, SI_L2),
            sellmeier_term(SI_S3, SI_L3),
        ])
        .expect("Frey silicon parameterisation is valid")
    }
}

impl<R> MaterialRecord for SiliconFrey2006<R>
where
    R: Float + FloatConst + FromPrimitive,
{
    type Real = R;

    fn name(&self) -> &'static str {
        "Silicon — Frey et al. (2006), 295 K"
    }

    fn identifier(&self) -> &'static str {
        "silicon-frey-2006-295k"
    }

    fn source(&self) -> &'static str {
        SOURCE
    }

    fn spectral_domain(&self) -> Option<SpectralDomain<R>> {
        SpectralDomain::from_spectral(
            Wavelength::micrometres(convert(1.1)),
            Wavelength::micrometres(convert(5.6)),
        )
    }

    fn notes(&self) -> &'static str {
        "295 K specialization of the CHARMS temperature-dependent Sellmeier fit. \
         Source fit covers 20–300 K and 1.1–5.6 µm."
    }
}

pub fn silicon_frey_2006<R>() -> SiliconFrey2006<R> {
    SiliconFrey2006::new()
}

/// Germanium at 295 K using the Frey–Leviton–Madison CHARMS fit.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct GermaniumFrey2006<R> {
    marker: PhantomData<R>,
}

impl<R> GermaniumFrey2006<R> {
    pub fn new() -> Self {
        Self {
            marker: PhantomData,
        }
    }
}

impl<R> GermaniumFrey2006<R>
where
    R: Float + FloatConst + FromPrimitive + Debug,
{
    pub fn model(&self) -> Sellmeier<R> {
        Sellmeier::with_unit_constant(vec![
            sellmeier_term(GE_S1, GE_L1),
            sellmeier_term(GE_S2, GE_L2),
            sellmeier_term(GE_S3, GE_L3),
        ])
        .expect("Frey germanium parameterisation is valid")
    }
}

impl<R> MaterialRecord for GermaniumFrey2006<R>
where
    R: Float + FloatConst + FromPrimitive,
{
    type Real = R;

    fn name(&self) -> &'static str {
        "Germanium — Frey et al. (2006), 295 K"
    }

    fn identifier(&self) -> &'static str {
        "germanium-frey-2006-295k"
    }

    fn source(&self) -> &'static str {
        SOURCE
    }

    fn spectral_domain(&self) -> Option<SpectralDomain<R>> {
        SpectralDomain::from_spectral(
            Wavelength::micrometres(convert(1.9)),
            Wavelength::micrometres(convert(5.5)),
        )
    }

    fn notes(&self) -> &'static str {
        "295 K specialization of the CHARMS temperature-dependent Sellmeier fit. \
         Source fit covers 20–300 K and 1.9–5.5 µm."
    }
}

pub fn germanium_frey_2006<R>() -> GermaniumFrey2006<R> {
    GermaniumFrey2006::new()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn silicon_constructs() {
        let material = silicon_frey_2006::<f64>();

        assert_eq!(material.model().terms().len(), 3);
    }

    #[test]
    fn germanium_constructs() {
        let material = germanium_frey_2006::<f64>();

        assert_eq!(material.model().terms().len(), 3);
    }

    #[test]
    fn domains_are_available() {
        assert!(silicon_frey_2006::<f64>().spectral_domain().is_some());

        assert!(germanium_frey_2006::<f64>().spectral_domain().is_some());
    }
}
