//! Lorentz–Drude metal parameterisations from Rakić et al. (1998).
//!
//! The source model is
//!
//! ε(ω) = 1
//!      - f₀ ωₚ² / (ω² + i Γ₀ ω)
//!      + Σⱼ fⱼ ωₚ² / (ωⱼ² - ω² - i Γⱼ ω).
//!
//! Rakić et al. tabulate ℏωₚ, ℏΓⱼ and ℏωⱼ in electron volts.
//! These are converted to ordinary frequency using E = h f and then
//! canonicalised by `lamina-units`.

use crate::{DrudeLorentz, LorentzOscillator, MaterialRecord, SpectralDomain};

use lamina_units::{Frequency, FrequencyUnit};
use num_traits::{Float, FloatConst, FromPrimitive};
use std::{fmt::Debug, marker::PhantomData};

const SOURCE: &str = "A. D. Rakić, A. B. Djurišić, J. M. Elazar, and M. L. Majewski, \
     Applied Optics 37, 5271–5283 (1998)";

const LD_NOTES: &str = "Lorentz–Drude fit from Rakić et al. (1998), fitted over 0.1–6 eV.";

const NOBLE_LD_NOTES: &str = "Lorentz–Drude fit from Rakić et al. (1998), fitted over 0.1–6 eV. \
     The authors report that their Brendel–Bormann model gives a better \
     description of interband absorption for noble metals.";

const ELECTRON_VOLT_JOULE: f64 = 1.602_176_634e-19;
const PLANCK_JOULE_SECOND: f64 = 6.626_070_15e-34;

#[derive(Clone, Copy)]
struct OscillatorParameters {
    f: f64,
    omega_ev: f64,
    gamma_ev: f64,
}

#[derive(Clone, Copy)]
struct RakicLdParameters {
    plasma_ev: f64,
    f0: f64,
    gamma0_ev: f64,
    oscillators: &'static [OscillatorParameters],
}

fn ev<R>(value: f64) -> Frequency<R>
where
    R: Float + FromPrimitive,
{
    let hz_per_ev = ELECTRON_VOLT_JOULE / PLANCK_JOULE_SECOND;

    Frequency::new(
        R::from_f64(value * hz_per_ev).expect("f64 material parameter converts to model scalar"),
        FrequencyUnit::Hertz,
    )
}

fn build_model<R>(parameters: RakicLdParameters) -> DrudeLorentz<R>
where
    R: Float + FloatConst + FromPrimitive + Debug,
{
    let plasma_ev = parameters.plasma_ev;

    let oscillators = parameters
        .oscillators
        .iter()
        .map(|oscillator| {
            // Rakić:
            //
            // S_j = f_j ω_p²
            //
            // Our dimensionless zero-frequency oscillator contribution is
            //
            // Δε_j = S_j / ω_j²
            //      = f_j (ω_p / ω_j)².
            let delta_epsilon = oscillator.f * (plasma_ev / oscillator.omega_ev).powi(2);

            LorentzOscillator::from_delta_epsilon(
                R::from_f64(delta_epsilon)
                    .expect("f64 material parameter converts to model scalar"),
                ev(oscillator.omega_ev),
                ev(oscillator.gamma_ev),
            )
            .expect("validated literature parameterisation")
        })
        .collect();

    // Rakić's free-electron term is f₀ ωₚ², so its equivalent
    // plasma frequency is sqrt(f₀) ωₚ.
    let drude_plasma_ev = parameters.f0.sqrt() * parameters.plasma_ev;

    DrudeLorentz::from_plasma_wavenumber(
        R::one(),
        ev(drude_plasma_ev),
        ev(parameters.gamma0_ev),
        oscillators,
    )
    .expect("validated literature parameterisation")
}

fn spectral_domain<R>() -> Option<SpectralDomain<R>>
where
    R: Float + FloatConst + FromPrimitive,
{
    SpectralDomain::from_spectral(ev(0.1), ev(6.0))
}

macro_rules! rakic_ld_record {
    (
        $(#[$meta:meta])*
        $type:ident,
        $function:ident,
        name = $name:literal,
        id = $id:literal,
        notes = $notes:expr,
        plasma = $plasma:literal,
        f0 = $f0:literal,
        gamma0 = $gamma0:literal,
        oscillators = [
            $(
                ($f:literal, $omega:literal, $gamma:literal)
            ),* $(,)?
        ]
    ) => {
        $(#[$meta])*
        #[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
        pub struct $type<R> {
            marker: PhantomData<R>,
        }

        impl<R> $type<R> {
            pub fn new() -> Self {
                Self {
                    marker: PhantomData,
                }
            }
        }

        impl<R> $type<R>
        where
            R: Float + FloatConst + FromPrimitive + Debug,
        {
            /// Construct the Rakić et al. Lorentz–Drude material model.
            pub fn model(&self) -> DrudeLorentz<R> {
                const OSCILLATORS: &[OscillatorParameters] = &[
                    $(
                        OscillatorParameters {
                            f: $f,
                            omega_ev: $omega,
                            gamma_ev: $gamma,
                        },
                    )*
                ];

                build_model(RakicLdParameters {
                    plasma_ev: $plasma,
                    f0: $f0,
                    gamma0_ev: $gamma0,
                    oscillators: OSCILLATORS,
                })
            }
        }

        impl<R> MaterialRecord for $type<R>
        where
            R: Float + FloatConst + FromPrimitive,
        {
            type Real = R;

            fn name(&self) -> &'static str {
                $name
            }

            fn identifier(&self) -> &'static str {
                $id
            }

            fn source(&self) -> &'static str {
                SOURCE
            }

            fn spectral_domain(
                &self,
            ) -> Option<SpectralDomain<R>> {
                spectral_domain()
            }

            fn notes(&self) -> &'static str {
                $notes
            }
        }

        pub fn $function<R>() -> $type<R> {
            $type::new()
        }
    };
}

rakic_ld_record!(
    /// Silver using the Rakić et al. (1998) Lorentz–Drude fit.
    SilverRakic1998Ld,
    silver_rakic_1998_ld,
    name = "Silver — Rakić LD (1998)",
    id = "silver-rakic-1998-ld",
    notes = NOBLE_LD_NOTES,
    plasma = 9.01,
    f0 = 0.845,
    gamma0 = 0.048,
    oscillators = [
        (0.065, 0.816, 3.886),
        (0.124, 4.481, 0.452),
        (0.011, 8.185, 0.065),
        (0.840, 9.083, 0.916),
        (5.646, 20.29, 2.419),
    ]
);

rakic_ld_record!(
    /// Gold using the Rakić et al. (1998) Lorentz–Drude fit.
    GoldRakic1998Ld,
    gold_rakic_1998_ld,
    name = "Gold — Rakić LD (1998)",
    id = "gold-rakic-1998-ld",
    notes = NOBLE_LD_NOTES,
    plasma = 9.03,
    f0 = 0.760,
    gamma0 = 0.053,
    oscillators = [
        (0.024, 0.415, 0.241),
        (0.010, 0.830, 0.345),
        (0.071, 2.969, 0.870),
        (0.601, 4.304, 2.494),
        (4.384, 13.32, 2.214),
    ]
);

rakic_ld_record!(
    /// Copper using the Rakić et al. (1998) Lorentz–Drude fit.
    CopperRakic1998Ld,
    copper_rakic_1998_ld,
    name = "Copper — Rakić LD (1998)",
    id = "copper-rakic-1998-ld",
    notes = NOBLE_LD_NOTES,
    plasma = 10.83,
    f0 = 0.575,
    gamma0 = 0.030,
    oscillators = [
        (0.061, 0.291, 0.378),
        (0.104, 2.957, 1.056),
        (0.723, 5.300, 3.213),
        (0.638, 11.18, 4.305),
    ]
);

rakic_ld_record!(
    /// Aluminium using the Rakić et al. (1998) Lorentz–Drude fit.
    AluminiumRakic1998Ld,
    aluminium_rakic_1998_ld,
    name = "Aluminium — Rakić LD (1998)",
    id = "aluminium-rakic-1998-ld",
    notes = LD_NOTES,
    plasma = 14.98,
    f0 = 0.523,
    gamma0 = 0.047,
    oscillators = [
        (0.227, 0.162, 0.333),
        (0.050, 1.544, 0.312),
        (0.166, 1.808, 1.351),
        (0.030, 3.473, 3.382),
    ]
);

rakic_ld_record!(
    /// Beryllium using the Rakić et al. (1998) Lorentz–Drude fit.
    BerylliumRakic1998Ld,
    beryllium_rakic_1998_ld,
    name = "Beryllium — Rakić LD (1998)",
    id = "beryllium-rakic-1998-ld",
    notes = LD_NOTES,
    plasma = 18.51,
    f0 = 0.084,
    gamma0 = 0.035,
    oscillators = [
        (0.031, 0.100, 1.664),
        (0.140, 1.032, 3.395),
        (0.530, 3.183, 4.454),
        (0.130, 4.604, 1.802),
    ]
);

rakic_ld_record!(
    /// Chromium using the Rakić et al. (1998) Lorentz–Drude fit.
    ChromiumRakic1998Ld,
    chromium_rakic_1998_ld,
    name = "Chromium — Rakić LD (1998)",
    id = "chromium-rakic-1998-ld",
    notes = LD_NOTES,
    plasma = 10.75,
    f0 = 0.168,
    gamma0 = 0.047,
    oscillators = [
        (0.151, 0.121, 3.175),
        (0.150, 0.543, 1.305),
        (1.149, 1.970, 2.676),
        (0.825, 8.775, 1.335),
    ]
);

rakic_ld_record!(
    /// Nickel using the Rakić et al. (1998) Lorentz–Drude fit.
    NickelRakic1998Ld,
    nickel_rakic_1998_ld,
    name = "Nickel — Rakić LD (1998)",
    id = "nickel-rakic-1998-ld",
    notes = LD_NOTES,
    plasma = 15.92,
    f0 = 0.096,
    gamma0 = 0.048,
    oscillators = [
        (0.100, 0.174, 4.511),
        (0.135, 0.582, 1.334),
        (0.106, 1.597, 2.178),
        (0.729, 6.089, 6.292),
    ]
);

rakic_ld_record!(
    /// Palladium using the Rakić et al. (1998) Lorentz–Drude fit.
    PalladiumRakic1998Ld,
    palladium_rakic_1998_ld,
    name = "Palladium — Rakić LD (1998)",
    id = "palladium-rakic-1998-ld",
    notes = LD_NOTES,
    plasma = 9.72,
    f0 = 0.330,
    gamma0 = 0.008,
    oscillators = [
        (0.649, 0.336, 2.950),
        (0.121, 0.501, 0.555),
        (0.638, 1.659, 4.621),
        (0.453, 5.715, 3.236),
    ]
);

rakic_ld_record!(
    /// Platinum using the Rakić et al. (1998) Lorentz–Drude fit.
    PlatinumRakic1998Ld,
    platinum_rakic_1998_ld,
    name = "Platinum — Rakić LD (1998)",
    id = "platinum-rakic-1998-ld",
    notes = LD_NOTES,
    plasma = 9.59,
    f0 = 0.333,
    gamma0 = 0.080,
    oscillators = [
        (0.191, 0.780, 0.517),
        (0.659, 1.314, 1.838),
        (0.547, 3.141, 3.668),
        (3.576, 9.249, 8.517),
    ]
);

rakic_ld_record!(
    /// Titanium using the Rakić et al. (1998) Lorentz–Drude fit.
    TitaniumRakic1998Ld,
    titanium_rakic_1998_ld,
    name = "Titanium — Rakić LD (1998)",
    id = "titanium-rakic-1998-ld",
    notes = LD_NOTES,
    plasma = 7.29,
    f0 = 0.148,
    gamma0 = 0.082,
    oscillators = [
        (0.899, 0.777, 2.276),
        (0.393, 1.545, 2.518),
        (0.187, 2.509, 1.663),
        (0.001, 19.43, 1.762),
    ]
);

rakic_ld_record!(
    /// Tungsten using the Rakić et al. (1998) Lorentz–Drude fit.
    TungstenRakic1998Ld,
    tungsten_rakic_1998_ld,
    name = "Tungsten — Rakić LD (1998)",
    id = "tungsten-rakic-1998-ld",
    notes = LD_NOTES,
    plasma = 13.22,
    f0 = 0.206,
    gamma0 = 0.064,
    oscillators = [
        (0.054, 1.004, 0.530),
        (0.166, 1.917, 1.281),
        (0.706, 3.580, 3.332),
        (2.590, 7.498, 5.836),
    ]
);

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn all_rakic_ld_records_construct() {
        let _ = silver_rakic_1998_ld::<f64>().model();
        let _ = gold_rakic_1998_ld::<f64>().model();
        let _ = copper_rakic_1998_ld::<f64>().model();
        let _ = aluminium_rakic_1998_ld::<f64>().model();
        let _ = beryllium_rakic_1998_ld::<f64>().model();
        let _ = chromium_rakic_1998_ld::<f64>().model();
        let _ = nickel_rakic_1998_ld::<f64>().model();
        let _ = palladium_rakic_1998_ld::<f64>().model();
        let _ = platinum_rakic_1998_ld::<f64>().model();
        let _ = titanium_rakic_1998_ld::<f64>().model();
        let _ = tungsten_rakic_1998_ld::<f64>().model();
    }

    #[test]
    fn silver_drude_parameters_follow_rakic_mapping() {
        let model = silver_rakic_1998_ld::<f64>().model();

        let expected_plasma = ev::<f64>(9.01 * 0.845_f64.sqrt())
            .vacuum_angular_wavenumber()
            .as_inverse_centimetres();

        let expected_damping = ev::<f64>(0.048)
            .vacuum_angular_wavenumber()
            .as_inverse_centimetres();

        assert_relative_eq!(
            model.plasma_wavenumber(),
            expected_plasma,
            max_relative = 1e-14,
        );

        assert_relative_eq!(
            model.drude_damping(),
            expected_damping,
            max_relative = 1e-14,
        );
    }

    #[test]
    fn rakic_domain_is_reported_fit_interval() {
        let domain = gold_rakic_1998_ld::<f64>().spectral_domain().unwrap();

        let lower = ev::<f64>(0.1)
            .vacuum_angular_wavenumber()
            .as_inverse_centimetres();

        let upper = ev::<f64>(6.0)
            .vacuum_angular_wavenumber()
            .as_inverse_centimetres();

        assert_relative_eq!(domain.minimum(), lower, max_relative = 1e-14,);

        assert_relative_eq!(domain.maximum(), upper, max_relative = 1e-14,);
    }
}
