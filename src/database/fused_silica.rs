use crate::{LengthUnit, MaterialRecord, Sellmeier, SellmeierTerm, SpectralDomain, Wavelength};
use num_traits::{Float, FloatConst, FromPrimitive};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FusedSilica<R> {
    marker: std::marker::PhantomData<R>,
}

impl<R> FusedSilica<R>
where
    R: Float + FloatConst + FromPrimitive + std::fmt::Debug,
{
    pub fn new() -> Self {
        Self {
            marker: std::marker::PhantomData,
        }
    }
    pub fn model(&self) -> Sellmeier<R> {
        let term = |b: f64, pole_um: f64| {
            SellmeierTerm::from_wavelength_pole(
                R::from_f64(b).unwrap(),
                Wavelength::new(R::from_f64(pole_um).unwrap(), LengthUnit::Micrometre).unwrap(),
            )
            .unwrap()
        };
        Sellmeier::with_unit_constant(vec![
            term(0.696_166_3, 0.068_404_3),
            term(0.407_942_6, 0.116_241_4),
            term(0.897_479_4, 9.896_161),
        ])
        .unwrap()
    }
}
impl<R> MaterialRecord for FusedSilica<R>
where
    R: Float + FloatConst + FromPrimitive,
{
    type Real = R;
    fn name(&self) -> &'static str {
        "Fused silica"
    }
    fn identifier(&self) -> &'static str {
        "fused-silica-malitson-1965"
    }
    fn source(&self) -> &'static str {
        "I. H. Malitson, JOSA 55, 1205-1209 (1965)"
    }
    fn spectral_domain(&self) -> Option<SpectralDomain<R>> {
        let lmin = R::from_f64(0.21e-4)?;
        let lmax = R::from_f64(3.71e-4)?;
        let two_pi = (R::one() + R::one()) * R::PI();
        SpectralDomain::new(two_pi / lmax, two_pi / lmin)
    }
    fn notes(&self) -> &'static str {
        "Validity range converted from 0.21-3.71 micrometres."
    }
}
pub fn fused_silica<R>() -> FusedSilica<R>
where
    R: Float + FloatConst + FromPrimitive + std::fmt::Debug,
{
    FusedSilica::new()
}
