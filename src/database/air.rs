use crate::{MaterialRecord, SpectralDomain};
use lamina_core::Constant;
use num_traits::Float;
use std::marker::PhantomData;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Air<R> {
    marker: PhantomData<R>,
}
impl<R> Air<R> {
    pub fn new() -> Self {
        Self {
            marker: PhantomData,
        }
    }
    pub fn model(&self) -> Constant<R>
    where
        R: Float,
    {
        Constant::new(R::one(), R::one())
    }
}
impl<R> MaterialRecord for Air<R> {
    type Real = R;
    fn name(&self) -> &'static str {
        "Air (idealised)"
    }
    fn identifier(&self) -> &'static str {
        "air-ideal"
    }
    fn source(&self) -> &'static str {
        "Idealised approximation: epsilon_r = mu_r = 1"
    }
    fn spectral_domain(&self) -> Option<SpectralDomain<R>> {
        None
    }
    fn notes(&self) -> &'static str {
        "Neglects refractivity, humidity, pressure, and absorption."
    }
}
pub fn air<R>() -> Air<R> {
    Air::new()
}
