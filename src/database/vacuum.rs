use crate::{MaterialRecord, SpectralDomain};
use lamina_core::Constant;
use num_traits::Float;
use std::marker::PhantomData;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Vacuum<R> {
    marker: PhantomData<R>,
}
impl<R> Vacuum<R> {
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
impl<R> MaterialRecord for Vacuum<R> {
    type Real = R;
    fn name(&self) -> &'static str {
        "Vacuum"
    }
    fn identifier(&self) -> &'static str {
        "vacuum"
    }
    fn source(&self) -> &'static str {
        "Definition: epsilon_r = mu_r = 1"
    }
    fn spectral_domain(&self) -> Option<SpectralDomain<R>> {
        None
    }
}
pub fn vacuum<R>() -> Vacuum<R> {
    Vacuum::new()
}
