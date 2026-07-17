//! Documented material parameter sets.
//!
//! This module will contain named literature parameterisations of the generic
//! models in [`crate::model`].
//!
//! Every database entry should document:
//!
//! - the original reference;
//! - the parameter convention;
//! - the valid spectral interval;
//! - units and conversions;
//! - temperature or sample conditions where relevant.
//!
//! The initial crate release focuses on validating the generic analytical
//! models before adding named material records.
//! Named and documented material records.

mod air;
mod fused_silica;
mod vacuum;

pub use air::{Air, air};
pub use fused_silica::{FusedSilica, fused_silica};
pub use vacuum::{Vacuum, vacuum};
