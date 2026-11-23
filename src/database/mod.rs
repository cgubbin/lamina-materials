//! Named material parameterisations.
//!
//! This module provides small, documented parameter sets built from the
//! generic models in [`crate::model`].
//!
//! Database records carry literature identity, validity information, and
//! construction details; the returned models themselves remain ordinary
//! `lamina-materials` model types.
//!
//! Each literature-backed record should document:
//!
//! - the original reference;
//! - the parameter convention;
//! - the valid spectral interval;
//! - units and conversions;
//! - temperature, sample, or fabrication conditions where relevant.
//!
//! This is intentionally a curated library rather than a comprehensive
//! optical-constants database.

mod air;
mod frey_2006;
mod fused_silica;
mod polar;
mod rakic_1998;
mod vacuum;

pub use air::{Air, air};
pub use frey_2006::{GermaniumFrey2006, SiliconFrey2006, germanium_frey_2006, silicon_frey_2006};
pub use fused_silica::{FusedSilica, fused_silica};
pub use polar::{
    AluminiumNitrideOrdinarySchoche2015, GalliumNitrideOrdinary, SiliconCarbideWang2012,
    aluminium_nitride_ordinary_schoche_2015, gallium_nitride_ordinary, silicon_carbide_wang_2012,
};
pub use rakic_1998::{
    AluminiumRakic1998Ld, BerylliumRakic1998Ld, ChromiumRakic1998Ld, CopperRakic1998Ld,
    GoldRakic1998Ld, NickelRakic1998Ld, PalladiumRakic1998Ld, PlatinumRakic1998Ld,
    SilverRakic1998Ld, TitaniumRakic1998Ld, TungstenRakic1998Ld, aluminium_rakic_1998_ld,
    beryllium_rakic_1998_ld, chromium_rakic_1998_ld, copper_rakic_1998_ld, gold_rakic_1998_ld,
    nickel_rakic_1998_ld, palladium_rakic_1998_ld, platinum_rakic_1998_ld, silver_rakic_1998_ld,
    titanium_rakic_1998_ld, tungsten_rakic_1998_ld,
};
pub use vacuum::{Vacuum, vacuum};
