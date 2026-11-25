//! Construction and validation errors.

use crate::CrystalAxis;
use std::fmt::Debug;

/// Validation errors shared by analytical models.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
pub enum MaterialModelError<R>
where
    R: Debug,
{
    /// A parameter is NaN or infinite.
    #[error("material parameter `{name}` must be finite, found {value:?}")]
    NonFiniteParameter {
        /// Parameter name.
        name: &'static str,
        /// Invalid value.
        value: R,
    },

    /// A parameter required to be nonnegative is negative.
    #[error("material parameter `{name}` must be nonnegative, found {value:?}")]
    NegativeParameter {
        /// Parameter name.
        name: &'static str,
        /// Invalid value.
        value: R,
    },

    /// A parameter required to be positive is zero or negative.
    #[error("material parameter `{name}` must be positive, found {value:?}")]
    NonPositiveParameter {
        /// Parameter name.
        name: &'static str,
        /// Invalid value.
        value: R,
    },

    /// A coefficient list is empty.
    #[error("material parameter `{name}` must contain at least one coefficient")]
    EmptyCoefficients {
        /// Parameter name.
        name: &'static str,
    },
}

/// Validation errors for tabulated constitutive data.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
pub enum TabulatedMaterialError<R>
where
    R: Debug,
{
    /// Fewer than two samples were supplied.
    #[error("tabulated material requires at least two samples, found {len}")]
    TooFewSamples {
        /// Number of samples.
        len: usize,
    },

    /// Coordinate and value lengths differ.
    #[error("tabulated material has {coordinate_len} coordinates but {value_len} values")]
    LengthMismatch {
        /// Coordinate count.
        coordinate_len: usize,
        /// Value count.
        value_len: usize,
    },

    /// A coordinate is non-finite.
    #[error("tabulated coordinate at index {index} must be finite, found {value:?}")]
    NonFiniteCoordinate {
        /// Sample index.
        index: usize,
        /// Invalid value.
        value: R,
    },

    /// A complex value has a non-finite component.
    #[error("tabulated value at index {index} must be finite")]
    NonFiniteValue {
        /// Sample index.
        index: usize,
    },

    /// Coordinates are not strictly increasing.
    #[error(
        "tabulated coordinates must be strictly increasing: index {right_index}          has {right:?}, following {left:?}"
    )]
    NonIncreasingCoordinates {
        /// Previous coordinate.
        left: R,
        /// Invalid coordinate.
        right: R,
        /// Invalid coordinate index.
        right_index: usize,
    },
}

/// Validation errors for constant tensor models.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum ConstantTensorError {
    /// An entry has a non-finite component.
    #[error("{quantity} entry ({row:?}, {column:?}) must be finite")]
    NonFiniteEntry {
        /// `permittivity` or `permeability`.
        quantity: &'static str,
        /// Entry row.
        row: CrystalAxis,
        /// Entry column.
        column: CrystalAxis,
    },

    /// The cc component vanishes.
    #[error("{quantity} cc component must be nonzero")]
    ZeroCcComponent {
        /// `permittivity` or `permeability`.
        quantity: &'static str,
    },
}
