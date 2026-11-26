//! Anisotropic models in the material's own (crystal) frame.
//!
//! Tensors are 3×3, row-major over the crystal axes a, b, c, and follow the
//! tensor capability traits of `lamina-core`. Orienting a model in a stack is
//! the job of a rotation applied after evaluation, `T = R T_crystal Rᵀ`; the
//! same rotation applies to the parameter rates reported here.

use crate::{ConstantTensorError, Sellmeier};
use lamina_core::{
    ComplexScalar, DerivativeOrder, DifferentiableMaterial, DifferentiableMeromorphicMaterial,
    DifferentiableMeromorphicTensorMaterial, DifferentiableTensorMaterial, Material,
    MeromorphicMaterial, MeromorphicTensorMaterial, Sampled, Scalar, SpectralSingularity,
    TensorMaterial, material::Tensor3,
};
use num_complex::Complex;
use num_traits::{Float, One, Zero};
use std::fmt::Debug;

/// A crystal axis: the rows and columns of a crystal-frame tensor.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CrystalAxis {
    A,
    B,
    C,
}

impl CrystalAxis {
    const ALL: [Self; 3] = [Self::A, Self::B, Self::C];
}

/// The scalar model of a [`PrincipalAxes`] medium that a parameter belongs to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PrincipalModel {
    /// The ordinary model of a uniaxial medium, on axes a and b.
    Ordinary,
    /// The extraordinary model of a uniaxial medium, on the optic axis c.
    Extraordinary,
    /// The model on one axis of a biaxial medium.
    Axis(CrystalAxis),
}

/// The parameter a [`TensorRate`] differentiates with respect to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TensorParameter {
    /// Parameter `name` of the scalar model `model`; `index` is the term for
    /// per-term parameters.
    Principal {
        model: PrincipalModel,
        name: &'static str,
        index: Option<usize>,
    },
    /// Entry `(row, column)` of a constant permittivity tensor.
    Permittivity {
        row: CrystalAxis,
        column: CrystalAxis,
    },
}

/// Rate `∂ε/∂p` of a crystal-frame permittivity tensor with respect to one
/// model parameter `p`, in the units of the scalar rate it comes from.
#[derive(Clone, Debug, PartialEq)]
pub struct TensorRate<R> {
    /// The parameter `p`.
    pub parameter: TensorParameter,
    /// `∂ε/∂p` in the crystal frame.
    pub value: Tensor3<Complex<R>>,
}

fn diagonal<C: Copy + Zero>([a, b, c]: [C; 3]) -> Tensor3<C> {
    let zero = C::zero();
    [[a, zero, zero], [zero, b, zero], [zero, zero, c]]
}

/// A medium whose permittivity and permeability are diagonal in its crystal
/// frame, with one scalar dispersion model per principal axis:
///
/// `ε = diag(ε_a, ε_b, ε_c)`, `μ = diag(μ_a, μ_b, μ_c)`.
///
/// A uniaxial medium has one ordinary model shared by a and b; a biaxial
/// medium has three independent models. Each capability trait is implemented
/// when all axis models have it, so a tabulated axis restricts the whole model
/// to the real axis. The model declares [`TensorMaterial::PRINCIPAL_AXES`].
#[derive(Clone, Debug, PartialEq)]
pub struct PrincipalAxes<X, Y, Z> {
    axes: Axes<X, Y, Z>,
}

#[derive(Clone, Debug, PartialEq)]
enum Axes<X, Y, Z> {
    Uniaxial { ordinary: X, extraordinary: Z },
    Biaxial { a: X, b: Y, c: Z },
}

/// A principal-axis model with a Sellmeier relation on every axis.
pub type SellmeierAxes<R> = PrincipalAxes<Sellmeier<R>, Sellmeier<R>, Sellmeier<R>>;

impl<X, Y, Z> PrincipalAxes<X, Y, Z> {
    /// A biaxial medium with independent models on the a, b and c axes.
    pub fn biaxial(a: X, b: Y, c: Z) -> Self {
        Self {
            axes: Axes::Biaxial { a, b, c },
        }
    }
}

impl<O, E> PrincipalAxes<O, O, E> {
    /// A uniaxial medium: the ordinary model on a and b, the extraordinary
    /// model on the optic axis c.
    pub fn uniaxial(ordinary: O, extraordinary: E) -> Self {
        Self {
            axes: Axes::Uniaxial {
                ordinary,
                extraordinary,
            },
        }
    }
}

impl<R: Float + Debug> SellmeierAxes<R> {
    /// Rates of the permittivity tensor with respect to every parameter of
    /// every scalar model, at real `k0` in `cm⁻¹`.
    ///
    /// Each rate is the model's [`Sellmeier::permittivity_parameter_rates`]
    /// entry on the diagonal elements of the axes it describes, zero
    /// elsewhere: `diag(∂ε_o, ∂ε_o, 0)` for an ordinary parameter,
    /// `diag(0, 0, ∂ε_e)` for an extraordinary one, and one diagonal element
    /// for a biaxial axis. Names are `constant`, `coefficient` and
    /// `pole_wavenumber`.
    pub fn permittivity_parameter_rates(&self, k0: R) -> Vec<TensorRate<R>> {
        let models = match &self.axes {
            Axes::Uniaxial {
                ordinary,
                extraordinary,
            } => vec![
                (ordinary, PrincipalModel::Ordinary, [true, true, false]),
                (
                    extraordinary,
                    PrincipalModel::Extraordinary,
                    [false, false, true],
                ),
            ],
            Axes::Biaxial { a, b, c } => vec![
                (
                    a,
                    PrincipalModel::Axis(CrystalAxis::A),
                    [true, false, false],
                ),
                (
                    b,
                    PrincipalModel::Axis(CrystalAxis::B),
                    [false, true, false],
                ),
                (
                    c,
                    PrincipalModel::Axis(CrystalAxis::C),
                    [false, false, true],
                ),
            ],
        };
        models
            .into_iter()
            .flat_map(|(sellmeier, model, axes)| {
                let rate = move |name, index, value| TensorRate {
                    parameter: TensorParameter::Principal { model, name, index },
                    value: diagonal(axes.map(|on| {
                        if on {
                            Complex::from(value)
                        } else {
                            Complex::zero()
                        }
                    })),
                };
                let rates = sellmeier.permittivity_parameter_rates(k0);
                let per_term = |name, values: Vec<R>| {
                    values
                        .into_iter()
                        .enumerate()
                        .map(move |(index, value)| rate(name, Some(index), value))
                };
                std::iter::once(rate("constant", None, rates.constant))
                    .chain(per_term("coefficient", rates.coefficients))
                    .chain(per_term("pole_wavenumber", rates.pole_wavenumbers))
            })
            .collect()
    }
}

macro_rules! principal_axes {
    ($tensor:ident for $scalar:ident, $k0:ty { $($item:item)* } $($method:ident => $axis_method:ident $(($order:ident))?;)+) => {
        impl<X, Y, Z> $tensor for PrincipalAxes<X, Y, Z>
        where
            X: $scalar,
            X::Real: Copy,
            Y: $scalar<Real = X::Real>,
            Z: $scalar<Real = X::Real>,
        {
            $($item)*
            $(
                fn $method<I, C>(&self, k0: I $(, $order: DerivativeOrder)?) -> I::Mapped<Tensor3<C>>
                where
                    C: ComplexScalar<RealField = Self::Real>,
                    I: Sampled<Elem = $k0>,
                {
                    k0.map(|k0| {
                        diagonal(match &self.axes {
                            Axes::Uniaxial { ordinary, extraordinary } => {
                                let ordinary = ordinary.$axis_method(Scalar(k0) $(, $order)?);
                                [ordinary, ordinary, extraordinary.$axis_method(Scalar(k0) $(, $order)?)]
                            }
                            Axes::Biaxial { a, b, c } => [
                                a.$axis_method(Scalar(k0) $(, $order)?),
                                b.$axis_method(Scalar(k0) $(, $order)?),
                                c.$axis_method(Scalar(k0) $(, $order)?),
                            ],
                        })
                    })
                }
            )+
        }
    };
}

principal_axes!(TensorMaterial for Material, Self::Real {
    type Real = X::Real;
    const PRINCIPAL_AXES: bool = true;
}
    relative_permittivity_tensor => relative_permittivity;
    relative_permeability_tensor => relative_permeability;
);
principal_axes!(DifferentiableTensorMaterial for DifferentiableMaterial, Self::Real {}
    relative_permittivity_tensor_derivative => relative_permittivity_derivative(order);
    relative_permeability_tensor_derivative => relative_permeability_derivative(order);
);
principal_axes!(MeromorphicTensorMaterial for MeromorphicMaterial, C {
    /// Each axis model's declaration on its entry, the ordinary model's of a
    /// uniaxial medium on a and b; nothing unless every axis model declares.
    fn tensor_spectral_singularities<C>(&self) -> Option<Vec<SpectralSingularity<C>>>
    where
        C: ComplexScalar<RealField = Self::Real>,
    {
        let on = |points: Vec<SpectralSingularity<C>>, axes: &[usize]| {
            axes.iter()
                .flat_map(|&axis| points.iter().map(move |p| p.on_axis(axis)))
                .collect::<Vec<_>>()
        };
        Some(match &self.axes {
            Axes::Uniaxial { ordinary, extraordinary } => {
                let mut points = on(ordinary.spectral_singularities()?, &[0, 1]);
                points.extend(on(extraordinary.spectral_singularities()?, &[2]));
                points
            }
            Axes::Biaxial { a, b, c } => {
                let mut points = on(a.spectral_singularities()?, &[0]);
                points.extend(on(b.spectral_singularities()?, &[1]));
                points.extend(on(c.spectral_singularities()?, &[2]));
                points
            }
        })
    }
}
    relative_permittivity_tensor_complex => relative_permittivity_complex;
    relative_permeability_tensor_complex => relative_permeability_complex;
);
principal_axes!(DifferentiableMeromorphicTensorMaterial for DifferentiableMeromorphicMaterial, C {}
    relative_permittivity_tensor_complex_derivative => relative_permittivity_complex_derivative(order);
    relative_permeability_tensor_complex_derivative => relative_permeability_complex_derivative(order);
);

/// A nondispersive medium with constant complex permittivity and permeability
/// tensors in its crystal frame.
///
/// The model declares no structure: it is a general tensor whatever its
/// values, and zero entries never make it diagonal.
#[derive(Clone, Debug, PartialEq)]
pub struct ConstantTensor<R> {
    permittivity: Tensor3<Complex<R>>,
    permeability: Tensor3<Complex<R>>,
}

impl<R: Float> ConstantTensor<R> {
    /// Construct from crystal-frame `ε` and `μ`.
    ///
    /// Every entry must be finite and the cc components nonzero.
    pub fn new(
        permittivity: Tensor3<Complex<R>>,
        permeability: Tensor3<Complex<R>>,
    ) -> Result<Self, ConstantTensorError> {
        for (quantity, tensor) in [
            ("permittivity", &permittivity),
            ("permeability", &permeability),
        ] {
            for (row, entries) in CrystalAxis::ALL.into_iter().zip(tensor) {
                for (column, entry) in CrystalAxis::ALL.into_iter().zip(entries) {
                    if !(entry.re.is_finite() && entry.im.is_finite()) {
                        return Err(ConstantTensorError::NonFiniteEntry {
                            quantity,
                            row,
                            column,
                        });
                    }
                }
            }
            if tensor[2][2].is_zero() {
                return Err(ConstantTensorError::ZeroCcComponent { quantity });
            }
        }
        Ok(Self {
            permittivity,
            permeability,
        })
    }

    /// Rates of the permittivity tensor with respect to each of its entries:
    /// `∂ε/∂ε_ij` is the unit tensor at `(i, j)`, at every wavenumber.
    pub fn permittivity_parameter_rates(&self) -> Vec<TensorRate<R>> {
        CrystalAxis::ALL
            .into_iter()
            .flat_map(|row| {
                CrystalAxis::ALL.map(|column| {
                    let mut value = [[Complex::zero(); 3]; 3];
                    value[row as usize][column as usize] = Complex::one();
                    TensorRate {
                        parameter: TensorParameter::Permittivity { row, column },
                        value,
                    }
                })
            })
            .collect()
    }
}

fn converted<R: Copy, C: ComplexScalar<RealField = R>>(tensor: &Tensor3<Complex<R>>) -> Tensor3<C> {
    tensor.map(|row| row.map(|entry| C::from_parts(entry.re, entry.im)))
}

impl<R: Float + Debug> TensorMaterial for ConstantTensor<R> {
    type Real = R;

    fn relative_permittivity_tensor<I, C>(&self, k0: I) -> I::Mapped<Tensor3<C>>
    where
        C: ComplexScalar<RealField = R>,
        I: Sampled<Elem = R>,
    {
        k0.map(|_| converted(&self.permittivity))
    }

    fn relative_permeability_tensor<I, C>(&self, k0: I) -> I::Mapped<Tensor3<C>>
    where
        C: ComplexScalar<RealField = R>,
        I: Sampled<Elem = R>,
    {
        k0.map(|_| converted(&self.permeability))
    }
}

impl<R: Float + Debug> DifferentiableTensorMaterial for ConstantTensor<R> {
    fn relative_permittivity_tensor_derivative<I, C>(
        &self,
        k0: I,
        _order: DerivativeOrder,
    ) -> I::Mapped<Tensor3<C>>
    where
        C: ComplexScalar<RealField = R>,
        I: Sampled<Elem = R>,
    {
        k0.map(|_| [[C::zero(); 3]; 3])
    }
}

impl<R: Float + Debug> MeromorphicTensorMaterial for ConstantTensor<R> {
    fn relative_permittivity_tensor_complex<I, C>(&self, k0: I) -> I::Mapped<Tensor3<C>>
    where
        C: ComplexScalar<RealField = R>,
        I: Sampled<Elem = C>,
    {
        k0.map(|_| converted(&self.permittivity))
    }

    fn relative_permeability_tensor_complex<I, C>(&self, k0: I) -> I::Mapped<Tensor3<C>>
    where
        C: ComplexScalar<RealField = R>,
        I: Sampled<Elem = C>,
    {
        k0.map(|_| converted(&self.permeability))
    }
}

impl<R: Float + Debug> DifferentiableMeromorphicTensorMaterial for ConstantTensor<R> {
    fn relative_permittivity_tensor_complex_derivative<I, C>(
        &self,
        k0: I,
        _order: DerivativeOrder,
    ) -> I::Mapped<Tensor3<C>>
    where
        C: ComplexScalar<RealField = R>,
        I: Sampled<Elem = C>,
    {
        k0.map(|_| [[C::zero(); 3]; 3])
    }
}
