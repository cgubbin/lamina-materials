//! A `PrincipalAxes` model with a tabulated axis placed in a lamina-core
//! `TensorStack` and evaluated on the real axis through both 2×2 backends, and
//! a scalar tabulated layer evaluated outside its table.

use std::{error::Error, f64::consts::FRAC_PI_2};

use lamina_core::{
    ConstitutiveClass, ConstitutiveError, CoordinateInput, Coordinates, FiniteLayerIndex,
    InPlaneCoordinate, IncidentSide, Material, MaterialRegion, MaterialStack, Polarisation,
    RealAxisEvaluator, Rotation, Scalar, Scatter2, SpectralCoordinate, Stack, StackRegion,
    TensorMaterial, TensorMaterialHandle, TensorStack, Transfer2,
    backend::{Scatter2Error, Transfer2Error},
    material::MaterialHandle,
    units::InverseLengthUnit,
};
use lamina_materials::{Extrapolation, PrincipalAxes, Sellmeier, SellmeierTerm, Tabulated};
use lamina_units::{Length, VacuumAngularWavenumber};
use ndarray::{Array1, Ix1, arr1};
use num_complex::Complex64;

const POLARISATIONS: [Polarisation; 2] = [
    Polarisation::TransverseElectric,
    Polarisation::TransverseMagnetic,
];

/// Table samples, `cm⁻¹`, between the Sellmeier poles.
const SAMPLES: [f64; 5] = [2.0e4, 2.4e4, 2.8e4, 3.2e4, 3.6e4];

/// In-plane wavenumber as a fraction of `k₀`: oblique and propagating in vacuum.
const IN_PLANE: f64 = 0.6;

/// Agreement at the table samples, relative to `max(1, |expected|)`.
const TOLERANCE: f64 = 1e-14;

/// Sellmeier constant and `(B, q)` terms, `q` in `cm⁻¹`.
type Parameters = (f64, [(f64, f64); 2]);

const ORDINARY: Parameters = (1.0, [(3.9, 1.3e5), (2.6, 2.6e3)]);
const EXTRAORDINARY: Parameters = (1.2, [(3.6, 1.2e5), (2.4, 2.7e3)]);
const CLADDING: Parameters = (1.1, [(1.0, 1.1e5), (0.4, 3.0e3)]);

fn sellmeier((constant, terms): Parameters) -> Sellmeier<f64> {
    let terms = terms
        .iter()
        .map(|&(coefficient, pole)| {
            SellmeierTerm::new(
                coefficient,
                VacuumAngularWavenumber::inverse_centimetres(pole),
            )
            .unwrap()
        })
        .collect();
    Sellmeier::new(constant, terms).unwrap()
}

/// The extraordinary Sellmeier curve sampled at `SAMPLES`.
fn table(extrapolation: Extrapolation) -> Tabulated<f64> {
    let curve = sellmeier(EXTRAORDINARY);
    let values = SAMPLES
        .iter()
        .map(|&k| curve.relative_permittivity(Scalar(k)))
        .collect::<Array1<Complex64>>();
    Tabulated::new(arr1(&SAMPLES), values, extrapolation).unwrap()
}

fn input(k0: &[f64]) -> CoordinateInput<f64, Ix1> {
    let in_plane: Vec<f64> = k0.iter().map(|k| IN_PLANE * k).collect();
    CoordinateInput::samples(
        Coordinates::new(
            SpectralCoordinate::VacuumAngularWavenumber(InverseLengthUnit::PerCentimetre),
            InPlaneCoordinate::ParallelAngularWavenumber(InverseLengthUnit::PerCentimetre),
        ),
        arr1(k0),
        arr1(&in_plane),
    )
    .unwrap()
}

/// Vacuum | the uniaxial layer, optic axis set by `rotation` | a dispersive
/// isotropic cladding exterior.
fn stack<M>(model: M, rotation: Rotation<f64>) -> TensorStack<Complex64>
where
    M: TensorMaterial<Real = f64> + Send + Sync + 'static,
{
    let scalar = |material| MaterialRegion::from(MaterialHandle::new(material));
    let layer = MaterialRegion::oriented(TensorMaterialHandle::new(model), rotation);
    assert_eq!(layer.class(), ConstitutiveClass::Diagonal);
    Stack::builder(
        scalar(Sellmeier::new(1.0, Vec::new()).unwrap()),
        scalar(sellmeier(CLADDING)),
    )
    .layer(layer, Length::micrometres(0.8))
    .finalise()
}

/// Optic axis along stack x, then along stack z, where only TM sees the
/// extraordinary axis; then along stack y, where only TE does.
fn rotations() -> [Rotation<f64>; 3] {
    [
        Rotation::from_optic_axis(FRAC_PI_2, 0.0),
        Rotation::identity(),
        Rotation::from_optic_axis(FRAC_PI_2, FRAC_PI_2),
    ]
}

/// `|actual − expected| ≤ TOLERANCE · max(1, |expected|)`, entrywise.
fn assert_close<'a>(
    actual: impl IntoIterator<Item = &'a Complex64>,
    expected: impl IntoIterator<Item = &'a Complex64>,
    label: &str,
) {
    for (index, (a, e)) in actual.into_iter().zip(expected).enumerate() {
        let error = (a - e).norm();
        assert!(
            error <= TOLERANCE * e.norm().max(1.0),
            "{label}[{index}]: expected {e}, got {a}, error {error:e}"
        );
    }
}

/// Run `$body` with `$evaluator` bound to each real-axis 2×2 backend.
macro_rules! for_each_backend {
    ($evaluator:ident, $body:block) => {{
        {
            let $evaluator = RealAxisEvaluator::new(Scatter2::new());
            $body
        }
        {
            let $evaluator = RealAxisEvaluator::new(Transfer2::new());
            $body
        }
    }};
}

/// Compare amplitudes and power of `actual` and `expected` at `k0`, both
/// polarisations, both incident sides, both backends.
fn assert_same_response(
    actual: &TensorStack<Complex64>,
    expected: &TensorStack<Complex64>,
    k0: &[f64],
) {
    for polarisation in POLARISATIONS {
        for_each_backend!(evaluator, {
            let a = evaluator.evaluate(input(k0), actual, polarisation).unwrap();
            let b = evaluator
                .evaluate(input(k0), expected, polarisation)
                .unwrap();
            for side in [IncidentSide::Left, IncidentSide::Right] {
                let label = format!("{polarisation:?} {side:?}");
                let (p, q) = (a.amplitudes(side).unwrap(), b.amplitudes(side).unwrap());
                assert_close(p.reflection(), q.reflection(), &format!("{label} r"));
                assert_close(p.transmission(), q.transmission(), &format!("{label} t"));
                let (p, q) = (a.power(side).unwrap(), b.power(side).unwrap());
                for (name, x, y) in [
                    ("R", p.reflectance(), q.reflectance()),
                    ("T", p.transmittance(), q.transmittance()),
                ] {
                    let (x, y) = (x.mapv(Complex64::from), y.mapv(Complex64::from));
                    assert!(
                        x.iter().all(|v| v.re.is_finite() && v.re > 0.0),
                        "{label} {name}"
                    );
                    assert_close(&x, &y, &format!("{label} {name}"));
                }
            }
        });
    }
}

#[test]
fn a_tabulated_principal_axis_in_a_stack_matches_its_sellmeier_at_the_samples() {
    for rotation in rotations() {
        let tabulated = stack(
            PrincipalAxes::uniaxial(sellmeier(ORDINARY), table(Extrapolation::Clamp)),
            rotation,
        );
        let analytical = stack(
            PrincipalAxes::uniaxial(sellmeier(ORDINARY), sellmeier(EXTRAORDINARY)),
            rotation,
        );
        assert_same_response(&tabulated, &analytical, &SAMPLES);
    }
}

#[test]
fn outside_the_table_a_clamped_axis_takes_its_endpoint_value() {
    // `Clamp` holds the endpoint value, so the stack evaluates as if the axis
    // were that constant.
    let beyond = [1.5e4, 4.0e4];
    let endpoint = |k: f64| {
        let curve = sellmeier(EXTRAORDINARY);
        let value: Complex64 = curve.relative_permittivity(Scalar(k));
        Sellmeier::new(value.re, Vec::new()).unwrap()
    };
    for rotation in rotations() {
        let tabulated = stack(
            PrincipalAxes::uniaxial(sellmeier(ORDINARY), table(Extrapolation::Clamp)),
            rotation,
        );
        for (&k0, &sample) in beyond.iter().zip([SAMPLES[0], SAMPLES[4]].iter()) {
            let constant = stack(
                PrincipalAxes::uniaxial(sellmeier(ORDINARY), endpoint(sample)),
                rotation,
            );
            assert_same_response(&tabulated, &constant, &[k0]);
        }
    }
}

/// Wavenumbers below and above the table, `cm⁻¹`.
const BEYOND: [f64; 2] = [1.5e4, 4.0e4];

/// The `ConstitutiveError` a 2×2 backend failed with, from an evaluation
/// error's source chain. The backends forward it transparently, so it is read
/// from the backend error rather than found in the chain.
fn constitutive_error(error: &(dyn Error + 'static)) -> Option<ConstitutiveError> {
    let mut current = Some(error);
    while let Some(error) = current {
        match error.downcast_ref() {
            Some(Scatter2Error::Constitutive(constitutive)) => return Some(*constitutive),
            Some(_) => return None,
            None => {}
        }
        match error.downcast_ref() {
            Some(Transfer2Error::Constitutive(constitutive)) => return Some(*constitutive),
            Some(_) => return None,
            None => {}
        }
        current = error.source();
    }
    None
}

#[test]
fn a_rejecting_principal_axis_matches_a_clamped_one_inside_the_table() {
    for rotation in rotations() {
        let axes =
            |extrapolation| PrincipalAxes::uniaxial(sellmeier(ORDINARY), table(extrapolation));
        let rejecting = stack(axes(Extrapolation::Reject), rotation);
        let clamped = stack(axes(Extrapolation::Clamp), rotation);
        assert_same_response(&rejecting, &clamped, &SAMPLES);
    }
}

#[test]
fn outside_the_table_a_rejecting_axis_is_a_non_finite_layer() {
    let expected = ConstitutiveError::NonFinite {
        region: StackRegion::Layer(FiniteLayerIndex::new(0)),
    };
    for rotation in rotations() {
        let tabulated = stack(
            PrincipalAxes::uniaxial(sellmeier(ORDINARY), table(Extrapolation::Reject)),
            rotation,
        );
        for k0 in BEYOND {
            for polarisation in POLARISATIONS {
                for_each_backend!(evaluator, {
                    let error = evaluator
                        .evaluate(input(&[k0]), &tabulated, polarisation)
                        .unwrap_err();
                    assert_eq!(
                        constitutive_error(&error),
                        Some(expected),
                        "{k0} {polarisation:?}: {error}"
                    );
                });
            }
        }
    }
}

/// lamina-core validates the stack-frame tensors of oriented regions only: a
/// scalar region off the table is not refused by name. The scattering-matrix
/// backend returns a NaN response and the transfer-matrix backend refuses a
/// non-finite final matrix.
#[test]
fn outside_the_table_a_rejecting_scalar_layer_gives_no_finite_response() {
    let stack: MaterialStack<Complex64> = Stack::from_materials(
        Sellmeier::new(1.0, Vec::new()).unwrap(),
        sellmeier(CLADDING),
    )
    .material_layer(table(Extrapolation::Reject), Length::micrometres(0.8))
    .finalise();
    for k0 in BEYOND {
        for polarisation in POLARISATIONS {
            let response = RealAxisEvaluator::new(Scatter2::new())
                .evaluate(input(&[k0]), &stack, polarisation)
                .unwrap();
            for side in [IncidentSide::Left, IncidentSide::Right] {
                let amplitudes = response.amplitudes(side).unwrap();
                assert!(
                    amplitudes
                        .reflection()
                        .iter()
                        .chain(amplitudes.transmission())
                        .all(|value| value.re.is_nan()),
                    "{k0} {polarisation:?} {side:?}"
                );
            }

            let error = RealAxisEvaluator::new(Transfer2::new())
                .evaluate(input(&[k0]), &stack, polarisation)
                .unwrap_err();
            let mut source = Some(&error as &(dyn Error + 'static));
            let mut transfer = None;
            while let Some(error) = source {
                transfer = transfer.or(error.downcast_ref::<Transfer2Error>());
                source = error.source();
            }
            assert!(
                matches!(transfer, Some(Transfer2Error::NonFiniteFinalMatrix { .. })),
                "{k0} {polarisation:?}: {error}"
            );
        }
    }
}
