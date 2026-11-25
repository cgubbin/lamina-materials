//! A `PrincipalAxes` model with a tabulated axis placed in a lamina-core
//! `TensorStack` and evaluated on the real axis through both 2×2 backends.

use std::f64::consts::FRAC_PI_2;

use lamina_core::{
    ConstitutiveClass, CoordinateInput, Coordinates, InPlaneCoordinate, IncidentSide, Material,
    MaterialRegion, Polarisation, RealAxisEvaluator, Rotation, Scalar, Scatter2,
    SpectralCoordinate, Stack, TensorMaterial, TensorMaterialHandle, TensorStack, Transfer2,
    material::MaterialHandle, units::InverseLengthUnit,
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
fn table() -> Tabulated<f64> {
    let curve = sellmeier(EXTRAORDINARY);
    let values = SAMPLES
        .iter()
        .map(|&k| curve.relative_permittivity(Scalar(k)))
        .collect::<Array1<Complex64>>();
    Tabulated::new(arr1(&SAMPLES), values, Extrapolation::Clamp).unwrap()
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
            PrincipalAxes::uniaxial(sellmeier(ORDINARY), table()),
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
    // `Tabulated` has no refusing extrapolation: `Clamp` holds the endpoint
    // value, so the stack evaluates as if the axis were that constant.
    let beyond = [1.5e4, 4.0e4];
    let endpoint = |k: f64| {
        let curve = sellmeier(EXTRAORDINARY);
        let value: Complex64 = curve.relative_permittivity(Scalar(k));
        Sellmeier::new(value.re, Vec::new()).unwrap()
    };
    for rotation in rotations() {
        let tabulated = stack(
            PrincipalAxes::uniaxial(sellmeier(ORDINARY), table()),
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
