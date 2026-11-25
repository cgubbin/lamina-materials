use approx::assert_relative_eq;
use lamina_core::{Material, Scalar};
use lamina_materials::{Extrapolation, Tabulated};
use ndarray::{array, s};
use num_complex::Complex64;

fn table(extrapolation: Extrapolation) -> Tabulated<f64> {
    Tabulated::new(
        array![1.0, 2.0, 4.0],
        array![
            Complex64::new(2.0, 0.1),
            Complex64::new(4.0, 0.3),
            Complex64::new(8.0, 0.7),
        ],
        extrapolation,
    )
    .unwrap()
}

#[test]
fn interpolation_is_linear() {
    let epsilon = table(Extrapolation::Clamp).relative_permittivity::<_, Complex64>(Scalar(3.0));

    assert_relative_eq!(epsilon.re, 6.0);
    assert_relative_eq!(epsilon.im, 0.5);
}

#[test]
fn clamp_uses_endpoint_values() {
    let epsilon = table(Extrapolation::Clamp).relative_permittivity::<_, Complex64>(Scalar(0.0));

    assert_relative_eq!(epsilon.re, 2.0);
    assert_relative_eq!(epsilon.im, 0.1);
}

#[test]
fn linear_extrapolation_continues_endpoint_segment() {
    let epsilon = table(Extrapolation::Linear).relative_permittivity::<_, Complex64>(Scalar(0.0));

    assert_relative_eq!(epsilon.re, 0.0);
    assert_relative_eq!(epsilon.im, -0.1);
}

fn permittivity(table: &Tabulated<f64>, k0: f64) -> Complex64 {
    table.relative_permittivity::<_, Complex64>(Scalar(k0))
}

#[test]
fn reject_matches_the_other_policies_inside_the_table() {
    let reject = table(Extrapolation::Reject);
    // Endpoints are inside the table.
    for k0 in [1.0, 1.5, 2.0, 3.0, 4.0] {
        let epsilon = permittivity(&reject, k0);
        assert_eq!(epsilon, permittivity(&table(Extrapolation::Clamp), k0));
        assert_eq!(epsilon, permittivity(&table(Extrapolation::Linear), k0));
    }
}

#[test]
fn reject_gives_nan_outside_the_table() {
    let reject = table(Extrapolation::Reject);
    for k0 in [0.0, 1.0 - 1e-12, 4.0 + 1e-12, 5.0, f64::INFINITY] {
        let epsilon = permittivity(&reject, k0);
        assert!(epsilon.re.is_nan() && epsilon.im.is_nan(), "k0 = {k0}");
    }
}

#[test]
fn a_nan_wavenumber_gives_nan_under_every_policy() {
    for extrapolation in [
        Extrapolation::Clamp,
        Extrapolation::Linear,
        Extrapolation::Reject,
    ] {
        let epsilon = permittivity(&table(extrapolation), f64::NAN);
        assert!(
            epsilon.re.is_nan() && epsilon.im.is_nan(),
            "{extrapolation:?}"
        );
    }
}

#[test]
fn a_strided_table_interpolates() {
    let coordinates = array![1.0, 0.0, 2.0, 0.0, 4.0].slice_move(s![..;2]);
    assert!(coordinates.as_slice().is_none());
    let strided = Tabulated::new(
        coordinates,
        table(Extrapolation::Clamp).values().clone(),
        Extrapolation::Clamp,
    )
    .unwrap();

    assert_eq!(
        permittivity(&strided, 3.0),
        permittivity(&table(Extrapolation::Clamp), 3.0)
    );
}
