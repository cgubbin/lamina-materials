use approx::assert_relative_eq;
use lamina_core::{Material, Scalar};
use lamina_materials::{Extrapolation, Tabulated};
use ndarray::array;
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
