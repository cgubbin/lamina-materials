use lamina_core::material::{
    AnalyticalMaterialHandle, EvaluateDifferentiableMaterial, EvaluateMaterial, MaterialHandle,
};
use lamina_materials::{Drude, Extrapolation, Tabulated};
use lamina_units::VacuumAngularWavenumber;
use ndarray::array;
use num_complex::Complex64;

fn spectral(value: f64) -> VacuumAngularWavenumber<f64> {
    VacuumAngularWavenumber::inverse_centimetres(value)
}

fn drude() -> Drude<f64> {
    Drude::new(1.0, spectral(2.0), spectral(0.1)).unwrap()
}

#[test]
fn erased_real_axis_evaluation_matches_direct_evaluation() {
    let model = drude();

    let input = array![0.8, 1.0, 1.2];

    let direct: ndarray::Array1<Complex64> = model.evaluate_relative_permittivity(input.clone());

    let handle = MaterialHandle::<Complex64>::new(model);

    let erased = handle.evaluate_relative_permittivity(input);

    assert_eq!(direct, erased);
}

#[test]
fn heterogeneous_real_axis_handles_share_one_type() {
    let drude = MaterialHandle::<Complex64>::new(drude());

    let tabulated = MaterialHandle::<Complex64>::new(
        Tabulated::new(
            array![1.0, 2.0],
            array![Complex64::new(2.0, 0.1), Complex64::new(3.0, 0.2),],
            Extrapolation::Clamp,
        )
        .unwrap(),
    );

    let materials = vec![drude, tabulated];

    assert_eq!(materials.len(), 2);
}

#[test]
fn analytical_handle_exposes_derivatives() {
    let handle = AnalyticalMaterialHandle::<Complex64>::new(drude());

    let result = handle.evaluate_relative_permittivity_derivative(
        array![1.0, 1.5],
        lamina_core::DerivativeOrder::First,
    );

    assert_eq!(result.shape(), &[2]);
}
