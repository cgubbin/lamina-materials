use ndarray::array;
use num_complex::Complex64;
use tmm_core::material::{
    AnalyticalMaterialHandle, EvaluateDifferentiableMaterial, EvaluateMaterial, MaterialHandle,
};
use tmm_materials::{Drude, Extrapolation, Tabulated};

#[test]
fn erased_real_axis_evaluation_matches_direct_evaluation() {
    let model = Drude::new(1.0, 2.0, 0.1).unwrap();

    let input = array![0.8, 1.0, 1.2,];

    let direct: ndarray::Array1<Complex64> = model.evaluate_relative_permittivity(input.clone());

    let handle = MaterialHandle::<f64, Complex64>::new(model);

    let erased = handle.evaluate_relative_permittivity(input);

    assert_eq!(direct, erased,);
}

#[test]
fn heterogeneous_real_axis_handles_share_one_type() {
    let drude = MaterialHandle::<f64, Complex64>::new(Drude::new(1.0, 2.0, 0.1).unwrap());

    let tabulated = MaterialHandle::<f64, Complex64>::new(
        Tabulated::new(
            array![1.0, 2.0],
            array![Complex64::new(2.0, 0.1), Complex64::new(3.0, 0.2),],
            Extrapolation::Clamp,
        )
        .unwrap(),
    );

    let materials = vec![drude, tabulated];

    assert_eq!(materials.len(), 2,);
}

#[test]
fn analytical_handle_exposes_derivatives() {
    let handle =
        AnalyticalMaterialHandle::<f64, Complex64>::new(Drude::new(1.0, 2.0, 0.1).unwrap());

    let result = handle.evaluate_relative_permittivity_derivative(
        array![1.0, 1.5],
        tmm_core::DerivativeOrder::First,
        tmm_core::SpectralVariable::VacuumWavenumber,
    );

    assert_eq!(result.shape(), &[2],);
}
