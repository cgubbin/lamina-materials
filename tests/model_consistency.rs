use approx::assert_relative_eq;
use num_complex::Complex64;
use tmm_core::{
    DerivativeOrder, DifferentiableMaterial, DifferentiableMeromorphicMaterial, Material,
    MeromorphicMaterial, Scalar, SpectralVariable,
};
use tmm_materials::{
    Cauchy, Drude, DrudeLorentz, Lorentz, LorentzOscillator, Sellmeier, SellmeierTerm,
};

fn assert_complex_close(left: Complex64, right: Complex64, relative: f64, absolute: f64) {
    assert_relative_eq!(
        left.re,
        right.re,
        max_relative = relative,
        epsilon = absolute,
    );
    assert_relative_eq!(
        left.im,
        right.im,
        max_relative = relative,
        epsilon = absolute,
    );
}

fn value<M>(material: &M, k0: f64) -> Complex64
where
    M: Material<Real = f64>,
{
    material.relative_permittivity::<_, Complex64>(Scalar(k0))
}

fn derivative<M>(
    material: &M,
    k0: f64,
    order: DerivativeOrder,
    variable: SpectralVariable,
) -> Complex64
where
    M: Material<Real = f64> + DifferentiableMaterial,
{
    material.relative_permittivity_derivative::<_, Complex64>(Scalar(k0), order, variable)
}

fn check_real_complex_consistency<M>(material: &M, k0: f64)
where
    M: Material<Real = f64>
        + DifferentiableMaterial
        + MeromorphicMaterial
        + DifferentiableMeromorphicMaterial,
{
    let real = value(material, k0);
    let complex = material.relative_permittivity_complex(Scalar(Complex64::new(k0, 0.0)));

    assert_complex_close(real, complex, 1.0e-12, 1.0e-12);

    for order in [DerivativeOrder::First, DerivativeOrder::Second] {
        for variable in [
            SpectralVariable::VacuumWavenumber,
            SpectralVariable::VacuumWavenumberSquared,
        ] {
            let real = derivative(material, k0, order, variable);
            let complex = material.relative_permittivity_complex_derivative(
                Scalar(Complex64::new(k0, 0.0)),
                order,
                variable,
            );

            assert_complex_close(real, complex, 1.0e-11, 1.0e-11);
        }
    }
}

fn check_derivatives<M>(material: &M, k0: f64)
where
    M: Material<Real = f64> + DifferentiableMaterial,
{
    let h = 2.0e-5 * k0.abs().max(1.0);
    let centre = value(material, k0);
    let plus = value(material, k0 + h);
    let minus = value(material, k0 - h);

    let first_fd = (plus - minus) / (2.0 * h);
    let second_fd = (plus - 2.0 * centre + minus) / (h * h);

    let first = derivative(
        material,
        k0,
        DerivativeOrder::First,
        SpectralVariable::VacuumWavenumber,
    );
    let second = derivative(
        material,
        k0,
        DerivativeOrder::Second,
        SpectralVariable::VacuumWavenumber,
    );

    assert_complex_close(first, first_fd, 2.0e-6, 2.0e-7);
    assert_complex_close(second, second_fd, 2.0e-4, 2.0e-5);

    let first_squared = derivative(
        material,
        k0,
        DerivativeOrder::First,
        SpectralVariable::VacuumWavenumberSquared,
    );
    let second_squared = derivative(
        material,
        k0,
        DerivativeOrder::Second,
        SpectralVariable::VacuumWavenumberSquared,
    );

    assert_complex_close(first_squared, first / (2.0 * k0), 1.0e-12, 1.0e-12);

    assert_complex_close(
        second_squared,
        second / (4.0 * k0 * k0) - first / (4.0 * k0 * k0 * k0),
        1.0e-12,
        1.0e-12,
    );
}

fn check_model<M>(material: &M)
where
    M: Material<Real = f64>
        + DifferentiableMaterial
        + MeromorphicMaterial
        + DifferentiableMeromorphicMaterial,
{
    for k0 in [0.7, 1.1, 1.6] {
        check_real_complex_consistency(material, k0);
        check_derivatives(material, k0);
    }
}

#[test]
fn analytical_models_are_consistent() {
    let drude = Drude::new(1.2, 2.0, 0.15).unwrap();

    let lorentz = Lorentz::new(
        1.5,
        vec![LorentzOscillator::from_delta_epsilon(0.8, 2.5, 0.2).unwrap()],
    )
    .unwrap();

    let combined = DrudeLorentz::from_plasma_wavenumber(
        1.1,
        1.7,
        0.12,
        vec![LorentzOscillator::new(3.0, 2.6, 0.18).unwrap()],
    )
    .unwrap();

    let sellmeier = Sellmeier::with_unit_constant(vec![
        SellmeierTerm::new(0.7, 3.0).unwrap(),
        SellmeierTerm::new(0.2, 5.0).unwrap(),
    ])
    .unwrap();

    let cauchy = Cauchy::new(vec![1.4, 0.02, 0.001]).unwrap();

    check_model(&drude);
    check_model(&lorentz);
    check_model(&combined);
    check_model(&sellmeier);
    check_model(&cauchy);
}
