use lamina_core::{DerivativeOrder, DifferentiableMaterial, DifferentiableMeromorphicMaterial};
use lamina_materials::{
    Cauchy, Drude, DrudeLorentz, Lorentz, LorentzOscillator, Sellmeier, SellmeierTerm,
};
use ndarray::{arr0, arr1};
use num_complex::Complex64;

const TOLERANCE: f64 = 1e-12;

fn assert_complex_close(actual: Complex64, expected: Complex64) {
    let error = (actual - expected).norm();
    let scale = expected.norm().max(1.0);

    assert!(
        error <= TOLERANCE * scale,
        "actual = {actual:?}, expected = {expected:?}, relative error = {}",
        error / scale
    );
}

fn assert_real_and_complex_derivatives_agree<M>(material: &M, k0: f64)
where
    M: DifferentiableMaterial<Real = f64> + DifferentiableMeromorphicMaterial<Real = f64>,
{
    for order in [
        DerivativeOrder::First,
        DerivativeOrder::Second,
        DerivativeOrder::Third,
    ] {
        let real = material.relative_permittivity_derivative::<_, Complex64>(arr0(k0), order);

        let complex = material.relative_permittivity_complex_derivative::<_, Complex64>(
            arr0(Complex64::new(k0, 0.0)),
            order,
        );

        assert_complex_close(real[()], complex[()]);
    }
}

fn assert_scalar_and_sampled_derivatives_agree<M>(material: &M, k0: f64)
where
    M: DifferentiableMaterial<Real = f64>,
{
    for order in [
        DerivativeOrder::First,
        DerivativeOrder::Second,
        DerivativeOrder::Third,
    ] {
        let scalar = material.relative_permittivity_derivative::<_, Complex64>(arr0(k0), order);

        let sampled = material.relative_permittivity_derivative::<_, Complex64>(arr1(&[k0]), order);

        assert_complex_close(sampled[0], scalar[()]);
    }
}

fn assert_scalar_and_sampled_complex_derivatives_agree<M>(material: &M, k0: Complex64)
where
    M: DifferentiableMeromorphicMaterial<Real = f64>,
{
    for order in [
        DerivativeOrder::First,
        DerivativeOrder::Second,
        DerivativeOrder::Third,
    ] {
        let scalar =
            material.relative_permittivity_complex_derivative::<_, Complex64>(arr0(k0), order);

        let sampled =
            material.relative_permittivity_complex_derivative::<_, Complex64>(arr1(&[k0]), order);

        assert_complex_close(sampled[0], scalar[()]);
    }
}

#[test]
fn drude_public_derivative_contract_is_consistent() {
    let material = Drude::new(1.2, 2.4, 0.3).unwrap();

    assert_real_and_complex_derivatives_agree(&material, 1.1);
    assert_scalar_and_sampled_derivatives_agree(&material, 1.1);
    assert_scalar_and_sampled_complex_derivatives_agree(&material, Complex64::new(1.1, 0.2));
}

#[test]
fn lorentz_public_derivative_contract_is_consistent() {
    let oscillator = LorentzOscillator::new(1.3, 2.7, 0.2).unwrap();

    let material = Lorentz::new(1.4, vec![oscillator]).unwrap();

    assert_real_and_complex_derivatives_agree(&material, 1.0);
    assert_scalar_and_sampled_derivatives_agree(&material, 1.0);
    assert_scalar_and_sampled_complex_derivatives_agree(&material, Complex64::new(1.0, 0.15));
}

#[test]
fn drude_lorentz_public_derivative_contract_is_consistent() {
    let oscillator = LorentzOscillator::new(0.9, 3.0, 0.25).unwrap();

    let material = DrudeLorentz::new(1.3, 1.8, 0.35, vec![oscillator]).unwrap();

    assert_real_and_complex_derivatives_agree(&material, 1.2);
    assert_scalar_and_sampled_derivatives_agree(&material, 1.2);
    assert_scalar_and_sampled_complex_derivatives_agree(&material, Complex64::new(1.2, 0.1));
}

#[test]
fn sellmeier_public_derivative_contract_is_consistent() {
    let material = Sellmeier::new(
        1.0,
        vec![
            SellmeierTerm::new(0.7, 2.8).unwrap(),
            SellmeierTerm::new(0.3, 4.1).unwrap(),
        ],
    )
    .unwrap();

    assert_real_and_complex_derivatives_agree(&material, 0.9);
    assert_scalar_and_sampled_derivatives_agree(&material, 0.9);
    assert_scalar_and_sampled_complex_derivatives_agree(&material, Complex64::new(0.9, 0.12));
}

#[test]
fn cauchy_public_derivative_contract_is_consistent() {
    let material = Cauchy::new(vec![1.4, 0.3, -0.07, 0.012]).unwrap();

    assert_real_and_complex_derivatives_agree(&material, 0.8);
    assert_scalar_and_sampled_derivatives_agree(&material, 0.8);
    assert_scalar_and_sampled_complex_derivatives_agree(&material, Complex64::new(0.8, 0.1));
}

#[test]
fn drude_public_derivative_orders_match_closed_form() {
    let epsilon_infinity = 1.5;
    let plasma_wavenumber: f64 = 2.0;
    let k0: f64 = 1.25;

    let material = Drude::new(epsilon_infinity, plasma_wavenumber, 0.0).unwrap();

    let strength = plasma_wavenumber.powi(2);

    let first =
        material.relative_permittivity_derivative::<_, Complex64>(arr0(k0), DerivativeOrder::First);

    let second = material
        .relative_permittivity_derivative::<_, Complex64>(arr0(k0), DerivativeOrder::Second);

    let third =
        material.relative_permittivity_derivative::<_, Complex64>(arr0(k0), DerivativeOrder::Third);

    assert_complex_close(first[()], Complex64::new(2.0 * strength / k0.powi(3), 0.0));

    assert_complex_close(
        second[()],
        Complex64::new(-6.0 * strength / k0.powi(4), 0.0),
    );

    assert_complex_close(third[()], Complex64::new(24.0 * strength / k0.powi(5), 0.0));
}
