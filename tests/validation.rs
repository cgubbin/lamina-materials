use lamina_materials::{
    Cauchy, Drude, Extrapolation, LorentzOscillator, MaterialModelError, SellmeierTerm, Tabulated,
    TabulatedMaterialError,
};
use lamina_units::{VacuumAngularWavenumber, VacuumAngularWavenumberSquared};
use ndarray::array;
use num_complex::Complex64;

fn spectral(value: f64) -> VacuumAngularWavenumber<f64> {
    VacuumAngularWavenumber::inverse_centimetres(value)
}

#[test]
fn negative_drude_plasma_wavenumber_is_rejected() {
    assert!(matches!(
        Drude::new(1.0, spectral(-1.0), spectral(0.1),),
        Err(MaterialModelError::NegativeParameter {
            name: "plasma_wavenumber",
            ..
        })
    ));
}

#[test]
fn zero_lorentz_resonance_is_rejected() {
    assert!(matches!(
        LorentzOscillator::from_delta_epsilon(1.0, spectral(0.0), spectral(0.1),),
        Err(MaterialModelError::NonPositiveParameter {
            name: "resonance",
            ..
        })
    ));
}

#[test]
fn negative_lorentz_strength_is_rejected() {
    assert!(matches!(
        LorentzOscillator::from_strength(
            VacuumAngularWavenumberSquared::inverse_centimetres_squared(-1.0),
            spectral(2.0),
            spectral(0.1),
        ),
        Err(MaterialModelError::NegativeParameter {
            name: "strength",
            ..
        })
    ));
}

#[test]
fn oscillator_delta_epsilon_constructor_is_consistent() {
    let oscillator =
        LorentzOscillator::from_delta_epsilon(2.0, spectral(3.0), spectral(0.1)).unwrap();

    assert_eq!(oscillator.strength(), 18.0);
    assert_eq!(oscillator.delta_epsilon(), 2.0);
}

#[test]
fn sellmeier_zero_pole_is_rejected() {
    assert!(matches!(
        SellmeierTerm::new(1.0, spectral(0.0),),
        Err(MaterialModelError::NonPositiveParameter {
            name: "pole_wavenumber",
            ..
        })
    ));
}

#[test]
fn empty_cauchy_coefficients_are_rejected() {
    assert!(matches!(
        Cauchy::<f64>::from_spectral_basis(Vec::new(), spectral(1.0),),
        Err(MaterialModelError::EmptyCoefficients {
            name: "coefficients",
        })
    ));
}

#[test]
fn non_finite_analytical_parameters_are_rejected() {
    assert!(matches!(
        Drude::new(f64::NAN, spectral(1.0), spectral(0.1),),
        Err(MaterialModelError::NonFiniteParameter {
            name: "epsilon_infinity",
            ..
        })
    ));

    assert!(matches!(
        LorentzOscillator::from_delta_epsilon(1.0, spectral(f64::INFINITY), spectral(0.1),),
        Err(MaterialModelError::NonPositiveParameter {
            name: "resonance",
            ..
        }) | Err(MaterialModelError::NonFiniteParameter {
            name: "resonance",
            ..
        })
    ));

    assert!(matches!(
        SellmeierTerm::new(f64::NAN, spectral(1.0),),
        Err(MaterialModelError::NonFiniteParameter {
            name: "coefficient",
            ..
        })
    ));

    assert!(matches!(
        Cauchy::from_spectral_basis(vec![1.0, f64::NAN], spectral(1.0),),
        Err(MaterialModelError::NonFiniteParameter {
            name: "coefficient",
            ..
        })
    ));
}

#[test]
fn negative_damping_is_rejected() {
    assert!(matches!(
        Drude::new(1.0, spectral(2.0), spectral(-0.1),),
        Err(MaterialModelError::NegativeParameter { .. })
    ));

    assert!(matches!(
        LorentzOscillator::from_delta_epsilon(1.0, spectral(2.0), spectral(-0.1),),
        Err(MaterialModelError::NegativeParameter {
            name: "damping",
            ..
        })
    ));
}

#[test]
fn tabulated_coordinates_must_increase() {
    let result = Tabulated::new(
        array![1.0, 1.0],
        array![Complex64::new(2.0, 0.0), Complex64::new(3.0, 0.0),],
        Extrapolation::Clamp,
    );

    assert!(matches!(
        result,
        Err(TabulatedMaterialError::NonIncreasingCoordinates { .. })
    ));
}

#[test]
fn tabulated_lengths_must_match() {
    let result = Tabulated::new(
        array![1.0, 2.0],
        array![Complex64::new(2.0, 0.0)],
        Extrapolation::Clamp,
    );

    assert!(matches!(
        result,
        Err(TabulatedMaterialError::LengthMismatch { .. })
    ));
}

#[test]
fn tabulated_requires_two_samples() {
    let result = Tabulated::new(
        array![1.0],
        array![Complex64::new(2.0, 0.0)],
        Extrapolation::Clamp,
    );

    assert!(matches!(
        result,
        Err(TabulatedMaterialError::TooFewSamples { len: 1 })
    ));
}

#[test]
fn tabulated_rejects_non_finite_coordinate() {
    let result = Tabulated::new(
        array![1.0, f64::NAN],
        array![Complex64::new(2.0, 0.0), Complex64::new(3.0, 0.0),],
        Extrapolation::Clamp,
    );

    assert!(matches!(
        result,
        Err(TabulatedMaterialError::NonFiniteCoordinate { index: 1, .. })
    ));
}

#[test]
fn tabulated_rejects_non_finite_value() {
    let result = Tabulated::new(
        array![1.0, 2.0],
        array![Complex64::new(2.0, 0.0), Complex64::new(f64::NAN, 0.0),],
        Extrapolation::Clamp,
    );

    assert!(matches!(
        result,
        Err(TabulatedMaterialError::NonFiniteValue { index: 1 })
    ));
}
