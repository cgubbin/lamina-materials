use lamina_materials::{
    Cauchy, Drude, Extrapolation, LorentzOscillator, MaterialModelError, SellmeierTerm, Tabulated,
    TabulatedMaterialError,
};
use ndarray::array;
use num_complex::Complex64;

#[test]
fn negative_drude_plasma_wavenumber_is_rejected() {
    assert!(matches!(
        Drude::new(1.0, -1.0, 0.1),
        Err(MaterialModelError::NegativeParameter {
            name: "plasma_wavenumber",
            ..
        })
    ));
}

#[test]
fn zero_lorentz_resonance_is_rejected() {
    assert!(matches!(
        LorentzOscillator::new(1.0, 0.0, 0.1),
        Err(MaterialModelError::NonPositiveParameter {
            name: "resonance",
            ..
        })
    ));
}

#[test]
fn oscillator_delta_epsilon_constructor_is_consistent() {
    let oscillator = LorentzOscillator::from_delta_epsilon(2.0, 3.0, 0.1).unwrap();

    assert_eq!(oscillator.strength(), 18.0);
    assert_eq!(oscillator.delta_epsilon(), 2.0);
}

#[test]
fn sellmeier_zero_pole_is_rejected() {
    assert!(SellmeierTerm::new(1.0, 0.0).is_err());
}

#[test]
fn empty_cauchy_coefficients_are_rejected() {
    assert!(matches!(
        Cauchy::<f64>::new(Vec::new()),
        Err(MaterialModelError::EmptyCoefficients {
            name: "coefficients"
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
