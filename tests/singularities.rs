//! Declared poles and zeros of the rational models.
use lamina_core::{
    AnalyticalMaterialHandle, AnalyticalTensorMaterialHandle, ConstitutiveQuantity,
    MeromorphicMaterial, Scalar, SingularityKind, SpectralSingularity,
};
use lamina_materials::{
    DrudeLorentz, Lorentz, LorentzOscillator, PrincipalAxes, Sellmeier, SellmeierTerm,
};
use lamina_units::{VacuumAngularWavenumber, VacuumWavenumber};
use num_complex::Complex64 as C;
use std::f64::consts::TAU;

/// ε(ν) = ε∞ (1 + (ν_L² − ν_T²)/(ν_T² − ν² − iγν)) at k0 = 2πν.
fn phonon(eps_inf: f64, to: f64, lo: f64, gamma: f64) -> Lorentz<f64> {
    let oscillator = LorentzOscillator::from_delta_epsilon(
        eps_inf * (lo * lo / (to * to) - 1.0),
        VacuumWavenumber::inverse_centimetres(to),
        VacuumWavenumber::inverse_centimetres(gamma),
    )
    .unwrap();
    Lorentz::new(eps_inf, vec![oscillator]).unwrap()
}

fn of_kind(points: &[SpectralSingularity<C>], kind: SingularityKind) -> Vec<C> {
    let mut found: Vec<C> = points
        .iter()
        .filter(|p| p.kind == kind)
        .map(|p| p.k0)
        .collect();
    found.sort_by(|a, b| a.re.total_cmp(&b.re));
    found
}

fn close(found: &[C], expected: &[C]) {
    assert_eq!(found.len(), expected.len(), "{found:?}");
    for (f, e) in found.iter().zip(expected) {
        assert!((f - e).norm() <= 1e-12 * e.norm(), "{f} vs {e}");
    }
}

#[test]
fn a_lorentz_film_declares_its_transverse_and_longitudinal_points() {
    let (to, lo, gamma) = (780.0, 830.0, 4.0);
    let points = phonon(2.95, to, lo, gamma)
        .spectral_singularities::<C>()
        .unwrap();
    assert!(
        points
            .iter()
            .all(|p| p.quantity == ConstitutiveQuantity::Permittivity)
    );
    assert!(points.iter().all(|p| p.axis.is_none()));
    let at = |nu: f64| {
        let re = (nu * nu - gamma * gamma / 4.0).sqrt();
        [
            C::new(-re, -gamma / 2.0) * TAU,
            C::new(re, -gamma / 2.0) * TAU,
        ]
    };
    close(&of_kind(&points, SingularityKind::Pole), &at(to));
    close(&of_kind(&points, SingularityKind::Zero), &at(lo));
}

#[test]
fn drude_lorentz_zeros_are_zeros_of_the_permittivity() {
    let oscillators = [(1.2, 900.0, 15.0), (0.4, 1500.0, 40.0)]
        .map(|(delta, resonance, damping)| {
            LorentzOscillator::from_delta_epsilon(
                delta,
                VacuumAngularWavenumber::inverse_centimetres(resonance),
                VacuumAngularWavenumber::inverse_centimetres(damping),
            )
            .unwrap()
        })
        .to_vec();
    let model = DrudeLorentz::from_plasma_wavenumber(
        2.5,
        VacuumAngularWavenumber::inverse_centimetres(3000.0),
        VacuumAngularWavenumber::inverse_centimetres(60.0),
        oscillators,
    )
    .unwrap();
    let points = model.spectral_singularities::<C>().unwrap();
    let poles = of_kind(&points, SingularityKind::Pole);
    let zeros = of_kind(&points, SingularityKind::Zero);
    assert_eq!((poles.len(), zeros.len()), (6, 6));
    assert!(poles.iter().any(|p| p.norm() < 1e-12));
    assert!(
        poles
            .iter()
            .any(|p| (p - C::new(0.0, -60.0)).norm() < 1e-12)
    );
    for zero in zeros {
        let epsilon: C = model.relative_permittivity_complex(Scalar(zero));
        assert!(epsilon.norm() < 1e-11, "ε({zero}) = {epsilon}");
    }
    // The erased handle forwards the declaration.
    assert_eq!(
        AnalyticalMaterialHandle::<C>::new(model).spectral_singularities(),
        Some(points)
    );
}

#[test]
fn a_uniaxial_medium_declares_each_axis_on_its_entries() {
    let (ordinary, extraordinary) = (
        phonon(4.87, 1370.0, 1610.0, 5.0),
        phonon(2.95, 780.0, 830.0, 4.0),
    );
    let expected = |model: &Lorentz<f64>, axes: &[usize]| {
        let points = model.spectral_singularities::<C>().unwrap();
        axes.iter()
            .flat_map(|&axis| points.iter().map(move |p| p.on_axis(axis)))
            .collect::<Vec<_>>()
    };
    let mut all = expected(&ordinary, &[0, 1]);
    all.extend(expected(&extraordinary, &[2]));
    let handle = AnalyticalTensorMaterialHandle::<C>::new(PrincipalAxes::uniaxial(
        ordinary.clone(),
        extraordinary,
    ));
    assert_eq!(handle.spectral_singularities(), Some(all));
    // A Sellmeier axis declares nothing, so neither does the medium.
    let sellmeier = Sellmeier::new(
        1.0,
        vec![SellmeierTerm::new(1.0, VacuumWavenumber::inverse_centimetres(10000.0)).unwrap()],
    )
    .unwrap();
    let mixed = AnalyticalTensorMaterialHandle::<C>::new(PrincipalAxes::biaxial(
        ordinary.clone(),
        ordinary,
        sellmeier,
    ));
    assert_eq!(mixed.spectral_singularities(), None);
}
