use lamina_core::{
    AnalyticalTensorMaterialHandle, DerivativeOrder, DifferentiableMaterial,
    DifferentiableMeromorphicTensorMaterial, DifferentiableTensorMaterial, Material,
    MeromorphicMaterial, MeromorphicTensorMaterial, Scalar, TensorMaterial, material::Tensor3,
};
use lamina_materials::{
    ConstantTensor, ConstantTensorError, CrystalAxis, Extrapolation, PrincipalAxes, PrincipalModel,
    Sellmeier, SellmeierAxes, SellmeierTerm, Tabulated, TensorParameter, TensorRate,
};
use lamina_units::VacuumAngularWavenumber;
use ndarray::array;
use num_complex::Complex64;

/// Whether `$type` implements `$trait`: the inherent constant shadows the
/// trait's fallback only when the bound holds, so one of the two is unused.
macro_rules! implements {
    ($type:ty: $trait:path) => {{
        struct Probe<T>(std::marker::PhantomData<T>);
        #[allow(dead_code)]
        trait Fallback {
            const IMPLEMENTS: bool = false;
        }
        impl<T> Fallback for Probe<T> {}
        impl<T: $trait> Probe<T> {
            #[allow(dead_code)]
            const IMPLEMENTS: bool = true;
        }
        Probe::<$type>::IMPLEMENTS
    }};
}

const K0: f64 = 4.0e4;

/// Sellmeier constant and `(B, q)` terms, `q` in `cm⁻¹`.
type Parameters = (f64, [(f64, f64); 2]);

const ORDINARY: Parameters = (1.0, [(3.9, 1.3e5), (2.6, 2.6e3)]);
const EXTRAORDINARY: Parameters = (1.2, [(3.6, 1.2e5), (2.4, 2.7e3)]);
const THIRD: Parameters = (1.1, [(3.0, 1.1e5), (1.9, 3.0e3)]);

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

fn biaxial(axes: [Parameters; 3]) -> SellmeierAxes<f64> {
    PrincipalAxes::biaxial(sellmeier(axes[0]), sellmeier(axes[1]), sellmeier(axes[2]))
}

fn diagonal([a, b, c]: [Complex64; 3]) -> Tensor3<Complex64> {
    let zero = Complex64::from(0.0);
    [[a, zero, zero], [zero, b, zero], [zero, zero, c]]
}

fn unit(row: usize, column: usize) -> Tensor3<Complex64> {
    let mut tensor = [[Complex64::from(0.0); 3]; 3];
    tensor[row][column] = Complex64::from(1.0);
    tensor
}

#[test]
fn uniaxial_sellmeier_is_diagonal_in_scalar_values_and_derivatives() {
    let (ordinary, extraordinary) = (sellmeier(ORDINARY), sellmeier(EXTRAORDINARY));
    let model = PrincipalAxes::uniaxial(ordinary.clone(), extraordinary.clone());
    const { assert!(<SellmeierAxes<f64> as TensorMaterial>::PRINCIPAL_AXES) };

    let epsilon: Tensor3<Complex64> = model.relative_permittivity_tensor(Scalar(K0));
    let o: Complex64 = ordinary.relative_permittivity(Scalar(K0));
    let e: Complex64 = extraordinary.relative_permittivity(Scalar(K0));
    assert_eq!(epsilon, diagonal([o, o, e]));
    let mu: Tensor3<Complex64> = model.relative_permeability_tensor(Scalar(K0));
    assert_eq!(mu, diagonal([Complex64::from(1.0); 3]));

    for order in [DerivativeOrder::First, DerivativeOrder::Second] {
        let derivative: Tensor3<Complex64> =
            model.relative_permittivity_tensor_derivative(Scalar(K0), order);
        let o: Complex64 = ordinary.relative_permittivity_derivative(Scalar(K0), order);
        let e: Complex64 = extraordinary.relative_permittivity_derivative(Scalar(K0), order);
        assert_eq!(derivative, diagonal([o, o, e]));
    }

    let k = Complex64::new(K0, 300.0);
    let epsilon: Tensor3<Complex64> = model.relative_permittivity_tensor_complex(Scalar(k));
    let o: Complex64 = ordinary.relative_permittivity_complex(Scalar(k));
    let e: Complex64 = extraordinary.relative_permittivity_complex(Scalar(k));
    assert_eq!(epsilon, diagonal([o, o, e]));

    let sampled = model.relative_permittivity_tensor::<_, Complex64>(array![1.0e4, K0]);
    assert_eq!(
        sampled[1],
        model.relative_permittivity_tensor::<_, Complex64>(Scalar(K0))
    );
}

#[test]
fn analytical_axes_give_an_analytical_tensor_handle() {
    let model = PrincipalAxes::uniaxial(sellmeier(ORDINARY), sellmeier(EXTRAORDINARY));
    let handle = AnalyticalTensorMaterialHandle::<Complex64>::new(model);
    assert!(format!("{handle:?}").contains("principal_axes: true"));
}

#[test]
fn a_tabulated_axis_restricts_the_model_to_the_real_axis() {
    let table = Tabulated::new(
        array![1.0e4, 5.0e4],
        array![Complex64::new(4.0, 0.1), Complex64::new(5.0, 0.3)],
        Extrapolation::Clamp,
    )
    .unwrap();
    let ordinary = sellmeier(ORDINARY);
    let model = PrincipalAxes::uniaxial(ordinary.clone(), table.clone());
    type Mixed = PrincipalAxes<Sellmeier<f64>, Sellmeier<f64>, Tabulated<f64>>;

    let epsilon: Tensor3<Complex64> = model.relative_permittivity_tensor(Scalar(K0));
    let o: Complex64 = ordinary.relative_permittivity(Scalar(K0));
    let e: Complex64 = table.relative_permittivity(Scalar(K0));
    assert_eq!(epsilon, diagonal([o, o, e]));

    assert!(implements!(Mixed: TensorMaterial));
    assert!(!implements!(Mixed: DifferentiableTensorMaterial));
    assert!(!implements!(Mixed: MeromorphicTensorMaterial));
    assert!(implements!(SellmeierAxes<f64>: DifferentiableMeromorphicTensorMaterial));
}

/// Check every rate of the model `build` makes from per-model parameters
/// against central differences of its full permittivity tensor; `slot` names
/// the entry of `models` a rate's scalar model comes from.
fn checked_rates<const N: usize>(
    models: [Parameters; N],
    build: impl Fn([Parameters; N]) -> SellmeierAxes<f64>,
    slot: impl Fn(PrincipalModel) -> usize,
) -> Vec<TensorRate<f64>> {
    let rates = build(models).permittivity_parameter_rates(K0);
    assert_eq!(rates.len(), N * 5);

    for rate in &rates {
        let TensorParameter::Principal { model, name, index } = rate.parameter else {
            panic!("unexpected parameter {:?}", rate.parameter);
        };
        let perturbed = |step: f64| {
            let mut models = models;
            let (constant, terms) = &mut models[slot(model)];
            let value = match (name, index) {
                ("constant", None) => constant,
                ("coefficient", Some(j)) => &mut terms[j].0,
                ("pole_wavenumber", Some(j)) => &mut terms[j].1,
                other => panic!("unexpected parameter {other:?}"),
            };
            let h = 1e-6 * value.abs().max(1.0);
            *value += step * h;
            (
                build(models).relative_permittivity_tensor::<_, Complex64>(Scalar(K0)),
                h,
            )
        };
        let ((plus, h), (minus, _)) = (perturbed(1.0), perturbed(-1.0));
        let difference = |i: usize, j: usize| (plus[i][j] - minus[i][j]) / (2.0 * h);
        let scale = (0..3).map(|i| difference(i, i).norm()).fold(0.0, f64::max);
        for i in 0..3 {
            for j in 0..3 {
                let error = (rate.value[i][j] - difference(i, j)).norm();
                assert!(
                    error <= 1e-7 * scale,
                    "{:?} entry ({i}, {j}): rate {}, difference {}",
                    rate.parameter,
                    rate.value[i][j],
                    difference(i, j),
                );
            }
        }
    }
    rates
}

fn coefficient_rate(rates: &[TensorRate<f64>], model: PrincipalModel) -> Tensor3<Complex64> {
    let parameter = TensorParameter::Principal {
        model,
        name: "coefficient",
        index: Some(0),
    };
    rates
        .iter()
        .find(|rate| rate.parameter == parameter)
        .unwrap()
        .value
}

#[test]
fn uniaxial_rates_are_one_per_parameter_and_match_central_differences() {
    let rates = checked_rates(
        [ORDINARY, EXTRAORDINARY],
        |[ordinary, extraordinary]| {
            PrincipalAxes::uniaxial(sellmeier(ordinary), sellmeier(extraordinary))
        },
        |model| match model {
            PrincipalModel::Ordinary => 0,
            PrincipalModel::Extraordinary => 1,
            PrincipalModel::Axis(_) => panic!("biaxial parameter in a uniaxial model"),
        },
    );

    let ordinary = coefficient_rate(&rates, PrincipalModel::Ordinary);
    let o = ordinary[0][0];
    assert!(o.norm() > 0.0);
    assert_eq!(ordinary, diagonal([o, o, 0.0.into()]));

    let extraordinary = coefficient_rate(&rates, PrincipalModel::Extraordinary);
    let e = extraordinary[2][2];
    assert!(e.norm() > 0.0);
    assert_eq!(extraordinary, diagonal([0.0.into(), 0.0.into(), e]));
}

#[test]
fn biaxial_rates_are_per_axis_and_match_central_differences() {
    let rates = checked_rates([ORDINARY, EXTRAORDINARY, THIRD], biaxial, |model| {
        let PrincipalModel::Axis(axis) = model else {
            panic!("uniaxial parameter in a biaxial model");
        };
        axis as usize
    });

    let c = coefficient_rate(&rates, PrincipalModel::Axis(CrystalAxis::C));
    let cc = c[2][2];
    assert!(cc.norm() > 0.0);
    assert_eq!(c, diagonal([0.0.into(), 0.0.into(), cc]));
}

#[test]
fn constant_tensor_evaluates_its_entries_with_no_dispersion() {
    let epsilon = [
        [
            Complex64::new(4.0, 0.1),
            Complex64::new(0.3, 0.0),
            0.0.into(),
        ],
        [
            Complex64::new(0.3, 0.0),
            Complex64::new(5.0, 0.2),
            0.0.into(),
        ],
        [0.0.into(), 0.0.into(), Complex64::new(6.0, 0.0)],
    ];
    let mu = diagonal([Complex64::from(1.0); 3]);
    let model = ConstantTensor::new(epsilon, mu).unwrap();
    const { assert!(!<ConstantTensor<f64> as TensorMaterial>::PRINCIPAL_AXES) };

    assert_eq!(
        model.relative_permittivity_tensor::<_, Complex64>(Scalar(K0)),
        epsilon
    );
    assert_eq!(
        model.relative_permeability_tensor::<_, Complex64>(Scalar(K0)),
        mu
    );
    let k = Complex64::new(K0, 300.0);
    assert_eq!(
        model.relative_permittivity_tensor_complex::<_, Complex64>(Scalar(k)),
        epsilon
    );
    assert_eq!(
        model.relative_permeability_tensor_complex::<_, Complex64>(Scalar(k)),
        mu
    );
    let zero = [[Complex64::from(0.0); 3]; 3];
    assert_eq!(
        model.relative_permittivity_tensor_derivative::<_, Complex64>(
            Scalar(K0),
            DerivativeOrder::First
        ),
        zero
    );
    assert_eq!(
        model.relative_permittivity_tensor_complex_derivative::<_, Complex64>(
            Scalar(k),
            DerivativeOrder::Second
        ),
        zero
    );
}

#[test]
fn constant_tensor_rejects_non_finite_entries_and_zero_cc_components() {
    let identity = diagonal([Complex64::from(1.0); 3]);

    let mut epsilon = identity;
    epsilon[2][2] = 0.0.into();
    assert_eq!(
        ConstantTensor::new(epsilon, identity),
        Err(ConstantTensorError::ZeroCcComponent {
            quantity: "permittivity"
        })
    );
    assert_eq!(
        ConstantTensor::new(identity, epsilon),
        Err(ConstantTensorError::ZeroCcComponent {
            quantity: "permeability"
        })
    );

    let mut epsilon = identity;
    epsilon[0][1] = Complex64::new(0.0, f64::NAN);
    assert_eq!(
        ConstantTensor::new(epsilon, identity),
        Err(ConstantTensorError::NonFiniteEntry {
            quantity: "permittivity",
            row: CrystalAxis::A,
            column: CrystalAxis::B,
        })
    );
}

#[test]
fn constant_tensor_component_rates_are_unit_tensors() {
    let identity = diagonal([Complex64::from(1.0); 3]);
    let rates = ConstantTensor::new(identity, identity)
        .unwrap()
        .permittivity_parameter_rates();
    assert_eq!(rates.len(), 9);

    let ab = TensorParameter::Permittivity {
        row: CrystalAxis::A,
        column: CrystalAxis::B,
    };
    assert!(
        rates
            .iter()
            .any(|rate| rate.parameter == ab && rate.value == unit(0, 1))
    );
    for rate in rates {
        let TensorParameter::Permittivity { row, column } = rate.parameter else {
            panic!("unexpected parameter {:?}", rate.parameter);
        };
        assert_eq!(rate.value, unit(row as usize, column as usize));
    }
}
