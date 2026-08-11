use ndarray::{Array1, Array2, array};
use num_complex::Complex64;

use lamina_core::{Material, MeromorphicMaterial};

use lamina_materials::Drude;

fn material() -> Drude<f64> {
    Drude::new(1.0, 2.0, 0.1).unwrap()
}

#[test]
fn real_vector_sampling_preserves_dimension() {
    let material = material();

    let input = array![1.0, 2.0, 3.0];

    let output: Array1<Complex64> = material.relative_permittivity(input);

    assert_eq!(output.shape(), &[3]);
}

#[test]
fn real_matrix_sampling_preserves_dimension() {
    let material = material();

    let input = array![[1.0, 2.0], [3.0, 4.0],];

    let output: Array2<Complex64> = material.relative_permittivity(input);

    assert_eq!(output.shape(), &[2, 2]);
}

#[test]
fn complex_vector_sampling_preserves_dimension() {
    let material = material();

    let input = array![Complex64::new(1.0, 0.1), Complex64::new(2.0, -0.2),];

    let output: Array1<Complex64> = material.relative_permittivity_complex(input);

    assert_eq!(output.shape(), &[2]);
}
