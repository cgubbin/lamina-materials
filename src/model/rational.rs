//! Poles and zeros of permittivities that are rational in k0.
//!
//! A model of the form `ε(k) = ε∞ + Σᵢ nᵢ / Dᵢ(k)` with monic quadratics `Dᵢ`
//! has its poles at the roots of the `Dᵢ` and its zeros at the roots of
//! `N = ε∞ Πᵢ Dᵢ + Σᵢ nᵢ Πₗ≠ᵢ Dₗ`. Terms sharing a denominator are merged
//! and vanishing terms dropped first, so `N` shares no root with `Πᵢ Dᵢ`.

use lamina_core::ComplexScalar;
use num_traits::Float;

/// One term `n / (k² + b k + c)`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Term<C> {
    pub numerator: C,
    /// `[c, b]`: the constant and linear coefficients of the monic quadratic.
    pub denominator: [C; 2],
}

/// Poles and zeros of `ε∞ + Σ terms`, the zeros polished by Newton on `ε`
/// itself; `None` when ε vanishes identically or the roots do not converge.
pub(crate) fn poles_and_zeros<C, R, E>(
    epsilon_infinity: R,
    terms: &[Term<C>],
    epsilon: E,
) -> Option<(Vec<C>, Vec<C>)>
where
    C: ComplexScalar<RealField = R>,
    R: Float,
    E: Fn(C) -> (C, C),
{
    // Merge equal denominators, drop vanishing numerators.
    let mut merged: Vec<Term<C>> = Vec::new();
    for term in terms {
        match merged
            .iter_mut()
            .find(|m| m.denominator == term.denominator)
        {
            Some(m) => m.numerator += term.numerator,
            None => merged.push(*term),
        }
    }
    merged.retain(|t| !t.numerator.is_zero());

    // Work in u = k / s, with s the scale of the coefficients.
    let one = R::one();
    let mut scale = R::zero();
    for t in &merged {
        scale = scale
            .max(t.denominator[0].modulus().sqrt())
            .max(t.denominator[1].modulus())
            .max(t.numerator.modulus().sqrt());
    }
    if !(scale.is_finite() && scale > R::zero()) {
        scale = one;
    }
    let s = C::from_real(scale);
    let quadratic = |t: &Term<C>| vec![t.denominator[0] / (s * s), t.denominator[1] / s, C::one()];

    let mut poles: Vec<C> = Vec::new();
    for t in &merged {
        let [c, b] = t.denominator;
        let two = C::from_real(one + one);
        let root = (b * b - c * two * two).sqrt();
        for pole in [(-b + root) / two, (-b - root) / two] {
            if !poles.contains(&pole) {
                poles.push(pole);
            }
        }
    }

    // N(u) = ε∞ Π Dᵢ + Σ (nᵢ/s²) Πₗ≠ᵢ Dₗ, in u.
    let mut numerator = vec![C::from_real(epsilon_infinity)];
    for t in &merged {
        numerator = multiply(&numerator, &quadratic(t));
    }
    for (i, t) in merged.iter().enumerate() {
        let mut product = vec![t.numerator / (s * s)];
        for (l, other) in merged.iter().enumerate() {
            if l != i {
                product = multiply(&product, &quadratic(other));
            }
        }
        add_into(&mut numerator, &product);
    }
    let magnitude = numerator.iter().fold(R::zero(), |m, c| m.max(c.modulus()));
    let negligible = R::epsilon() * R::from(64.0)? * magnitude;
    while numerator.last().is_some_and(|c| c.modulus() <= negligible) {
        numerator.pop();
    }
    if numerator.is_empty() {
        return None;
    }
    let zeros = aberth(&numerator)?
        .into_iter()
        .map(|u| polish(u * s, &epsilon))
        .collect::<Option<Vec<_>>>()?;
    Some((poles, zeros))
}

fn multiply<C: ComplexScalar>(a: &[C], b: &[C]) -> Vec<C> {
    let mut product = vec![C::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            product[i + j] += *x * *y;
        }
    }
    product
}

fn add_into<C: ComplexScalar>(sum: &mut Vec<C>, term: &[C]) {
    if sum.len() < term.len() {
        sum.resize(term.len(), C::zero());
    }
    for (s, t) in sum.iter_mut().zip(term) {
        *s += *t;
    }
}

fn finite<C: ComplexScalar>(z: C) -> bool
where
    C::RealField: Float,
{
    z.real().is_finite() && z.imaginary().is_finite()
}

/// The value and derivative of a polynomial with coefficients low to high.
fn horner<C: ComplexScalar>(coefficients: &[C], z: C) -> (C, C) {
    let mut value = C::zero();
    let mut derivative = C::zero();
    for c in coefficients.iter().rev() {
        derivative = derivative * z + value;
        value = value * z + *c;
    }
    (value, derivative)
}

/// Every root of a polynomial (coefficients low to high, leading non-zero)
/// by Aberth–Ehrlich iteration; `None` if it does not converge.
fn aberth<C, R>(coefficients: &[C]) -> Option<Vec<C>>
where
    C: ComplexScalar<RealField = R>,
    R: Float,
{
    let degree = coefficients.len() - 1;
    if degree == 0 {
        return Some(vec![]);
    }
    let leading = coefficients[degree];
    // Cauchy bound on the root moduli.
    let radius = coefficients[..degree]
        .iter()
        .fold(R::zero(), |m, c| m.max((*c / leading).modulus()))
        + R::one();
    let n = R::from(degree)?;
    let mut roots: Vec<C> = (0..degree)
        .map(|i| {
            let angle = R::from(std::f64::consts::TAU)? * (R::from(i)? + R::from(0.25)?) / n;
            Some(C::from_parts(radius * angle.cos(), radius * angle.sin()))
        })
        .collect::<Option<_>>()?;
    let tolerance = R::epsilon() * R::from(16.0)?;
    for _ in 0..500 {
        let mut largest = R::zero();
        for i in 0..degree {
            let z = roots[i];
            let (value, derivative) = horner(coefficients, z);
            if value.is_zero() {
                continue;
            }
            let ratio = value / derivative;
            let repulsion = roots
                .iter()
                .enumerate()
                .filter(|&(j, _)| j != i)
                .fold(C::zero(), |sum, (_, w)| sum + C::one() / (z - *w));
            let step = ratio / (C::one() - ratio * repulsion);
            if !finite(step) {
                return None;
            }
            roots[i] = z - step;
            largest = largest.max(step.modulus() / roots[i].modulus().max(R::one()));
        }
        if largest <= tolerance {
            return Some(roots);
        }
    }
    None
}

/// Two Newton steps on `ε` from `k`; `None` if a step is not finite.
fn polish<C, E>(mut k: C, epsilon: &E) -> Option<C>
where
    C: ComplexScalar,
    C::RealField: Float,
    E: Fn(C) -> (C, C),
{
    for _ in 0..2 {
        let (value, derivative) = epsilon(k);
        if value.is_zero() {
            break;
        }
        let step = value / derivative;
        if !finite(step) {
            return None;
        }
        k -= step;
    }
    Some(k)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_complex::Complex;

    #[test]
    fn aberth_finds_the_roots_of_a_cubic() {
        let c = |re: f64, im: f64| Complex::new(re, im);
        let roots = [c(1.0, 2.0), c(-3.0, 0.5), c(0.25, -1.0)];
        let mut coefficients = vec![c(1.0, 0.0)];
        for r in roots {
            coefficients = multiply(&coefficients, &[-r, c(1.0, 0.0)]);
        }
        let found: Vec<Complex<f64>> = aberth(&coefficients).unwrap();
        for r in roots {
            assert!(found.iter().any(|f| (f - r).norm() < 1e-12), "{found:?}");
        }
    }
}
