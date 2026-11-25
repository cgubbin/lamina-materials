# lamina-materials

Optical material models for planar multilayer calculations.

`lamina-materials` implements the material capability traits defined by
`lamina-core`, so a stack can evaluate permittivity and permeability, and their
spectral derivatives, on the real axis or in the complex plane.

## What it provides

- **Analytic models:** `Drude`, `Lorentz`, `DrudeLorentz`, `Sellmeier` and
  `Cauchy`. They are meromorphic in the spectral coordinate, so they support
  complex-frequency evaluation and analytic derivatives.
- **Anisotropic models:** `PrincipalAxes` describes a birefringent medium by one
  scalar model per principal axis, `ε = diag(ε_a, ε_b, ε_c)` in the crystal frame,
  with `uniaxial` (ordinary on a and b, extraordinary on the optic axis c) and
  `biaxial` constructors. It has each capability that all three axis models have,
  so a tabulated axis makes it real-axis only. `ConstantTensor` holds constant
  complex 3×3 permittivity and permeability tensors. Both report the rates of
  their permittivity tensor with respect to their own parameters. Orientation is
  applied by a rotation from `lamina-core`; the crate holds no crystal data.
- **Tabulated data:** `Tabulated` interpolates measured optical constants on the
  real axis only. A real-axis interpolation defines no physical complex
  continuation, so tabulated materials don't implement the complex-plane traits.
  Outside the table (endpoints are inside it) the permittivity is the nearest
  endpoint value (`Clamp`), the nearest segment continued (`Linear`), or NaN
  (`Reject`). `Reject` is the safe choice when the table should cover every
  evaluated wavenumber: lamina-core refuses an oriented tensor region with a
  non-finite entry with a typed error naming the region.
- **A small curated set of named parameterisations** in `database`: vacuum, air,
  fused silica, silicon and germanium (Frey 2006), Rakić 1998 metals and polar
  dielectrics. Each record states its reference, parameter convention, valid
  spectral interval and units, for examples and tests.

Material evaluation uses vacuum angular wavenumber `k₀ = 2π/λ₀` in `cm⁻¹` as its
canonical spectral coordinate. Constructors accept typed spectral quantities
from `lamina-units` and convert them before validation.

## Scope

lamina-materials provides constitutive models for planar multilayers: Drude, Lorentz, Drude–Lorentz,
Sellmeier, Cauchy and tabulated optical constants. Each model evaluates its permittivity and permeability, and
their spectral derivatives, on the real axis or, where the model defines one, in the complex plane, with
validation and a stated validity domain. Anisotropic media are built from these as principal-axis tensors in
the crystal frame, or given as constant tensors. The Sellmeier model and the tensor models also give the rates of
their permittivity with respect to their own parameters.

## Family

```text
lamina-units       typed quantities and units
lamina-core        stacks, transfer and scattering matrices, fields
lamina-materials   material models implementing lamina-core's material traits (this crate)
lamina-mode        modes of planar multilayers, with numerical evidence
```

The crates are versioned together until 1.0; each changelog names the
conventions a release changes.

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or
  <https://www.apache.org/licenses/LICENSE-2.0>)
- MIT license ([LICENSE-MIT](LICENSE-MIT) or <https://opensource.org/licenses/MIT>)

at your option.

### Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.
