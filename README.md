# lamina-materials

Optical material models for planar multilayer calculations.

`lamina-materials` implements the material capability traits defined by
`lamina-core`, so a stack can evaluate permittivity and permeability, and their
spectral derivatives, on the real axis or in the complex plane.

## What it provides

- **Analytic models:** `Drude`, `Lorentz`, `DrudeLorentz`, `Sellmeier` and
  `Cauchy`. They are meromorphic in the spectral coordinate, so they support
  complex-frequency evaluation and analytic derivatives.
- **Tabulated data:** `Tabulated` interpolates measured optical constants on the
  real axis only. A real-axis interpolation defines no physical complex
  continuation, so tabulated materials don't implement the complex-plane traits.
- **A small curated set of named parameterisations** in `database`: vacuum, air,
  fused silica, silicon and germanium (Frey 2006), Rakić 1998 metals and polar
  dielectrics. Each record states its reference, parameter convention, valid
  spectral interval and units. It is a convenience for examples and tests, not
  an optical-constants database.

Material evaluation uses vacuum angular wavenumber `k₀ = 2π/λ₀` in `cm⁻¹` as its
canonical spectral coordinate. Constructors accept typed spectral quantities
from `lamina-units` and convert them before validation.

## Boundary

| lamina-materials owns | lamina-materials never owns |
| --- | --- |
| Constitutive models and their validation | Layer stacks and matrix methods (`lamina-core`) |
| Spectral derivatives of ε and μ | Units and spectral quantities (`lamina-units`) |
| Validity domains and literature identity of named records | Mode searches and fields (`lamina-mode`) |
| A small set of documented reference parameterisations | Comprehensive material databases, foundry or process data |

Consumers supply the parameters of the materials they model. A change that
turns this crate into a material database, or makes it depend on a solver,
crosses the boundary.

## Family

```text
lamina-units       typed quantities and units
lamina-core        stacks, transfer and scattering matrices, fields
lamina-materials   material models implementing lamina-core's material traits (this crate)
lamina-mode        modes of planar multilayers, with numerical evidence
```

The crates are versioned together until 1.0; each changelog names the
conventions a release changes.
