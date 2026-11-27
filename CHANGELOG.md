# Changelog

The lamina crates (lamina-units, lamina-core, lamina-materials, lamina-mode) are
released in step until 1.0. Each entry names any shared convention it changes.

## Unreleased

## 0.2.0 (2026-11-27)

Added:

- Declared spectral singularities: `Drude`, `Lorentz` and `DrudeLorentz`
  declare every pole of ε and every zero (roots of its numerator polynomial,
  polished on ε), and `PrincipalAxes` declares its axis models' points on
  their entries when every axis model declares. Other models keep the
  undeclared default.
- `PrincipalAxes`: a tensor model diagonal in its crystal frame, with one scalar
  dispersion model per principal axis (`uniaxial`, `biaxial`). It implements each
  tensor capability shared by its three axis models.
- `ConstantTensor`: constant complex permittivity and permeability tensors,
  rejecting non-finite entries and a zero cc component.
- Tensor parameter rates (`TensorRate`) in the crystal frame: one rate per
  Sellmeier parameter, on the diagonal entries of the axes its model describes
  (`diag(∂ε_o, ∂ε_o, 0)` for an ordinary parameter, `diag(0, 0, ∂ε_e)` for an
  extraordinary one, one entry for a biaxial axis), and unit rates for each
  constant-tensor component.
- `Extrapolation::Reject`: a tabulated permittivity is NaN outside its table
  (the endpoints are inside it) and identical to `Clamp` and `Linear` inside.
  lamina-core refuses an oriented tensor region with such an entry with
  `ConstitutiveError::NonFinite` naming the region; a scalar region is not yet
  refused by name (the scattering-matrix backend returns a NaN response, the
  transfer-matrix backend a non-finite matrix error). `Extrapolation` gains a
  variant, so exhaustive matches on it need a new arm.

Fixed:

- `Tabulated` no longer panics at a NaN wavenumber (the permittivity is NaN) or
  on a table given in a non-contiguous array.

## 0.1.0

First release on crates.io.

Conventions fixed by this release:

- the canonical spectral coordinate is vacuum angular wavenumber `k₀ = 2π/λ₀`
  in `cm⁻¹`; in-plane wavenumbers use the same unit.
