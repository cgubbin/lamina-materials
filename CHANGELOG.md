# Changelog

The lamina crates (lamina-units, lamina-core, lamina-materials, lamina-mode) are
released in step until 1.0. Each entry names any shared convention it changes.

## Unreleased

Added:

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

Fixed:

- `Tabulated` no longer panics at a NaN wavenumber (the permittivity is NaN) or
  on a table given in a non-contiguous array.

## 0.1.0 (unreleased)

First release on crates.io.

Conventions fixed by this release:

- the canonical spectral coordinate is vacuum angular wavenumber `k₀ = 2π/λ₀`
  in `cm⁻¹`; in-plane wavenumbers use the same unit.
