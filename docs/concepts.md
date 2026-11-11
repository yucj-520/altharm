---
layout: default
title: "Concepts"
nav_order: 4
---

# Concepts

## Harmonics of the spin splitting

On a circle of radius k₀ around an axis, the spin splitting of one band is a periodic
function of the angle θ and can be written as

```
ΔE(θ) = Σ_l [ a_l cos(lθ) + b_l sin(lθ) ]
```

`altharm` fits the coefficients a_l and b_l by least squares and reports the amplitude
√(a_l² + b_l²) of each l. For a splitting that varies in three dimensions the same is
done on a sphere with real spherical harmonics.

## Allowed and forbidden harmonics

The magnetic symmetry of the crystal fixes which l can be non-zero. `altharm` derives this
numerically: it takes the symmetry operations of the structure from spglib, labels each
one by whether it keeps or exchanges the two magnetic sublattices, and tests every
harmonic against them. No table is consulted, so strained or distorted cells are handled
in the same way.

Example, Cr₂Se₂As at k_z = 0.25 b₃: the three-fold axis keeps the sublattices, so l must
be a multiple of 3; the horizontal mirror exchanges them, so l must be odd. The allowed
set is l = 3, 9, 15, …. In the plane k_z = 0 nothing is allowed.

A compensated antiferromagnet with an inversion centre between the sublattices gives the
empty set: it is not an altermagnet.

## The noise floor

All harmonics up to `lmax` are fitted, not only the allowed ones. A forbidden harmonic
must vanish exactly, so any amplitude found there is numerical noise of the calculation,
measured in the same way as the signal. This is the noise floor.

- An allowed harmonic is resolved only if it stands well above the floor.
- The floor answers whether `lmax` is high enough, whether the circle is sampled densely
  enough and whether a given harmonic is resolved.
- The floor is not the error bar of a coefficient. See
  [How far to trust the numbers](reliability.md).

## d-, g- and i-wave labels

The usual labels count nodal surfaces, not the index l. The two agree only when there is
no horizontal nodal plane.

| Material | Label | Leading l on the circle | Nodal surfaces |
|---|---|---|---|
| CoF₂ | d-wave | 2 | two vertical planes |
| MnTe | g-wave | 3 | three vertical planes and the plane k_z = 0 |

Reading "g-wave" as "l = 4" would make the correct MnTe result look wrong.

## Why directly from DFT

Harmonics are usually taken from a Wannier model fitted to DFT. A hopping of range R
contributes to harmonic l on a circle of radius k₀ with a weight given by the Bessel
function J_l(k₀R), which is very small once l exceeds k₀R. A Wannier model with a finite
hopping range therefore suppresses high harmonics, whether they are physical or not.
Plane-wave DFT has no hopping range and no such suppression, so it can decide whether a
high harmonic is real. In exchange it has numerical noise, which the noise floor measures.
