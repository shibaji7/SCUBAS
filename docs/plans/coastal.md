# Plan: FDTD/FEM modelling of coastal effects

**Roadmap:** planned. **Stub:** `scubas.future.coastal.numerical`. **Extends:** `scubas.cables` (the end potentials). **Model pages:** [the DSTL model](../model/dstl.md) and [numerical coast models](../model/fdtd.md).

## Goal

A numerical solution of the electric field and the Earth potentials where a cable comes ashore, to check, and where it fails to improve, the end potentials $U_0$ and $U_1$ that
the transmission-line model gives with active terminations.

## What SCUBAS has today

- The coast effect enters only through the end potentials of the distributed-source transmission-line (DSTL) model and its active terminations; the land beyond the shore is an equivalent
  circuit (`TransmissionLine.term_params`, see [active terminations](../tutorial/active.md)). Nothing in SCUBAS solves a two- or three-dimensional problem.
- The model is 1-D in depth for each section, with lateral changes between sections only.

## What the papers say

- The mechanism: current in the sea is denser than on land, charge accumulates at the shore, the potential gradient raises the field on the resistive side and lowers it on the
  conductive side, and the cable end sees an Earth potential that partly cancels the induced emf (Boteler et al., 2024, Section 2.3, after Price, 1973).
- DSTL is equivalent to generalised thin-sheet modelling (Wang et al., 2023, as described in Boteler et al., 2024, Section 3.2).
- The segmented 1-D framework "does not explicitly resolve" transverse two-dimensional conductivity gradients or localised edge effects at abrupt bathymetric transitions, and SCUBAS is
  "not intended to replace" the three-dimensional electromagnetic models (Chakraborty et al., 2026, Sections 1 and 2.3).

From the wider literature (found by search, **not read in full**): Liu et al. (2018, *Space Weather and Space Climate*) and a two-dimensional finite-element study of the coast effect
(frequency, the shape of the seafloor boundary and the ocean depth matter); a regional three-dimensional study of the British Isles finding the coast effect on land almost everywhere;
a global 3-D study finding peak fields at least twice those of 1-D models in many coastal regions. Cite only after reading them.

## Open questions

1. FEM or FDTD? The problem is quasi-static (displacement currents are negligible at these frequencies, as the finite-element papers note), which suits a frequency-domain FEM; an FDTD
   scheme needs a time-domain diffusion formulation. This plan does not decide; the choice needs the authors' view and a short trial.
2. 2-D (a section across the coast) first, or 3-D? 2-D is cheap and tests the DSTL end potentials directly.
3. Dependencies: a mesh generator and a solver are heavy. They go behind an optional extra, never into the runtime requirements (ADR 0003).
4. What is the output that SCUBAS uses: a table of correction factors for $U_0$, $U_1$, or a replacement of the active termination?
5. Benchmarks: the vertical-boundary problem of Price (1973) and the coast cases of Wang et al. (2023).

## How it will be tested

A uniform half-space must reproduce the 1-D surface impedance (the half-space check already in the compatibility suite); two homogeneous quarter-spaces far from the interface must
reproduce their own 1-D fields; the DSTL end potential of a simple two-region case is the benchmark, with a tolerance set by the first trial rather than assumed.
