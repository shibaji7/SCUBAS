# Plan: auroral electrojet design

**Roadmap:** planned. **Stub:** `scubas.future.sources.electrojet`. **Extends:** `scubas.sources`. **Model page:** [Auroral electrojet currents](../model/electrojet.md).

## Goal

A modelled ionospheric source (an auroral electrojet with a position, width, height and time-varying strength) that gives the magnetic field, and from it the electric field, along a
cable route where no magnetometer exists. It removes the main limit the papers name: the missing magnetic data over the ocean.

## What SCUBAS has today

- `scubas.sources.current_sources.Current` is a single line current at height `h` with a Cauchy-distributed width `a`, over a `Site` (a layered Earth). `compute_B_field` returns the
  magnetic field components and an electric field along the current for given distances, in the frequency domain. It is documented in [Current source modelling](../dev/sources.md).
- It is not connected to `Cable`; the cable is driven by observatory records only. The module has [known defects](../dev/testing.md#known-defects) that are pinned and not fixed.
- A design needs more than one filament, a finite length, a time dependence and a way to fit the parameters to observations; none exists.

## What the papers say

- Boteler et al. (2024, Section 5 and Figure 10) infer electrojet positions from the magnetic variations at the three observatories for two times of the 1989 storm, and show that an
  electrojet over the whole route gives the largest voltage while opposite signs at the two ends cancel. They state that future work will examine induction by modelled ionospheric and
  magnetospheric current sources.
- Chakraborty et al. (2026, Sections 6 and 7) name small-scale ionospheric currents and sparse magnetometers as the main limit, and say SCUBAS does not interpolate magnetic structure.

From the wider literature (found by search, **not read in full**): Pirjola and Viljanen (1998, *Annales Geophysicae* 16, 1434) extend the complex image method to an electrojet of finite
length, a horizontal line current with vertical currents at its ends over a layered Earth, and report that it is accurate and much faster than the exact approach; they build on
Thomson and Weaver (1975) and Boteler and Pirjola (1998). Pirjola and Boteler (2002, *Radio Science*) review line-current methods. Verify the equations in the originals before use.

## Open questions

1. Geometry: infinite line, finite line with end currents (Pirjola and Viljanen), or a sheet? One electrojet or an eastward and a westward one?
2. Parameters and their provenance: height, width, position, strength, time dependence, each a [parameter with a status](../dev/parameters.md); fitted from observatory data (an inverse
   problem, which brings in [uncertainty](uncertainty.md)) or given?
3. Earth response: the existing `Current` uses the skin depth of one layer; the complex image method represents the layered Earth more completely. Which does the design need?
4. Validation: the 1989 event and the electrojet positions of Boteler et al. (2024, Figure 10) are the natural first case.
5. Relation to `Current`: extend it, or add a new class beside it and leave 0.1.9 untouched (the compatibility suite pins the existing one).

## How it will be tested

Against closed forms: the field of an infinite line current at height $h$ and horizontal distance $x$ has the horizontal component $\mu_0 I h / (2\pi(h^2+x^2))$ and the vertical
component $\mu_0 I x / (2\pi(h^2+x^2))$, which the existing code reproduces for the external part; a finite line must tend to it as its length grows; and the total (external plus
internal) horizontal component must tend to twice the external part when the skin depth goes to zero, the limit of a perfectly conducting Earth.
