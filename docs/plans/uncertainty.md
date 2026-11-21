# Plan: uncertainty quantification

**Roadmap:** planned. **Stub:** `scubas.future.uncertainty.propagation`. **Extends:** `scubas.conductivity` and `scubas.cables`. **Model pages:** [the SCUBAS model](../model/scubas.md) and [the DSTL model](../model/dstl.md).

## Goal

Give every voltage SCUBAS computes an honest range: propagate what is uncertain in the inputs to `Vt`, `U0` and `U1`, and say which input matters most (sensitivity).

## What SCUBAS has today

- `ConductivityProfile.compile_mcmc_bined_profiles` draws random realisations of binned LITHO1.0 profiles (defaults: 100 draws, seed 0, end-bin seawater thickness between
  0.01 and 1 km; see the [parameter list](../dev/parameters.md)). It produces profiles, not a distribution of voltages.
- The [uncertainty tutorial](../tutorial/uq.md) describes the intended Monte Carlo scheme (N sections times V virtual points) and says that no error estimate exists yet for the input
  magnetic field.
- Nothing computes a spread of `Vt`, and nothing computes a sensitivity index.

## What the papers say

Sources of uncertainty named in the SCUBAS studies (all read in full or in the relevant sections):

- The magnetic field along the route is not measured over the ocean; the nearest observatory stands in for it, and the choice matters when the disturbance is not uniform
  (Boteler et al., 2024, Section 5; Chakraborty et al., 2026, Sections 6 and 7).
- The section length relative to the adjustment distance $1/\gamma$ sets what the model resolves; sections shorter than a few hundred kilometres strain the thin-wire assumption
  (Chakraborty et al., 2026, Section 2.3).
- The 1-D Earth profile of each section (LITHO1.0, ETOPO1 depth) is itself a model (Boteler et al., 2024, Section 4.1).
- Digitised magnetograms and digitised voltage curves carry extraction error (Chakraborty et al., 2026, Sections 3 and 7).

From the wider literature (found by search, **not read in full**): Liu et al. (2020, IEEE) treat Earth conductivity and substation resistances as uncertain inputs of the GIC in a
power grid, with polynomial chaos and Sobol indices checked against Monte Carlo; Love and co-authors (USGS, 2018 and 2021) show that 1-D conductivity applied over wide regions and
sparse magnetotelluric sampling limit geoelectric hazard estimates. Verify these before citing them in a paper.

## Open questions

1. Which inputs first: conductivity layers, section boundaries, observatory choice, B-field error? (The papers suggest all four; the cheapest to start with are the layers.)
2. Method: plain Monte Carlo (simple, any model, slow), or polynomial chaos / surrogate (fewer runs, needs smooth response)? SCUBAS runs are cheap enough that Monte Carlo is
   probably the baseline either way.
3. What distribution for each input, and where does it come from (LITHO1.0 spread, measurement error of an observatory)? Each distribution is a [parameter with a status](../dev/parameters.md);
   an unknown one must raise, not be guessed.
4. What is reported: percentiles of `Vt(t)`, of its peak, of `U0 - U1`? How is the existing `compile_mcmc_bined_profiles` seed handled (it sets numpy's global seed)?

## How it will be tested

Against closed forms, not snapshots. For an electrically short section $V_t = LE/2$ is linear in $E$, so a Gaussian $E$ of known mean and variance gives a known mean and variance of
the voltage; for an electrically long one the end potential is $-E/\gamma$. A fixed seed gives identical output; the existing `compile_mcmc_bined_profiles` is not changed.
