# Plan: automatic segmenter

**Roadmap:** planned. **Stub:** `scubas.future.segmentation.auto_segmenter`. **Extends:** [cable line and network](../model/components/cable-network.md),
[Earth conductivity](../model/components/conductivity.md).

## Goal

Turn a cable route (a list of latitude and longitude points) and a bathymetry into the list of sections SCUBAS needs: for each section its ends, depth, Earth model and the
magnetic observatory that drives it, by a rule that can be stated and repeated.

## What SCUBAS has today

- `CableSection`, `TransmissionLine` and `Cable` take sections that the user has already chosen; `CableSection.compute_lengths` turns end points into north and east lengths
  (WGS84-based constants, see the [parameter list](../dev/parameters.md)).
- `ConductivityProfile` can build a profile at a point or for binned points from LITHO1.0, and `PROFILES` holds ready-made sections (continental shelf, deep ocean, Mid-Atlantic Ridge).
- Nothing chooses the section boundaries or assigns observatories.

## What the papers say

- The boundaries are chosen by hand, "primarily" from comparable water depth and section length; the seafloor is approximated by flat steps (Chakraborty et al., 2026, Section 2.1;
  Boteler et al., 2024, Section 4). The TAT-8 route is nine sections (Boteler et al., 2024, Table 1), TAT-1 eight (Chakraborty et al., 2026, Table 1).
- Depths come from ETOPO1 and the Earth model of each section from LITHO1.0 (both papers).
- Each section is assigned the observatory closest to its midpoint (Boteler et al., 2024, Section 4; Chakraborty et al., 2026, Section 2.2).
- A lower bound on the section length: below a few hundred kilometres the thin-wire approximation of DSTL is strained, while more sections follow the route better
  (Chakraborty et al., 2026, Section 2.3). The adjustment distance of the TAT-8 sections is 490 to 1750 km (Boteler et al., 2024, Table 3).

## Open questions

1. The rule: merge consecutive bathymetry points while the depth stays within a tolerance? A change-point search? Both give a different answer on a ridge such as the Mid-Atlantic Ridge,
   which the papers treat as its own region.
2. The minimum length: a fixed number of kilometres, or a multiple of the section's own $1/\gamma$ (which depends on the depth, so on the segmentation itself)?
3. Bathymetry source and resolution (ETOPO1 in the papers); a download would be new network access and must stay out of the tests and behind an extra.
4. Observatory list: from the user, or from a catalogue (and which)? Nearest to the midpoint is the papers' rule; it ignores data gaps and storm-time coherence.
5. The output: a ready list of `TransmissionLine` arguments, so that it feeds `Cable` without changing it.

## How it will be tested

A synthetic bathymetry with known steps must return exactly those breakpoints; a constant depth returns one section; the result does not depend on the order of the points; and
the section lengths add up to the great-circle length of the route within the tolerance of `compute_lengths`.
