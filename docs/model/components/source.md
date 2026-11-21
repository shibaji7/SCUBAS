# Component: source (auroral electrojet and the magnetic record)

:::{note}
Draft stub. The writer filled in what the code and the papers say; the authors complete the physics below. Section 3 ("To fill in") is yours.
:::

## 1. In the code today

- The model is driven by a **recorded magnetic field**: IAGA-2002 `.txt` / `.min` files or CSV with a `Date` column and `X, Y, Z, F` (`OceanModel.read_iaga`, `read_Bfield_data`; gaps in a
  component are filled by a second-order spline). One record drives each section; the user chooses which.
- `scubas.sources.current_sources.Current` models a line current at height `h` with a Cauchy-distributed width over a layered Earth and returns the magnetic and electric field
  at given distances. It is **not** connected to `Cable`; see [Current source modelling](../../dev/sources.md).
- Nothing models an auroral electrojet as a source for a cable. That is the [auroral electrojet plan](../../plans/electrojet.md).

## 2. In the papers

- The name of the model points at the source: auroral streams. In practice the 2024 and 2026 studies use observatory records only: Fredericksburg, St John's and Hartland for TAT-8,
  Eskdalemuir for TAT-1 (Boteler et al., 2024, Table 1; Chakraborty et al., 2026, Tables 1 and 2).
- The electrojet shows only through the records: an eastward electrojet in the afternoon at one end and a westward one at the other give opposite fields and a small voltage; a
  westward electrojet over the whole route gives a large one. Electrojet positions were estimated from the magnetic variations (Boteler et al., 2024, Section 5, Figure 10).
- Magnetometer coverage over the ocean is the main limit; small-scale current structure is not captured (Chakraborty et al., 2026, Sections 6 and 7).

## 3. To fill in (authors)

- What the model should say about the electrojet: geometry, altitude, width, relation to magnetic latitude.
- How a source field enters the transfer function (uniform field, or a finite source; see the assumption flagged **to confirm** on the model page).
- How the 1958 and 1989 events differ in source terms (AE of about 2500 nT in 1958 and 1800 nT in 1989, and the magnetic latitude of the two routes, Chakraborty et al., 2026, Section 6).

Back to [the model](../index.md).
