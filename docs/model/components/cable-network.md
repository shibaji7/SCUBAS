# Component: cable line and network

:::{note}
Draft stub. The writer filled in what the code and the papers say; the authors complete the physics below. Section 3 ("To fill in") is yours.
:::

## 1. In the code today

- `scubas.cables`: `CableSection` (geometry: north and east lengths from end points), `TransmissionLine` (the line parameters $Z$, $Y$, $\gamma$, $Z_0$, the equivalent-$\pi$ circuit,
  active terminations), `Cable` (joins the sections, solves the admittance matrix, returns the potentials in `tot_params` and writes `est_cable_props.json` and `sim-params.csv`).
- Equations: [theory](../../tutorial/theory.md), [electrically long and short cables](../../tutorial/elsc.md), [active terminations](../../tutorial/active.md), [network modelling](../../tutorial/netmodel.md).

## 2. In the papers

- Per unit length $Z$ is the series impedance of the seawater and sediment layers and $Y$ the shunt admittance through the crust and mantle lithosphere (Boteler et al., 2024, Section 3.2).
- A check anyone can repeat: the "Ocean 2" section of the TAT-8 route has $Z = 5.56\times10^{-2}$ ohm/km and $Y = 5.88\times10^{-6}$ S/km, so
  $\gamma = \sqrt{ZY} = 5.72\times10^{-4}$ /km, $Z_0 = \sqrt{Z/Y} = 97.2$ ohm and $1/\gamma \approx 1750$ km, as in Table 3 of Boteler et al. (2024).
- Active terminations represent the land at each end (Boteler et al., 2023); the potential inside a section follows Equation 15 of Boteler et al. (2024).
- Limit: the thin-wire assumption of DSTL; sections shorter than a few hundred kilometres strain it (Chakraborty et al., 2026, Section 2.3).

## 3. To fill in (authors)

- The ground (electrode) model at the ends and the power-feed configuration of a real system.
- Branching units, such as the TAT-8 branch to France and the UK (the model uses the main trunk, Chakraborty et al., 2026, Section 2.2).
- The conventions for sign and direction of the voltage.

Back to [the model](../index.md).
