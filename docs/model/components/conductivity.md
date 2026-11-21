# Component: Earth conductivity

:::{note}
Draft stub. The writer filled in what the code and the papers say; the authors complete the physics below. Section 3 ("To fill in") is yours.
:::

## 1. In the code today

- `scubas.datasets` holds `Layer`, `Site` and `PROFILES`, 19 ready-made profiles (continental shelf, shallow and deep ocean, land, the nine TAT-8 sections, and others).
- `scubas.conductivity.ConductivityProfile` builds a layered profile from LITHO1.0 (downloaded on first use; see `SECURITY.md`): at one point, for binned points, or as Monte Carlo
  realisations. Default resistivities (ohm m): seawater 0.3, sediments 3, crust 3000, mantle lithosphere 1000, upper mantle 100, transition zone 10, lower mantle 1
  ([parameter list](../../dev/parameters.md)).
- The layer table of a ready profile is in [Network modelling](../../tutorial/netmodel.md).

## 2. In the papers

- The depth of each section comes from ETOPO1 and the layers below from LITHO1.0 (Pasyanos et al., 2014); the TAT-8 profiles are in Table 2 of Boteler et al. (2024) and the TAT-1 ones in
  Table 3 of Chakraborty et al. (2026). The mid-Atlantic ridge is its own region because of the thinner mantle lithosphere.
- Lateral change is handled by a different 1-D profile for each section, not by a 3-D structure.
- **Seawater conductivity is quoted differently:** 4 S/m in the text of Boteler et al. (2024, Section 3.1), 3 S/m in Chakraborty et al. (2026, Section 2.3), and 0.3 ohm m (3.3 S/m) in the
  layer tables and as the default of the code. The authors should say which value the published runs used and which one the documentation should show.

## 3. To fill in (authors)

- The layer model: what each layer represents, and why those resistivities.
- Where the conductivity of land at the ends comes from, and how the profiles are chosen for new cables.
- The seawater value (above), and which parameters have a measured source.

Back to [the model](../index.md).
