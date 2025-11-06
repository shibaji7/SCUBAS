<!-- 
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is under the MIT license found in the root directory LICENSE.md 
Everyone is permitted to copy and distribute verbatim copies of this license 
document.

This version of the MIT Public License incorporates the terms
and conditions of MIT General Public License.
-->
#### Example Modelling Example Codes
---
To illustrate the application of the modelling technique and to examine the voltages produced in typical Ocean/Earth conductivity structures we present results for three models representing 1) a shallow sea/continental shelf, 2) a shallow ocean, and 3) a deep ocean.

1. Shallow Sea/Continental Shelf
> This model consists of a shallow sea 600 km wide, with a depth of 100 m, with land at either end  as shown in figure (a). The seawater section has the same properties as the continental shelf in table listed in [network modeling](netmodel.md). This is a simplified model approximation to the North Sea.

![Case 1 – Continental shelf potential](figures/conduct_case1_continental_shelf.png "Potential profile along a 600 km continental-shelf cable")

2. Shallow Ocean
> This model consists of a deep ocean section 1800 km wide, with a depth of 1 km, with continental shelf 100 km wide and 100 m deep on either side as shown in figure (b). The total width from coast to coast is 2000 km. The parameters for each section are the same as those shown in table listed in [network modeling](netmodel.md), except the ocean section has a seawater depth of 1 km. This is a simplified model approximation of the Tasman Sea.

![Case 2 – Shallow ocean potential](figures/conduct_case2_shallow_ocean.png "Potential profile for the shallow-ocean composite cable")

3. Deep Ocean
> This model consists of a deep ocean section 7800 km wide, with a depth of 4 km, with continental shelf 100 km wide and 100 m deep on either side as shown in Figure (c). The total width from coast to coast is 8000 km. The parameters for each section are as shown in table listed in [network modeling](netmodel.md). This is a simplified model approximation of the Pacific Ocean.![Alt text](../figures/Case-Studies.png)

![Case 3 – Deep ocean potential](figures/conduct_case3_deep_ocean.png "Potential profile across the deep-ocean model with continental shelves on both ends")

##### Transmission Line Parameters 
For each model the layer thicknesses and resistivities from table in [network modeling](netmodel.md) are used to calculate the transmission line series impedance, $Z$, and parallel admittance, $Y$. These values are then used to calculate the transmission line propagation constant, $\gamma$, characteristic impedance, $Z_0$, and adjustment distance, $\frac{1}{\gamma}$ for each section. The values obtained are shown in table below. The 'Shallow Sea' parameters are used for the continental shelf sections in models 2 and 3 as well as for the seawater section in model 1.

|                    | Depth (m) | $Z$ ($\Omega/km$) | $Y$ (S/km) | $\gamma$ ($km^{-1}$) | $Z_0$ ($\Omega$) | Adj. Distance (km) |
| ------------------ | --------- | ----------------- | ---------- | -------------------- | ---------------- | ------------------ |
| Continental Shelf  | 100       | $7.5\times 10^{-1}$ | $5\times 10^{-6}$ | $1.9365\times 10^{-3}$ | $3.873\times 10^{2}$ | 516.4 |
| Shalow Ocean  | 1000       | $2.5\times 10^{-1}$ | $1.0\times 10^{-5}$ | $1.58\times 10^{-3}$ | $1.58\times 10^{2}$ | 632.5 |
| Deep Ocean  | 4000       | $7.14\times 10^{-1}$ | $1.0\times 10^{-5}$ | $8.45\times 10^{-4}$ | $8.45\times 10^{1}$ | 1183.2 |
| Land  | 0       | 3.0 | $5.0\times 10^{-6}$ | $3.873\times 10^{-3}$ | $7.746\times 10^{2}$ | 258.2 |


!!! Example
    ``` py
    from pathlib import Path

    import pandas as pd

    from scubas.datasets import PROFILES
    from scubas.cables import Cable, TransmissionLine
    from scubas.plotlib import cable_potential, update_rc_params

    figures_dir = Path("docs/tutorial/figures")
    figures_dir.mkdir(parents=True, exist_ok=True)

    update_rc_params(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"],
        },
        science=True,
    )

    def create_induced_field(millivolts_per_km: float) -> pd.DataFrame:
        return pd.DataFrame({"X": [millivolts_per_km]}, index=pd.RangeIndex(1, name="Time"))

    def build_transmission_line(
        sec_id: str,
        length_km: float,
        profile,
        efield: pd.DataFrame,
        terminations=None,
    ) -> TransmissionLine:
        tl = TransmissionLine(
            sec_id=sec_id,
            directed_length={"length_north": length_km},
            elec_params={"site": profile, "width": 1.0, "flim": [1e-6, 1.0]},
            active_termination=terminations,
        )
        tl.compute_eqv_pi_circuit(Efield=efield, components=["X"])
        return tl

    def simulate_case(title: str, slug: str, sections, ylim):
        tlines = [build_transmission_line(**section) for section in sections]
        cable = Cable(tlines, components=["X"])
        potentials, distances = cable._pot_along_cable_(timestamp=0)

        plot = cable_potential(potentials, distances, ylim=ylim)
        plot.axes.set_title(title)
        plot.axes.text(
            0.05,
            0.85,
            rf"$L_{\mathrm{total}}$ = {sum(section['length_km'] for section in sections):.0f} km",
            ha="left",
            va="center",
            transform=plot.axes.transAxes,
        )
        plot.figure.savefig(
            figures_dir / f"{slug}.png",
            dpi=300,
            bbox_inches="tight",
        )
        plot.figure.clf()

    efield_cs = create_induced_field(300.0)
    efield_ocean = create_induced_field(100.0)

    simulate_case(
        title="Case 1 – Continental Shelf",
        slug="conduct_case1_continental_shelf",
        sections=[
            {"sec_id": "CS", "length_km": 600.0, "profile": PROFILES.CS, "efield": efield_cs},
        ],
        ylim=(-200.0, 200.0),
    )

    simulate_case(
        title="Case 2 – Shallow Ocean",
        slug="conduct_case2_shallow_ocean",
        sections=[
            {
                "sec_id": "CS-west",
                "length_km": 100.0,
                "profile": PROFILES.CS,
                "efield": efield_cs,
                "terminations": {"left": {"site": PROFILES.LD, "width": 1.0}},
            },
            {
                "sec_id": "SO",
                "length_km": 1800.0,
                "profile": PROFILES.SO,
                "efield": efield_ocean,
            },
            {
                "sec_id": "CS-east",
                "length_km": 100.0,
                "profile": PROFILES.CS,
                "efield": efield_cs,
                "terminations": {"right": {"site": PROFILES.LD, "width": 1.0}},
            },
        ],
        ylim=(-60.0, 60.0),
    )

    simulate_case(
        title="Case 3 – Deep Ocean",
        slug="conduct_case3_deep_ocean",
        sections=[
            {
                "sec_id": "CS-west",
                "length_km": 100.0,
                "profile": PROFILES.CS,
                "efield": efield_cs,
                "terminations": {"left": {"site": PROFILES.LD, "width": 1.0}},
            },
            {
                "sec_id": "DO",
                "length_km": 7800.0,
                "profile": PROFILES.DO,
                "efield": efield_ocean,
            },
            {
                "sec_id": "CS-east",
                "length_km": 100.0,
                "profile": PROFILES.CS,
                "efield": efield_cs,
                "terminations": {"right": {"site": PROFILES.LD, "width": 1.0}},
            },
        ],
        ylim=(-120.0, 120.0),
    )
    ```

#### Boundary Conditions: Calculating Potentials
---
