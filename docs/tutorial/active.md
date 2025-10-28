<!-- 
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is under the MIT license found in the root directory LICENSE.md 
Everyone is permitted to copy and distribute verbatim copies of this license 
document.

This version of the MIT Public License incorporates the terms
and conditions of MIT General Public License.
-->
#### Transmission Line with 'Active' Terminations
---

When studying electromagnetic induction in a transmission line with multiple sections, it can be convenient to simplify the analysis by considering each section individually. To do this, we can represent the sections on either side of a specific section using their Thevenin equivalent circuits.

The Thevenin equivalent circuit is a simplified model that replaces a complex network of components with a single voltage source and a single equivalent impedance. It allows us to analyze a portion of a larger circuit as if it were a self-contained system.In the context of electromagnetic induction in a transmission line, the Thevenin equivalent circuit represents the behavior of the adjacent sections on either side of a particular section. It captures the combined effect of those sections, providing a simplified representation that can be easily analyzed.

The figure here likely illustrates the transmission line with multiple sections, and it shows how each section can be represented by its Thevenin equivalent circuit. This approach allows for a more manageable analysis of the electromagnetic induction phenomenon. By breaking down the transmission line into smaller sections and representing them with their Thevenin equivalent circuits, we can study the induction process in a step-by-step manner. This simplification enables us to analyze the impact of each section on the overall behavior of the transmission line, making the analysis more tractable and facilitating a better understanding of the electromagnetic phenomena involved. Transmission Line with distributed series impedance, $Z$, parallel admittance, $Y$, and voltage sources representing the electric field, $E$, with 'active' terminations at each end represented by Thevenin equivalent circuits.
![Alt text](../figures/AT-Cable-CS.png)


In this case, electromagnetic induction in the transmission line in 'Active Termination' figure can be represented by the circuit shown in following figure. ![Alt text](../figures/Short-Cable-AT.png)

The current along the transmission line is then given by:

$$
I=\frac{V_1+EL+V_2}{Z_1+ZL+Z_2}
$$

The potential at the left end of the transmission line is then:

$$
U_1=V_1-IZ_1
$$

And the potential at the right end of the transmission line is:

$$
U_2=IZ_2-V_2
$$

!!! Example
    ``` py
    # Core scientific stack
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    # SCUBAS dependencies
    from scubas.datasets import PROFILES
    from scubas.cables import Cable, TransmissionLine
    from scubas.models import OceanModel
    from scubas.plotlib import plot_transfer_function, update_rc_params

    # Configure Matplotlib using the helper (gracefully falls back if SciencePlots is missing)
    update_rc_params(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"],
        },
        science=True,
    )

    # Render the ocean model transfer function for the deep ocean profile
    ocean_site = PROFILES.DO_3
    ocean_model = OceanModel(ocean_site)
    tf_artifacts = plot_transfer_function(ocean_model.get_TFs())
    tf_artifacts.figure.suptitle("Deep Ocean Transfer Function")

    ####################################################################
    # Simulating the case: Induced electric field 0.3 V/km on a
    # shallow continental shelf with depth 100 m, length 600 km
    ####################################################################

    length = 600.0
    induced_field = pd.DataFrame({"X": np.array([300.0])}, index=pd.RangeIndex(1, name="Time"))

    tl_active = TransmissionLine(
        sec_id="CS-active",
        directed_length={"length_north": length},
        elec_params={"site": PROFILES.CS, "width": 1.0, "flim": [1e-6, 1.0]},
        active_termination={"right": {"site": PROFILES.LD, "width": 1.0}, "left": {"site": PROFILES.LD, "width": 1.0}},
    )
    tl_passive = TransmissionLine(
        sec_id="CS-passive",
        directed_length={"length_north": length},
        elec_params={"site": PROFILES.CS, "width": 1.0, "flim": [1e-6, 1.0]},
    )

    tl_active.compute_eqv_pi_circuit(Efield=induced_field, components=["X"])
    tl_passive.compute_eqv_pi_circuit(Efield=induced_field, components=["X"])

    cable_active = Cable([tl_active], components=["X"])
    cable_passive = Cable([tl_passive], components=["X"])

    potentials_active, distances = cable_active._pot_along_cable_(timestamp=0)
    potentials_passive, _ = cable_passive._pot_along_cable_(timestamp=0)

    fig, ax = plt.subplots(figsize=(6, 3), dpi=300)
    ax.plot(distances, potentials_active, label="Active terminations", color="tab:red")
    ax.plot(distances, potentials_passive, label="Passive", color="tab:blue", linestyle="--")
    ax.set_xlabel("Distance along cable (km)")
    ax.set_ylabel("Potential (V)")
    ax.set_ylim([-100, 100])
    ax.set_xlim([0, length])
    ax.legend(loc="upper right")
    ax.set_title("Cable potential comparison")
    plt.show()
    ```
