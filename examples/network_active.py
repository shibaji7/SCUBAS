"""Example showcasing active-termination modelling on a continental shelf.

This script mirrors the narrative from the developer documentation.  It
demonstrates how the refactored API pieces come together:

* configure plotting using ``scubas.plotlib.update_rc_params``;
* build an ``OceanModel`` and fetch its transfer function;
* convert an induced electric field into equivalent transmission-line
  parameters; and
* evaluate the voltage induced along a short cable with active terminations.
"""

# NOTE: Matplotlib is used for plotting the transfer function and the cable
# potentials.  ``mpl`` is unused directly but kept for clarity if additional
# configuration is desired.
import matplotlib as mpl  # noqa: F401  (kept for parity with docs)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Core SCUBAS imports: cable modelling primitives, canonical conductivity
# profiles, ocean-model utilities, and plotting helpers.
from scubas.cables import Cable, TransmissionLine
from scubas.datasets import PROFILES
from scubas.models import OceanModel
from scubas.plotlib import plot_transfer_function, update_rc_params


def main() -> None:
    """Run the active termination demonstration step-by-step."""

    # Configure Matplotlib so the output mirrors the documentation styling.
    # Passing ``science=False`` keeps the default theme while still applying
    # the sans-serif fonts to make the plot consistent across environments.
    update_rc_params(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"],
        },
        science=True,
        latex=True,
    )

    # Instantiate the ocean model for the deep-ocean profile and visualise
    # the associated transfer function (how magnetic fields map to electric
    # fields in the layered medium).
    ocean_site = PROFILES.DO_3
    ocean_model = OceanModel(ocean_site)
    transfer_function = ocean_model.get_TFs()
    tf_artifacts = plot_transfer_function(transfer_function)
    tf_artifacts.figure.suptitle("Deep Ocean Transfer Function")

    ####################################################################
    # Simulated case:
    #   * Induced electric field of 0.3 V/km (X component)
    #   * Continental shelf 600 km in length, 100 m depth (via profile CS)
    ####################################################################

    # Continental shelf length (km) and induced electric field (0.3 V/km ->
    # 300 mV/km).  The DataFrame mirrors the structure expected by
    # ``compute_eqv_pi_circuit``: component columns plus an index representing
    # time samples.
    length_km = 600.0
    induced_e_field = pd.DataFrame(
        {"X": np.array([300.0])},
        index=pd.RangeIndex(1, name="Time"),
    )

    # Configure the transmission line, pointing both terminations to the land
    # model (PROFILES.LD).  The ``directed_length`` only needs the northward
    # component here, but other keys (e.g. ``length_east``) are also accepted.
    transmission_line_active = TransmissionLine(
        sec_id="CS-active",
        directed_length={"length_north": length_km},
        elec_params={"site": PROFILES.CS, "width": 1.0, "flim": [1e-6, 1.0]},
        active_termination={
            "right": {"site": PROFILES.LD, "width": 1.0},
            "left": {"site": PROFILES.LD, "width": 1.0},
        },
    )

    # For comparison, repeat the configuration without any active terminations.
    transmission_line_passive = TransmissionLine(
        sec_id="CS-passive",
        directed_length={"length_north": length_km},
        elec_params={"site": PROFILES.CS, "width": 1.0, "flim": [1e-6, 1.0]},
    )

    # Convert the synthetic electric field into π-circuit parameters (Ye, Yp2,
    # Ie).  These feed directly into the nodal analysis executed by ``Cable``.
    transmission_line_active.compute_eqv_pi_circuit(
        Efield=induced_e_field, components=["X"]
    )
    transmission_line_passive.compute_eqv_pi_circuit(
        Efield=induced_e_field, components=["X"]
    )

    # Assemble a single-section cable; the constructor triggers ``compile()``
    # and therefore runs the nodal analysis immediately.
    # Assemble both cables (one with active termination, one without) and
    # evaluate the potential distributions side-by-side.
    cable_active = Cable([transmission_line_active], components=["X"])
    cable_passive = Cable([transmission_line_passive], components=["X"])

    # Sample the potential distribution along the cable for timestamp zero and
    # render the results.  The helper returns a dictionary with the figure and
    # axes so users can further customise the plot.
    potentials_active, distances = cable_active._pot_along_cable_(timestamp=0)
    potentials_passive, _ = cable_passive._pot_along_cable_(timestamp=0)

    # Produce an overlaid plot comparing active vs passive configurations.
    fig, ax = plt.subplots(figsize=(3, 3), dpi=180)
    ax.plot(distances, potentials_active, label="Active terminations", color="tab:red")
    ax.plot(
        distances, potentials_passive, label="Passive", color="tab:blue", linestyle="--"
    )
    ax.set_xlabel("Distance along cable (km)")
    ax.set_ylabel("Potential (V)")
    ax.set_ylim([-100, 100])
    ax.set_xlim([0, length_km])
    ax.legend(loc="upper right")
    ax.set_title("Cable potential comparison")

    plt.show()


if __name__ == "__main__":
    main()
