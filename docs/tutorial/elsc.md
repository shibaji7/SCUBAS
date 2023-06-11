<!-- 
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is under the MIT license found in the root directory LICENSE.md 
Everyone is permitted to copy and distribute verbatim copies of this license 
document.

This version of the MIT Public License incorporates the terms
and conditions of MIT General Public License.
-->
#### Electrically Long Cable
---
Consider an Ocean-Earth section with physical length, $L$, and propagation constant, $\gamma$. The section has an adjustment distance $\frac{1}{\gamma}$. For an electrically-long transmission line, where $L>\frac{4}{\gamma}$ (i.e. length is greater than four times of adjustment distance) we have the following scenario,

$$
e^{\gamma L}>>1>>e^{-\gamma L}
$$

When the transmission line length is considerably shorter than the adjustment distance it is referred to as 'electrically-long' and the equivalent-$\pi$ components reduce to (Boteler et a, 2013):

$$
Y_E=\frac{2}{Z_0e^{\gamma L}}
$$

$$
\frac{Y'}{2}=\frac{1}{Z_0}
$$

$$
I_E=\frac{E}{Z}
$$

$$
V_i=-U_k=-\frac{E}{\gamma}
$$

$$
V(x)=Ve^{-\gamma (L-x)}-Ve^{-\gamma x}
$$

where $V_i=-V_k=V=\frac{E}{\gamma}$.

!!! Example    
    ``` py
    # Import all required libs
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    # Set matplotlib styles
    plt.style.use(["science", "ieee"])
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans",
                                       "Lucida Grande", "Verdana"]
    import pandas as pd

    # Import SCUBAS dependencies
    from scubas.datasets import PROFILES
    from scubas.models import OceanModel
    from scubas.plotlib import plot_transfer_function,\
            potential_along_section, cable_potential, update_rc_params
    from scubas.cables import TransmissionLine, Cable
    from scubas.conductivity import ConductivityProfile as CP

    # Ocean-Earth conductivity profile of a Deep Ocean
    site = PROFILES.DO_3
    # Rended ocean model
    om = OceanModel(site)
    # Generate transfer function
    tf = om.get_TFs()
    # Transfer function of a Deep Ocean.
    # We are going to use this tranfer function to compute differen electrical cases
    _ = plot_transfer_function(tf)

    ####################################################################
    # Simulating the case: Induced electric field 0.3 V/km on a 
    # shallow continental shelf with depth 100 m, length 600 km
    ####################################################################
    e_CS = pd.DataFrame()
    e_CS["X"], e_CS["dTime"] = [300], [0]
    length=600
    tl = TransmissionLine(
        sec_id="CS",
        directed_length=dict(
            length_north=length,
        ),
        elec_params=dict(
            site=PROFILES.CS,
            width=1.0,
            flim=[1e-6, 1e0],
        ),
        active_termination=dict(
            right=None,
            left=None,
        ),
    )
    tl.compute_eqv_pi_circuit(e_CS, ["X"])
    cable = Cable([tl], ["X"])
    Vc, Lc = cable._pot_along_cable_(0)
    tag = cable_potential(Vc, Lc, ylim=[-200, 200])
    _ = tag["axes"].text(0.05, 0.85, r"$L_{cs}$=%dkm"%length, ha="left",\
                    va="center", transform=tag["axes"].transAxes)
    ```

#### Electrically Short Cable
For an electrically-short section, where the physical length is less than the adjustment distance, i.e., $L<\frac{1}{\gamma}$. When the transmission line length is considerably shorter than the adjustment distance it is referred to as 'electrically-short' and the equivalent-$\pi$ components reduce to (Boteler et a, 2013):

$$
e^{\pm\gamma L}\approx 1\pm\gamma L
$$

$$
Y_E=\frac{1}{Z_0\gamma L}
$$

$$
\frac{Y'}{2}=0
$$

$$
I_E=\frac{E}{Z}
$$

$$
V_i=-V_k=-\frac{LE}{2}=V
$$

$$
V(x)=\frac{V}{2L}(2+\gamma L)(2x-L)
$$

!!! Example
    ``` py
    ####################################################################
    # Simulating the case: Induced electric field 0.3 V/km on a 
    # shallow continental shelf with depth 100 m, length 4000 km
    ####################################################################
    length=4000
    tl = TransmissionLine(
        sec_id="CS",
        directed_length=dict(
            length_north=length,
        ),
        elec_params=dict(
            site=PROFILES.CS,
            width=1.0,
            flim=[1e-6, 1e0],
        ),
        active_termination=dict(
            right=None,
            left=None,
        ),
    )
    tl.compute_eqv_pi_circuit(e_CS, ["X"])
    cable = Cable([tl], ["X"])
    Vc, Lc = cable._pot_along_cable_(0)
    tag = cable_potential(Vc, Lc, ylim=[-200, 200])
    _ = tag["axes"].text(0.05, 0.85, r"$L_{cs}$=%dkm"%length, ha="left",\
            va="center", transform=tag["axes"].transAxes)
    ```