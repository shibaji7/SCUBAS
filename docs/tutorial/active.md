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


In this case, electromagnetic induction in the transmission line in 'Active Termination' figure can be represented by the circuit shown in following figure.

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