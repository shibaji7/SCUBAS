<!-- 
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is under the MIT license found in the root directory LICENSE.md 
Everyone is permitted to copy and distribute verbatim copies of this license 
document.

This version of the MIT Public License incorporates the terms
and conditions of MIT General Public License.
-->

# Cable Modelling Utilities

The refactored cable module provides a clear separation between geometry
(`CableSection`), electrical properties (`TransmissionLine`), and the nodal
solution (`Cable`).  Nested configuration dictionaries are normalised using
`RecursiveNamespace`, and constructors now validate required arguments (for
example at least one component must be supplied when instantiating a `Cable`).

## Quick start

```python
from scubas.cables import Cable, TransmissionLine
from scubas.datasets import PROFILES

# Build two transmission-line sections backed by the default profiles
section_a = TransmissionLine(
    "CS-001",
    directed_length={"length": 2.5},
    elec_params={"site": PROFILES.CS},
)
section_b = TransmissionLine("CS-002", directed_length={"length": 1.8})

# Populate the ocean model and generate equivalent π-circuit parameters
for section in (section_a, section_b):
    section.compile_oml()
    section.compute_eqv_pi_circuit()

# Assemble the cable – the constructor runs compile() automatically
cable = Cable([section_a, section_b], components=["X", "Y"])
print(cable.tot_params[["V(v)", "E.X"]].head())
```

This snippet mirrors the new workflow: a section handles geometry and
conductivity, `TransmissionLine` converts B-fields to E-fields via
`OceanModel`, and `Cable` orchestrates nodal analysis and stores the derived
statistics in `tot_params` and `result`.

## API reference

### CableSection

::: scubas.cables.CableSection
    handler: python
    options:
      members:
        - check_location
        - compute_lengths
        - _pot_alongCS_
      show_root_heading: true
      show_source: true

### TransmissionLine

::: scubas.cables.TransmissionLine
    handler: python
    options:
      members:
        - to_str
        - compile_oml
        - add_active_termination
        - calc_trasmission_line_parameters
        - compute_eqv_pi_circuit
        - compute_Vj
      show_root_heading: true

### Cable

::: scubas.cables.Cable
    handler: python
    options:
      members:
        - compile
        - run_nodal_analysis
        - solve_admitance_matrix
        - consolidate_final_result
        - save
        - _pot_endCS_byComp_
        - _pot_endCS_
        - _pot_end_cable_byComp_
        - _pot_end_cable_
        - _pot_along_cable_
      show_root_heading: true
      show_source: true
