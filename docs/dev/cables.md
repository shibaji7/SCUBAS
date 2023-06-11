<!-- 
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is under the MIT license found in the root directory LICENSE.md 
Everyone is permitted to copy and distribute verbatim copies of this license 
document.

This version of the MIT Public License incorporates the terms
and conditions of MIT General Public License.
-->

::: scubas.cables.CableSection
    handler: python
    options:
      members:
        - check_location
        - compute_lengths
      show_root_heading: true
      show_source: true
      
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
        - _pot_alongCS_
      show_root_heading: true
      
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