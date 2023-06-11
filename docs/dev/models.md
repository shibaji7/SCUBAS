<!-- 
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is under the MIT license found in the root directory LICENSE.md 
Everyone is permitted to copy and distribute verbatim copies of this license 
document.

This version of the MIT Public License incorporates the terms
and conditions of MIT General Public License.
-->

::: scubas.models.OceanModel
    handler: python
    options:
      members:
        - calcZ
        - calcTF
        - get_TFs
        - read_iaga
        - read_Bfield_data
        - to_Efields
      show_root_heading: true
      show_source: true
      
::: scubas.models.Preprocess
    handler: python
    options:
      members:
        - get_tapering_function
        - detrend_magnetic_field
      show_root_heading: true
      show_source: true