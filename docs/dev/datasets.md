<!-- 
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is under the MIT license found in the root directory LICENSE.md 
Everyone is permitted to copy and distribute verbatim copies of this license 
document.

This version of the MIT Public License incorporates the terms
and conditions of MIT General Public License.
-->

::: scubas.datasets.Layer
    handler: python
    options:
      show_root_heading: true
      show_source: true
      
::: scubas.datasets.Site
    handler: python
    options:
      members:
        - get_thicknesses
        - get_conductivities
        - get_resistivities
        - get_names
        - get
        - calcZ
        - init
      show_root_heading: true
      show_source: true
      
::: scubas.datasets.PROFILES
    handler: python
    options:
      show_root_heading: true
      show_source: true