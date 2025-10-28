<!-- 
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is under the MIT license found in the root directory LICENSE.md 
Everyone is permitted to copy and distribute verbatim copies of this license 
document.

This version of the MIT Public License incorporates the terms
and conditions of MIT General Public License.
-->

# Datasets module

`scubas.datasets` now exposes lightweight data classes powered by `dataclass`
and offers convenience helpers for synthesising layered conductivity profiles.
The modernised layout makes it easier to construct and inspect a `Site`, while
remaining compatible with existing utilities such as the conductivity and
cable modules.

Highlights:

* `Layer` is an immutable dataclass that computes resistivity lazily via a
  property.
* `Site` accepts any iterable of `Layer` objects and provides helpers for
  thickness/conductivity/resistivity extraction as lists or individual values.
* The `Site.init` static method builds a `Site` from parallel sequences of
  conductivities, thicknesses, and names—a pattern used widely in
  `PROFILES`.

## Quick start

```python
from scubas.datasets import Layer, Site, PROFILES

layers = [
    Layer(name="Seawater", thickness=1.0, conductivity=3.33),
    Layer(name="Sediments", thickness=2.5, conductivity=0.2),
]
site = Site(layers=layers, description="Example", name="Demo Site")

print(site.get_resistivities())

# Build from sequences
alt_site = Site.init(
    conductivities=[3.33, 0.2, 0.0003],
    thicknesses=[1.0, 2.5, float("inf")],
    names=["Seawater", "Sediments", "Mantle"],
    description="Three-layer example",
    site_name="Alt Site",
)

# Reuse predefined profiles
print(PROFILES.CS.get_names())
```

## API reference

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
        - set_thicknesses
        - calcZ
        - calcP
        - init
      show_root_heading: true
      show_source: true

::: scubas.datasets.PROFILES
    options:
      show_root_heading: true
      show_source: true
