<!--
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is licensed under the GNU Lesser General Public License v3.0 or later
(LGPL-3.0-or-later); see LICENSE and COPYING in the root directory. Versions up to
and including 0.1.9 were released under the MIT license (see NOTICE.md).
-->

# Datasets module

- **Module:** `scubas.datasets`
- **Class:** `Layer`, `Site`, `PROFILES`

:::{admonition} Why this module matters
:class: note

Most SCUBAS workflows pass `Site` objects into conductivity, model, and
cable utilities. Defining these correctly improves downstream stability.
:::

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

```{eval-rst}
.. autoclass:: scubas.datasets.Layer
   :members:
```

```{eval-rst}
.. autoclass:: scubas.datasets.Site
   :members: get_thicknesses, get_conductivities, get_resistivities, get_names, get, set_thicknesses, calcZ, calcP, init
```

```{py:data} scubas.datasets.PROFILES
:type: types.SimpleNamespace

Ready-made `Site` objects, one attribute per profile (`PROFILES.CS`, `PROFILES.DO`, `PROFILES.LD`, ...). The keys and what each profile is are listed in the [glossary](../glossary.md#profile-and-site).
```
