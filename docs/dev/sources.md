# Current source modelling

<span class="api-badge api-package">Package</span>
<span class="api-badge api-module">Module</span> `scubas.sources.current_sources`
<span class="api-badge api-class">Class</span> `Direction`, `Current`

This module models ionospheric line-current sources and computes associated
magnetic and electric perturbations.

!!! important "Unit conventions"
    Default inputs use `km` for length and `MA` for current. Internally these
    values are converted to SI units (`m`, `A`) before field calculations.

## Quick start

```python
import numpy as np
from scubas.sources.current_sources import Current, Direction

x = np.linspace(-500, 500, 201)  # km
source = Current(
    h=110.0,
    I=1.0,
    current_direction=Direction.GEO_EAST_WEST,
)
source.compute_B_field(x=x, length_unit="km")

print(source.B_fields[0]["direction"], source.B_field_units)
print(source.E_fields[0]["direction"], source.E_field_units)
```

## API reference

::: scubas.sources.current_sources.Direction
    handler: python
    options:
      members:
        - enumerate_key
      show_root_heading: true
      show_source: true

::: scubas.sources.current_sources.Current
    handler: python
    options:
      members:
        - compute_B_field
      show_root_heading: true
      show_source: true
