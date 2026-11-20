# Current source modelling

- **Module:** `scubas.sources.current_sources`
- **Class:** `Direction`, `Current`

This module models ionospheric line-current sources and computes associated
magnetic and electric perturbations.

:::{admonition} Unit conventions
:class: important

Default inputs use `km` for length and `MA` for current. Internally these
values are converted to SI units (`m`, `A`) before field calculations.
:::

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

```{eval-rst}
.. autoclass:: scubas.sources.current_sources.Direction
   :members: enumerate_key
```

```{eval-rst}
.. autoclass:: scubas.sources.current_sources.Current
   :members: compute_B_field
```
