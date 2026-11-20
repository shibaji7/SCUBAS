<!--
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is licensed under the GNU Lesser General Public License v3.0 or later
(LGPL-3.0-or-later); see LICENSE and COPYING in the root directory. Versions up to
and including 0.1.9 were released under the MIT license (see NOTICE.md).
-->

# Utility helpers

- **Module:** `scubas.utils`
- **Class:** `RecursiveNamespace`, `GreatCircle`
- **Method / Function:** `fft`, `ifft`, `frexp102str`, `component_mappings`, `component_sign_mappings`

:::{admonition} Core dependency module
:class: important

Utilities here are reused across cable, conductivity, and model workflows.
Changes in this module can impact multiple parts of SCUBAS.
:::

`scubas.utils` provides the small supporting functions and classes used across
the codebase.  The refactor introduced dataclass-friendly namespaces, more
robust FFT scaling, and guards against invalid input values.

Highlights:

* `RecursiveNamespace` recursively converts mapping/list entries, which is
  used by the cable/conductivity modules to normalise configuration dictionaries.
* `frexp102str` now accepts zero and negative values, returning a formatted
  scientific-notation string.
* `fft`/`ifft` wrap NumPy’s real FFT operations with consistent scaling and
  input validation.
* `GreatCircle` exposes both great-circle and haversine distance utilities with
  input checking.

## Quick start

```python
import numpy as np
from scubas.utils import RecursiveNamespace, fft, ifft, GreatCircle

cfg = RecursiveNamespace(alpha=1, nested={"beta": 2})
print(cfg.nested.beta)

t = np.linspace(0, 1, 128, endpoint=False)
signal = np.sin(2 * np.pi * 5 * t)
spectrum, freqs = fft(signal, dT=t[1] - t[0])
reconstructed = ifft(spectrum)

start = RecursiveNamespace(lat=0.0, lon=0.0)
end = RecursiveNamespace(lat=0.0, lon=90.0)
gc = GreatCircle(start, end)
print(gc.haversine())
```

## API reference

```{eval-rst}
.. autoclass:: scubas.utils.RecursiveNamespace
   :members:
```

```{eval-rst}
.. autofunction:: scubas.utils.frexp102str
```

```{eval-rst}
.. autofunction:: scubas.utils.fft
```

```{eval-rst}
.. autofunction:: scubas.utils.ifft
```

```{eval-rst}
.. autofunction:: scubas.utils.component_mappings
```

```{eval-rst}
.. autofunction:: scubas.utils.component_sign_mappings
```

```{eval-rst}
.. autoclass:: scubas.utils.GreatCircle
   :members: great_circle, haversine
```
