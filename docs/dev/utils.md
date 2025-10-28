<!-- 
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is under the MIT license found in the root directory LICENSE.md 
Everyone is permitted to copy and distribute verbatim copies of this license 
document.

This version of the MIT Public License incorporates the terms
and conditions of MIT General Public License.
-->

# Utility helpers

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

::: scubas.utils.RecursiveNamespace
    handler: python
    options:
      show_root_heading: true
      show_source: true

::: scubas.utils.frexp102str
    options:
      show_root_heading: true
      show_source: true

::: scubas.utils.fft
    options:
      show_root_heading: true
      show_source: true

::: scubas.utils.ifft
    options:
      show_root_heading: true
      show_source: true

::: scubas.utils.component_mappings
    options:
      show_root_heading: true
      show_source: true

::: scubas.utils.component_sign_mappings
    options:
      show_root_heading: true
      show_source: true

::: scubas.utils.GreatCircle
    handler: python
    options:
      members:
        - great_circle
        - haversine
      show_root_heading: true
      show_source: true
