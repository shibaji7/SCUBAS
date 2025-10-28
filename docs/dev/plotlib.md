<!-- 
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is under the MIT license found in the root directory LICENSE.md 
Everyone is permitted to copy and distribute verbatim copies of this license 
document.

This version of the MIT Public License incorporates the terms
and conditions of MIT General Public License.
-->

# Plotting helpers

`scubas.plotlib` now offers lightweight plotting utilities that return a
typed `PlotArtifacts` dataclass containing the `Figure` and `Axes` objects.
The helpers respect optional “science” styling, accept numpy arrays directly,
and summarise the most common inspection plots used in the cable/conductivity
workflows.

## Quick start

```python
import numpy as np
from scubas.plotlib import PlotArtifacts, plot_transfer_function

class DummyTF:
    freq = np.logspace(-4, -2, 5)
    E2B = np.exp(1j * np.linspace(0, np.pi / 2, 5))

artifacts = plot_transfer_function(DummyTF())
assert isinstance(artifacts, PlotArtifacts)
artifacts.figure.savefig("transfer_function.png")
```

Similar helpers are available for plotting along-section potentials and whole
Cable potentials.

## API reference

::: scubas.plotlib.update_rc_params
    options:
      show_root_heading: true
      show_source: true

::: scubas.plotlib.plot_transfer_function
    options:
      show_root_heading: true
      show_source: true

::: scubas.plotlib.potential_along_section
    options:
      show_root_heading: true
      show_source: true

::: scubas.plotlib.cable_potential
    options:
      show_root_heading: true
      show_source: true
