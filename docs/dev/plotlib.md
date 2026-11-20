<!--
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is licensed under the GNU Lesser General Public License v3.0 or later
(LGPL-3.0-or-later); see LICENSE and COPYING in the root directory. Versions up to
and including 0.1.9 were released under the MIT license (see NOTICE.md).
-->

# Plotting helpers

- **Module:** `scubas.plotlib`
- **Class:** `PlotArtifacts`
- **Method / Function:** `update_rc_params`, `plot_transfer_function`, `potential_along_section`, `cable_potential`

:::{admonition} Return type
:class: note

Plot helpers return a `PlotArtifacts` object to make figure/axes reuse and
custom post-processing easier.
:::

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

```{eval-rst}
.. autoclass:: scubas.plotlib.PlotArtifacts
   :members:
```

```{eval-rst}
.. autofunction:: scubas.plotlib.update_rc_params
```

```{eval-rst}
.. autofunction:: scubas.plotlib.plot_transfer_function
```

```{eval-rst}
.. autofunction:: scubas.plotlib.potential_along_section
```

```{eval-rst}
.. autofunction:: scubas.plotlib.cable_potential
```
