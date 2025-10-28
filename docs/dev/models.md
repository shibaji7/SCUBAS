<!-- 
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is under the MIT license found in the root directory LICENSE.md 
Everyone is permitted to copy and distribute verbatim copies of this license 
document.

This version of the MIT Public License incorporates the terms
and conditions of MIT General Public License.
-->

# Ocean and electric field modelling

The `scubas.models` module orchestrates the conversion between magnetic and
electric fields using the modernised `OceanModel` class.  The refactor
introduced tighter validation (e.g. positive frequency limits), better error
handling when reading IAGA/CSV data, and a typed `Preprocess` helper that
encapsulates detrending and tapering.

Highlights:

* `OceanModel` now requires a `Site` instance and validates the frequency
  range (`flim`) at construction time.
* Transfer functions are cached after a call to `calcTF`, and the optional
  `compile_oml` flow used by transmission lines leverages these helpers.
* `read_Bfield_data` reads both IAGA text files and simple CSV exports,
  returning a single time-aligned dataframe ready for FFT conversion.
* `to_Efields` handles time stamps expressed either as datetimes or seconds,
  storing the results in `self.Efield` for downstream usage.

## Quick start

```python
import pandas as pd
import numpy as np

from scubas.datasets import PROFILES
from scubas.models import OceanModel

model = OceanModel(PROFILES.CS, flim=(1e-4, 1e-2))

# Create a synthetic B-field dataframe (two components, 1-second cadence)
times = pd.date_range("2024-01-01", periods=32, freq="s")
bfield = pd.DataFrame(
    {
        "X": np.sin(np.linspace(0, 2 * np.pi, len(times))),
        "Y": np.cos(np.linspace(0, 2 * np.pi, len(times))),
    },
    index=times,
)

# Convert to electric fields using the cached transfer functions
model.to_Efields(bfield)
print(model.Efield.head())
```

## API reference

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
