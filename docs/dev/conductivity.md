<!-- 
Author(s): Shibaji Chakraborty, Xueling Shi

Disclaimer:
SCUBAS is under the MIT license found in the root directory LICENSE.md 
Everyone is permitted to copy and distribute verbatim copies of this license 
document.

This version of the MIT Public License incorporates the terms
and conditions of MIT General Public License.
-->

# Conductivity profiles

The conductivity utilities were modernised to provide type-safe configuration,
download-on-demand handling of the LITHO1.0 model, and a consistent interface
that can emit either raw `pandas.DataFrame` outputs or ready-to-use
`scubas.datasets.Site` objects.  All helper methods now accept numpy arrays
transparently and raise informative exceptions when model data is missing.

Highlights:

* `ConductivityProfile` merges user-provided dictionaries with sensible
  defaults (grid interpolation method, resistivities, NetCDF URI, etc.).
* The NetCDF file is cached in `.scubas_config/`; download and IO failures are
  surfaced as `RuntimeError`, making automation friendlier.
* Batch helpers (`compile_profiles`, `compile_bined_profiles`,
  `compile_mcmc_bined_profiles`) now operate on numpy/native lists and can
  return `Site` instances via `to_site=True`.

## Quick start

```python
from scubas.conductivity import ConductivityProfile

profile = ConductivityProfile()

# Build a site object for a single latitude/longitude pair
site = profile.compile_profile([45.0, -110.0])
print(site.get_names())

# Generate dataframe outputs for a list of coordinates
locations = [[40.0, -120.0], [42.5, -118.0]]
df_profiles = profile.compile_profiles(locations, to_site=False)
for df in df_profiles:
    print(df[["name", "thickness", "resistivity"]])

# Interpolate along a binned track and keep Site output
bined = np.array([[40.0, -120.0], [41.0, -119.0], [42.0, -118.5]])
site_bins = profile.compile_bined_profiles(bined, to_site=True)
```

## API reference

::: scubas.conductivity.ConductivityProfile
    handler: python
    options:
      members:
        - load_earth_model
        - get_interpolation_points
        - get_water_layer
        - get_sediment_layer
        - get_crust_layer
        - get_lithosphere_layer
        - get_upper_mantle_layer
        - get_transition_zone_layer
        - get_lower_mantle_layer
        - _compile_profile_
        - compile_profile
        - compile_profiles
        - compile_bined_profile
        - compile_bined_profiles
        - compile_mcmc_bined_profiles
      show_root_heading: true
      show_source: true
