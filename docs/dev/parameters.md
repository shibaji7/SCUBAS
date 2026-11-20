# Parameters and provenance

Every number a model depends on should say where it came from. SCUBAS 0.1.9 predates that rule: its constants are written into the code. So there are two parts, and they do
not mix.

## For new code: `scubas.provenance`

A parameter is a table `{value, unit, status, source}` in a TOML file, loaded as a frozen `Param`:

```toml
[line]
length_km = { value = 750.8, unit = "km", status = "assumed", source = "the example link of the tutorial" }
tempco    = { value = "TBD", unit = "1/K", status = "unknown", source = "needed for the temperature model; not published" }
```

| Status | Meaning |
|---|---|
| `sourced` | taken from a document; the source names it |
| `assumed` | a placeholder, or a legacy default whose origin is not documented |
| `unknown` | the value is `"TBD"`; asking for it raises `MissingParameterError` with the reason, instead of a guess |
| `override` | changed for one run (a bare value in a user file, or an explicit override) |

```python
from scubas.provenance import load_parameters

ps = load_parameters("defaults.toml", path="mine.toml", overrides={"line.length_km": 600.0})
ps.get("line", "length_km")      # 600.0, status "override"
ps.get("line", "tempco")         # raises MissingParameterError: line.tempco is unknown (TBD): needed for ...
ps.counts(), ps.fingerprint()    # how many of each status, and a short hash to trace a result to its inputs
```

A file shipped with the package must give all four fields. A user file may give a bare value, and a key that is not in the defaults is an error (a typo, not a new setting). The `unit` is a
label in this first version; nothing is converted. No existing module reads these tables: the rule applies to **new** code.

## For the existing code: the list of legacy constants

`scubas/defaults/legacy_constants.toml` lists the constants that the existing modules use (`legacy_constants()` loads it), each with its status. It is **read-only**: the code does not read
it. It is a record, and `tests/test_provenance.py` fails if a value differs from the value in the code, so the list cannot drift. A status is honest: where the code gives no citation the
status is `assumed` and the source says so; no reference is invented. Where a number can be traced to a document the status is `sourced` and names it.

To add a constant: add the entry with its unit, status and source, add its check to `CHECKS` in `tests/test_provenance.py` (the test fails until you do), and run
`python scripts/gen_param_table.py` to refresh the table below. Moving a constant out of the code, so that the code reads it from a file, is a change of the legacy modules and
needs its own approved pull request with the compatibility suite unchanged ([ADR 0001](adr/0001-legacy-numerics-are-pinned.md)).

## The legacy constants

<!-- PARAMS:START -->
| key | value | unit | status | source |
|---|---|---|---|---|
| `conductivity.seawater_resistivity_ohm_m` | 0.3 | ohm m | sourced | default of ConductivityProfile; the layer table of docs/tutorial/netmodel.md, from the LITHO1.0-based study of Chakraborty et al. (2022), Front. Phys. 10:1022475 |
| `conductivity.sediment_resistivity_ohm_m` | 3.0 | ohm m | sourced | default of ConductivityProfile; the layer table of docs/tutorial/netmodel.md (Chakraborty et al. 2022) |
| `conductivity.crust_resistivity_ohm_m` | 3000.0 | ohm m | sourced | default of ConductivityProfile; the layer table of docs/tutorial/netmodel.md (Chakraborty et al. 2022) |
| `conductivity.lithosphere_resistivity_ohm_m` | 1000.0 | ohm m | sourced | default of ConductivityProfile; the layer table of docs/tutorial/netmodel.md (Chakraborty et al. 2022) |
| `conductivity.asthenosphere_resistivity_ohm_m` | 100.0 | ohm m | sourced | default of ConductivityProfile (the 'Upper Mantle' layer); the layer table of docs/tutorial/netmodel.md (Chakraborty et al. 2022) |
| `conductivity.transition_zone_resistivity_ohm_m` | 10.0 | ohm m | sourced | default of ConductivityProfile; the layer table of docs/tutorial/netmodel.md (Chakraborty et al. 2022) |
| `conductivity.lower_mantle_resistivity_ohm_m` | 1.0 | ohm m | sourced | default of ConductivityProfile; the layer table of docs/tutorial/netmodel.md (Chakraborty et al. 2022) |
| `conductivity.transition_zone_top_km` | 410.0 | km | assumed | default of ConductivityProfile: the depth at which the transition zone starts (the conventional 410 km seismic discontinuity); no citation in the code |
| `conductivity.transition_zone_bottom_km` | 660.0 | km | assumed | default of ConductivityProfile: the depth at which the transition zone ends (the conventional 660 km seismic discontinuity); no citation in the code |
| `conductivity.profile_max_depth_km` | 1000.0 | km | assumed | default of ConductivityProfile: the profile ends at this depth, so the lower mantle is 1000 - 660 km thick; no citation in the code |
| `conductivity.grid_interpolation_method` | nearest |  | assumed | default of ConductivityProfile: how the LITHO1.0 grid is interpolated ('nearest' or 'linear') |
| `conductivity.earth_model_file` | LITHO1.0.nc |  | sourced | default of ConductivityProfile: the file name of the LITHO1.0 netCDF product distributed by the IRIS Earth Model Collaboration; cached in .scubas_config/ of the working directory |
| `conductivity.earth_model_uri` | http://ds.iris.edu/files/products/emc/emc-files/LITHO1.0.nc |  | sourced | default of ConductivityProfile: the IRIS Earth Model Collaboration copy of LITHO1.0; plain HTTP and no integrity check (see SECURITY.md) |
| `preprocess.taper_fraction` | 0.1 |  | assumed | default p of Preprocess: the fraction of the record that the cosine taper covers before the transform to the frequency domain |
| `preprocess.baseline_samples` | 120 | samples | assumed | hard-coded in Preprocess.detrend_magnetic_field: the median of the first 120 samples is removed as the baseline |
| `models.flim_min_Hz` | 1e-06 | Hz | assumed | default of OceanModel: the lowest frequency of the transfer-function grid, which is also its spacing (round(fmax/fmin) + 1 points) |
| `models.flim_max_Hz` | 1.0 | Hz | assumed | default of OceanModel: the highest frequency of the transfer-function grid |
| `cables.width_m` | 1.0 | m | assumed | default width of TransmissionLine (the docstring says metres): the cable's cross-section that scales the shunt admittance |
| `cables.flim_min_Hz` | 1e-06 | Hz | assumed | default of TransmissionLine elec_params, passed to OceanModel |
| `cables.flim_max_Hz` | 1.0 | Hz | assumed | default of TransmissionLine elec_params, passed to OceanModel |
| `geometry.north_east_split_divisor` | 1.4142135623730951 |  | assumed | CableSection: a total length given without a direction is split equally, length_north = length_east = length / sqrt(2) |
| `geometry.lat_km_per_degree_a` | 111.133 | km/degree | assumed | CableSection.compute_lengths: length of a degree of latitude is (a - b cos 2 lambda) with lambda the mean latitude; the origin of a and b is not documented in the code |
| `geometry.lat_km_per_degree_b` | 0.56 | km/degree | assumed | CableSection.compute_lengths: see lat_km_per_degree_a |
| `geometry.lon_km_per_degree_a` | 111.5065 | km/degree | assumed | CableSection.compute_lengths: length of a degree of longitude is (a - b cos 2 lambda) cos lambda; the origin of a and b is not documented in the code |
| `geometry.lon_km_per_degree_b` | 0.1872 | km/degree | assumed | CableSection.compute_lengths: see lon_km_per_degree_a |
| `utils.earth_radius_km` | 6371.0 | km | assumed | default Re of GreatCircle: the mean radius of a spherical Earth, rounded |
| `sources.current_height_km` | 110.0 | km | assumed | default h of Current: the height of the ionospheric line current (given in km by default) |
| `sources.cauchy_width` | 0.0 | km | assumed | default a of Current: the width of the Cauchy-distributed current (km by default) |
| `sources.current_amplitude` | 1.0 | MA | assumed | default I of Current: the current (mega-amperes by default) |
| `readers.iaga_missing_values` | [99999.9, 99999.0, 88888.8, 88888.0] |  | sourced | OceanModel.read_iaga: the IAGA-2002 data exchange format marks a missing value with 99999 and an element that is not recorded with 88888; these become gaps |
| `readers.gap_spline_order` | 2 |  | assumed | OceanModel.read_Bfield_data: gaps in a component are filled by Series.interpolate(method='spline', order=2) |
| `monte_carlo.default_draws` | 100 |  | assumed | default n of ConductivityProfile.compile_mcmc_bined_profiles |
| `monte_carlo.default_random_seed` | 0 |  | assumed | default random_seed of ConductivityProfile.compile_mcmc_bined_profiles (numpy's global seed is set to it) |
| `monte_carlo.shelf_thickness_min_km` | 0.01 | km | assumed | default lower bound of the seawater thickness drawn for the end bins in compile_mcmc_bined_profiles |
| `monte_carlo.shelf_thickness_max_km` | 1.0 | km | assumed | default upper bound of the seawater thickness drawn for the end bins in compile_mcmc_bined_profiles |
<!-- PARAMS:END -->
