# Glossary

Terms and symbols used in these pages and in the code, with where each one comes from. Equations and sources are in the
[theory](tutorial/theory.md) and [electrically long and short cables](tutorial/elsc.md) tutorials; the code names are those of
`scubas.cables`, `scubas.datasets` and `scubas.models`.

## Symbols

| Symbol | Meaning | Unit | In the code |
|---|---|---|---|
| $E$ | geoelectric field along a cable section | mV/km in the data frames, V/m inside the section equations (the code multiplies by $10^{-6}$) | `Efield` |
| $B$ | geomagnetic field of the observatory record | nT | `Bfield` |
| $Z$ | series impedance of a section per unit length | $\Omega$/m (the printed summary of `TransmissionLine.to_str` is per km) | `TransmissionLine.Z` |
| $Y$ | shunt admittance of a section per unit length | S/m (summary per km) | `TransmissionLine.Y` |
| $\gamma=\sqrt{ZY}$ | propagation constant | 1/m | `gma` |
| $1/\gamma$ | [adjustment distance](#adjustment-distance) | m (the summary prints km as `Ad`) | |
| $Z_0=\sqrt{Z/Y}$ | characteristic impedance | $\Omega$ | `Z0` |
| $L$ | length of a section | km in the inputs, m in the equations | `length`, `length_north`, `length_east` |
| $V_i$, $V_k$ | potentials at the two ends of a section | V | `end_pot.Vi`, `end_pot.Vk` |
| $Y_E$, $Y'/2$, $I_E$ | [equivalent-pi](#equivalent-pi-circuit) branch admittance, admittance to ground at each end, and equivalent current source | S, S, A | `Ye`, `Yp2`, `Ie` |
| $U_0$, $U_1$ | potentials at the two ends of the whole cable | V | columns `U0`, `U1` of `Cable.tot_params` |
| `V(v)` | sum of the section potentials | V | column `V(v)` |
| `Vt(v)` | `U0 - U1 + V(v)` as computed from end potentials rounded to 0.01 V (see [known defects](dev/testing.md#known-defects)) | V | column `Vt(v)` |

## Terms

### Active termination
Instead of grounding the end of a section directly, the structure beyond it (for example land) is represented by an equivalent
circuit that carries its own induced field. The code builds it from a `Site` and the same B-field files, and stores its
admittance `Yn = 1/Z0` and current source `Jn = E/Z` in `TransmissionLine.term_params`. See [active
terminations](tutorial/active.md).

### Adjustment distance
$1/\gamma$: the distance over which an end effect decays in a section. It decides whether a section is electrically long or short.

### Alpha
The development status of SCUBAS (version 0.x, classifier "Development Status :: 3 - Alpha"): published method, results compared
with observations in the cited studies, but an API and numbers that may still change.

### B2B and E2B
The two transfer functions of `OceanModel` (`calcTF`, `get_TFs`): `E2B` gives the electric field per unit magnetic field
(mV/km per nT) and `B2B` the magnetic-field-to-magnetic-field transfer of the ocean layer; the code writes them as
`Zd/(cosh kd + (Zd/Z) sinh kd)` and `1/(cosh kd + (Zd/Z) sinh kd)`.

### Compatibility suite
`tests/compat/`: the tests that pin the public API, the file formats and the numbers of SCUBAS 0.1.9, so that refactoring and
tooling work cannot change what users get. See [testing](dev/testing.md).

### Distributed-source transmission line (DSTL)
The model of the ocean and the Earth as a transmission line with a series impedance $Z$ from the conducting layers, a parallel admittance $Y$ from the resistive layers, and the induced electric field as a distributed source (Boteler et al., 2023). It gives the Earth potentials at the ends of a cable. See [the DSTL model](model/dstl.md).

### Electrically long and electrically short
By the thresholds of the [tutorial](tutorial/elsc.md): a section is electrically long when its length exceeds four times the
adjustment distance ($L>4/\gamma$), where the end potential tends to $-E/\gamma$; and electrically short when it is shorter than
the adjustment distance ($L<1/\gamma$), where the end potential tends to $-LE/2$. Both limits are checked against the code in
`tests/compat/test_analytic.py`.

### Equivalent-pi circuit
The reduction of a transmission-line section to a series branch and an admittance to ground at each end, driven by a voltage
source $E'$ or, equivalently, by a current source $I_E$ (Boteler, 2013). It is what lets sections with different conductivity
profiles be joined in one nodal network.

### Extreme value analysis
The statistical treatment of the 1989 TAT-8 voltages in the [EVA tutorial](tutorial/eva.md): return periods of the largest
voltages at the cable ends and along the span.

### IAGA-2002
The text format of geomagnetic observatory records (`.txt` and `.min` files). `OceanModel.read_iaga` reads it; missing values
(`99999.00` and similar) become gaps that `read_Bfield_data` fills by second-order spline interpolation.

### INTERMAGNET
The global network of geomagnetic observatories whose records drive SCUBAS; see [INTERMAGNET data access](user/intermagnet.md).

### Provenance
The unit, status (`sourced`, `assumed`, `unknown`, `override`) and source that a parameter carries, so that a number is never separated from where it came from. The pattern is `scubas.provenance`
for new code; the constants of the existing code are listed with their provenance in `legacy_constants.toml`. See [parameters and provenance](dev/parameters.md).

### Future slot
A planned capability kept as a stub in `scubas.future` that raises `NotImplementedYetError` naming its roadmap item. The package is empty until a capability
is named. See [planned capabilities](dev/future.md).

### Known defect
Behaviour of 0.1.9 that looks wrong but that users or published results may depend on. It is pinned by a test as it is and
reported; a fix is a separate, approved change. See [testing](dev/testing.md#known-defects).

### LITHO1.0
The global model of the lithosphere (a netCDF file) from which `ConductivityProfile` builds layered profiles of seawater,
sediments, crust, lithosphere and mantle. It is downloaded on first use into `.scubas_config/` of the working directory.

### Nodal analysis
The solution of the equivalent-pi network of a multi-section cable: a nodal admittance matrix and the injected currents give the
potential at every node (`Cable.run_nodal_analysis`, `Cable.solve_admitance_matrix`). See [network modelling](tutorial/netmodel.md).

### Norton equivalent
A current source in parallel with an admittance. The land at each end of a cable is represented by its Norton equivalent in the
[network model](tutorial/netmodel.md).

### Profile and site
A layered Earth model: a `Site` is a list of `Layer` objects (name, thickness in metres, conductivity in S/m), and
`scubas.datasets.PROFILES` holds the shipped ones.

| Key | Profile |
|---|---|
| `CS`, `SO`, `DO` | continental shelf, shallow ocean and deep ocean reference profiles used in the IGS analyses |
| `LD0`, `LD` | land profiles from the IGS study (`LD` is the one used for terminations) |
| `CS_W`, `CS_E`, `DO_1` to `DO_6`, `MAR` | the sections of the March 1989 TAT-8 study, including the Mid-Atlantic Ridge |
| `BM`, `OM`, `DB`, `Quebec`, `UN` | benchmark and study profiles (a benchmark model, the same with a 1 km ocean layer, a Quebec model, and a uniform model) |

### Skin depth
The depth $p$ at which an electromagnetic field has decayed, computed from the impedance of the layered Earth:
$p=Z/(i\omega\mu_0)$ (`Site.calcP`).

### Thin sheet and transmission-line model (TLM)
A conductive layer over a resistive one is represented as a transmission line with series impedance $Z$ and shunt admittance $Y$
(generalised thin-sheet analysis, Ranganayaki and Madden, 1980; as a transmission line, Wang et al., 2023). See the
[theory](tutorial/theory.md) page.

### TAT-8
The trans-Atlantic telecommunications cable whose March 1989 voltages are reproduced in the [1989 event study](tutorial/1989.md).
Its observatory records come from the FRD, STJ and HAD stations.
