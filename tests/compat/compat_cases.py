# SPDX-License-Identifier: MIT
# Copyright (c) 2023 Shibaji Chakraborty
"""Deterministic characterization cases: they run the legacy public API and return plain arrays and JSON-able metadata.

``scripts/gen_compat_baseline.py`` runs :func:`build_cases` against the unchanged 0.1.9 code and writes
``baseline/behaviour.npz`` and ``baseline/behaviour_meta.json``; ``test_behaviour.py`` repeats the run and compares.
Nothing here uses the network or random numbers (the Monte Carlo case fixes its seed through the public argument).
"""

from __future__ import annotations

import contextlib
import os
import warnings
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
from loguru import logger

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = HERE / "data"
DATASET_1989 = ROOT / "examples" / "datasets" / "1989"
SEP = "|"

LITHO_FIELDS = [
    "asthenospheric_mantle_top_depth",
    "lid_bottom_depth",
    "lid_top_depth",
    "lower_crust_bottom_depth",
    "lower_crust_top_depth",
    "middle_crust_bottom_depth",
    "middle_crust_top_depth",
    "upper_crust_bottom_depth",
    "upper_crust_top_depth",
    "lower_sediments_bottom_depth",
    "lower_sediments_top_depth",
    "middle_sediments_bottom_depth",
    "middle_sediments_top_depth",
    "upper_sediments_bottom_depth",
    "upper_sediments_top_depth",
    "water_bottom_depth",
    "water_top_depth",
]


# ---- helpers --------------------------------------------------------------------------------------------------------
def put(arrays: dict[str, np.ndarray], name: str, value: Any) -> None:
    """Store a value under ``name`` as an array."""
    arrays[name] = np.asarray(value)


def put_frame(arrays: dict[str, np.ndarray], prefix: str, frame: pd.DataFrame) -> None:
    """Store every numeric column of a DataFrame and its index."""
    idx = frame.index
    if isinstance(idx, pd.DatetimeIndex):
        arrays[f"{prefix}{SEP}__index_ns__"] = idx.asi8
    else:
        arrays[f"{prefix}{SEP}__index__"] = np.asarray(idx, dtype=float)
    for col in frame.columns:
        values = frame[col].to_numpy()
        if values.dtype.kind in "fiuc":
            arrays[f"{prefix}{SEP}{col}"] = values


def put_tree(arrays: dict[str, np.ndarray], prefix: str, obj: Any) -> None:
    """Flatten nested dicts and lists of numbers into arrays."""
    if isinstance(obj, dict):
        for key, value in obj.items():
            put_tree(arrays, f"{prefix}{SEP}{key}", value)
    else:
        arrays[prefix] = np.asarray(obj, dtype=float)


@contextlib.contextmanager
def quiet() -> Any:
    """No log lines and no library warnings while a case runs."""
    logger.disable("scubas")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            yield
        finally:
            logger.enable("scubas")


@contextlib.contextmanager
def chdir(path: Path) -> Any:
    """Change directory for the duration (the LITHO1.0 cache lives in ``./.scubas_config``)."""
    old = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


def write_synthetic_litho(path: Path) -> None:
    """A small LITHO1.0-shaped netCDF (10 degree grid) with smooth, made-up depths in km; for code-path tests only."""
    from scipy.io import netcdf_file

    lat = np.arange(-90.0, 91.0, 10.0)
    lon = np.arange(-180.0, 181.0, 10.0)
    la, lo = np.meshgrid(lat, lon, indexing="ij")
    tilt = 0.01 * la + 0.002 * lo
    depth = {
        "water_top_depth": 0.0 * la,
        "water_bottom_depth": 2.0 + tilt,
        "upper_sediments_top_depth": 2.0 + tilt,
        "upper_sediments_bottom_depth": 2.5 + tilt,
        "middle_sediments_top_depth": 2.5 + tilt,
        "middle_sediments_bottom_depth": 3.0 + tilt,
        "lower_sediments_top_depth": 3.0 + tilt,
        "lower_sediments_bottom_depth": 3.4 + tilt,
        "upper_crust_top_depth": 3.4 + tilt,
        "upper_crust_bottom_depth": 8.0 + tilt,
        "middle_crust_top_depth": 8.0 + tilt,
        "middle_crust_bottom_depth": 15.0 + tilt,
        "lower_crust_top_depth": 15.0 + tilt,
        "lower_crust_bottom_depth": 30.0 + tilt,
        "lid_top_depth": 30.0 + tilt,
        "lid_bottom_depth": 120.0 + tilt,
        "asthenospheric_mantle_top_depth": 120.0 + tilt,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with netcdf_file(str(path), "w") as f:
        f.createDimension("latitude", len(lat))
        f.createDimension("longitude", len(lon))
        v = f.createVariable("latitude", "f8", ("latitude",))
        v[:] = lat
        v = f.createVariable("longitude", "f8", ("longitude",))
        v[:] = lon
        for name in LITHO_FIELDS:
            v = f.createVariable(name, "f8", ("latitude", "longitude"))
            v[:] = depth[name]


# ---- cases ----------------------------------------------------------------------------------------------------------
def case_sites(arrays: dict[str, np.ndarray], meta: dict[str, Any]) -> None:
    """Layered impedance and skin depth of every shipped profile; transfer functions of the ocean profiles."""
    from scubas.datasets import PROFILES, Layer, Site
    from scubas.models import OceanModel

    freqs = np.logspace(-6, 0, 25)
    for key, site in vars(PROFILES).items():
        z = site.calcZ(freqs)
        put(arrays, f"sites{SEP}{key}{SEP}Z", z)
        put(arrays, f"sites{SEP}{key}{SEP}P", site.calcP(freqs))
        put(arrays, f"sites{SEP}{key}{SEP}resistivities", site.get_resistivities())
    put_frame(arrays, f"sites{SEP}CS{SEP}get", PROFILES.CS.get())
    meta["site_type"] = Site.__name__
    meta["layer_resistivity"] = Layer("x", 1.0, 4.0).resistivity
    for key in ("CS", "SO", "DO", "OM"):
        model = OceanModel(getattr(PROFILES, key), flim=(1e-4, 1e-2))
        tfs = model.calcTF()
        put(arrays, f"tf{SEP}{key}{SEP}freqs", model.freqs)
        put(arrays, f"tf{SEP}{key}{SEP}E2B", tfs["E2B"])
        put(arrays, f"tf{SEP}{key}{SEP}B2B", tfs["B2B"])
        put(arrays, f"tf{SEP}{key}{SEP}Zocean", model.Z["ocean"])
        put(arrays, f"tf{SEP}{key}{SEP}Zfloor", model.Z["floor"])
    meta["default_flim_points"] = len(OceanModel(PROFILES.CS).freqs)


def case_sections(arrays: dict[str, np.ndarray], meta: dict[str, Any]) -> None:
    """Section geometry, line parameters, the equivalent-pi circuit and the full ingest from a B-field file."""
    from scubas.cables import CableSection, TransmissionLine
    from scubas.datasets import PROFILES

    for tag, directed in {
        "total": {"length": 4.0},
        "components": {"length_north": 3.0, "length_east": 4.0},
        "edges": {
            "edge_locations": {
                "initial": {"lat": 39.6, "lon": -74.33},
                "final": {"lat": 38.79, "lon": -72.62},
            }
        },
    }.items():
        sec = CableSection(f"len-{tag}", directed_length=directed)
        put(
            arrays,
            f"geometry{SEP}{tag}",
            [sec.length, sec.length_north, sec.length_east],
        )
    for key in ("CS", "SO", "DO", "LD"):
        tl = TransmissionLine(
            f"p-{key}",
            directed_length={"length": 10.0},
            elec_params={"site": getattr(PROFILES, key)},
        )
        put(arrays, f"line{SEP}{key}", [tl.C, tl.R, tl.Z, tl.Y, tl.gma, tl.Z0])
    meta["to_str"] = TransmissionLine("s", directed_length={"length": 10.0}).to_str()
    times = pd.date_range("2020-01-01", periods=5, freq="60s", name="Time")
    efield = pd.DataFrame(
        {"X": [1.0, 0.5, 0.25, -0.5, 2.0], "Y": [0.1, -0.2, 0.3, 0.0, 0.4]}, index=times
    )
    tl = TransmissionLine(
        "pi", directed_length={"length_north": 30.0, "length_east": 40.0}
    )
    tl.compute_eqv_pi_circuit(Efield=efield, components=["X", "Y"])
    for comp in ("X", "Y"):
        put(arrays, f"pi{SEP}Ye{SEP}{comp}", tl.Ye[comp])
        put(arrays, f"pi{SEP}Yp2{SEP}{comp}", tl.Yp2[comp])
        put(arrays, f"pi{SEP}Ie{SEP}{comp}", tl.Ie[comp])
    put_frame(arrays, f"pi{SEP}V", tl.V)
    tl2 = TransmissionLine(
        "ingest",
        directed_length={
            "edge_locations": {
                "initial": {"lat": 39.6, "lon": -74.33},
                "final": {"lat": 38.79, "lon": -72.62},
            }
        },
        elec_params={"site": PROFILES.CS_W, "width": 1.0, "flim": [1e-6, 1.0]},
    )
    tl2.compile_oml(
        bfield_data_files=[DATA / "tst20200101vmin.min"], csv_file_date_name="Date"
    )
    put_frame(arrays, f"ingest{SEP}Efield", tl2.Efield)
    put_frame(arrays, f"ingest{SEP}V", tl2.V)
    meta["ingest_components"] = tl2.components


def _cable_case(
    name: str,
    sections: list[Any],
    components: list[str],
    arrays: dict[str, np.ndarray],
    meta: dict[str, Any],
) -> None:
    from scubas.cables import Cable

    cable = Cable(sections, components=components)
    put_frame(arrays, f"{name}{SEP}tot_params", cable.tot_params)
    meta[f"{name}_tot_params_columns"] = list(cable.tot_params.columns)
    for comp in components:
        put(arrays, f"{name}{SEP}V{SEP}{comp}", cable.V[comp])
    put_tree(arrays, f"{name}{SEP}result", dict(cable.result.items()))
    v_along, l_along = cable._pot_along_cable_(timestamp=0)
    put(arrays, f"{name}{SEP}along_V", v_along)
    put(arrays, f"{name}{SEP}along_L", l_along)
    for i in range(len(sections)):
        u0, u1 = cable._pot_endCS_(i)
        put(arrays, f"{name}{SEP}endCS{SEP}{i}{SEP}U0", u0)
        put(arrays, f"{name}{SEP}endCS{SEP}{i}{SEP}U1", u1)
    u0, u1 = cable._pot_end_cable_(unit="mV")
    put(arrays, f"{name}{SEP}end_cable_mV{SEP}U0", u0)
    put(arrays, f"{name}{SEP}end_cable_mV{SEP}U1", u1)


def case_tutorial_cables(arrays: dict[str, np.ndarray], meta: dict[str, Any]) -> None:
    """The three conductivity tutorial cases (``examples/network_model.py``) driven by the 1989 FRD time axis."""
    from scubas.cables import TransmissionLine
    from scubas.datasets import PROFILES

    bfile = DATASET_1989 / "FRD.csv"
    index = pd.DatetimeIndex(
        pd.read_csv(bfile, parse_dates=["Date"])["Date"], name="Time"
    )

    def efield(mv_per_km: float) -> pd.DataFrame:
        return pd.DataFrame({"X": [mv_per_km] * len(index)}, index=index)

    def section(
        sec_id: str,
        length: float,
        profile: Any,
        ef: pd.DataFrame,
        terminations: Any = None,
    ) -> Any:
        tl = TransmissionLine(
            sec_id=sec_id,
            directed_length={"length_north": length},
            elec_params={"site": profile, "width": 1.0, "flim": [1e-6, 1.0]},
            active_termination=terminations,
        )
        tl.bfield_data_files = [bfile]
        tl.csv_file_date_name = "Date"
        tl.p = None
        tl.compute_eqv_pi_circuit(Efield=ef, components=["X"])
        return tl

    cs, ocean = efield(300.0), efield(100.0)
    cases = {
        "case1": [section("CS", 600.0, PROFILES.CS, cs)],
        "case2": [
            section("CS-west", 100.0, PROFILES.CS, cs, {"left": PROFILES.LD}),
            section("SO", 1800.0, PROFILES.SO, ocean),
            section("CS-east", 100.0, PROFILES.CS, cs, {"right": PROFILES.LD}),
        ],
        "case3": [
            section("CS-west", 100.0, PROFILES.CS, cs, {"left": PROFILES.LD}),
            section("DO", 7800.0, PROFILES.DO, ocean),
            section("CS-east", 100.0, PROFILES.CS, cs, {"right": PROFILES.LD}),
        ],
    }
    for name, sections in cases.items():
        _cable_case(f"tutorial_{name}", sections, ["X"], arrays, meta)


def case_event_1989(arrays: dict[str, np.ndarray], meta: dict[str, Any]) -> None:
    """The March 1989 TAT-8 study (``examples/network_1989.py``): nine sections, three stations, X and Y."""
    from scubas.cables import TransmissionLine
    from scubas.datasets import PROFILES

    definitions = [
        (
            "CS-W",
            PROFILES.CS_W,
            "FRD",
            {"left": PROFILES.LD},
            (39.60, -74.33, 38.79, -72.62),
        ),
        ("DO-1", PROFILES.DO_1, "FRD", None, (38.79, -72.62, 37.11, -68.94)),
        ("DO-2", PROFILES.DO_2, "STJ", None, (37.11, -68.94, 39.80, -48.20)),
        ("DO-3", PROFILES.DO_3, "STJ", None, (39.80, -48.20, 40.81, -45.19)),
        ("DO-4", PROFILES.DO_4, "STJ", None, (40.81, -45.19, 43.15, -39.16)),
        ("DO-5", PROFILES.DO_5, "STJ", None, (43.15, -39.16, 44.83, -34.48)),
        ("MAR", PROFILES.MAR, "STJ", None, (44.83, -34.48, 46.51, -22.43)),
        ("DO-6", PROFILES.DO_6, "HAD", None, (46.51, -22.43, 47.85, -9.05)),
        (
            "CS-E",
            PROFILES.CS_E,
            "HAD",
            {"right": PROFILES.LD},
            (47.85, -9.05, 50.79, -4.55),
        ),
    ]
    sections = []
    for sec_id, profile, station, terminations, (la0, lo0, la1, lo1) in definitions:
        tl = TransmissionLine(
            sec_id=sec_id,
            directed_length={
                "edge_locations": {
                    "initial": {"lat": la0, "lon": lo0},
                    "final": {"lat": la1, "lon": lo1},
                }
            },
            elec_params={"site": profile, "width": 1.0, "flim": [1e-6, 1.0]},
            active_termination=terminations,
        )
        tl.compile_oml(
            bfield_data_files=[DATASET_1989 / f"{station}.csv"],
            csv_file_date_name="Date",
        )
        tl.compute_eqv_pi_circuit(components=["X", "Y"])
        sections.append(tl)
        put(
            arrays,
            f"event1989{SEP}section{SEP}{sec_id}",
            [tl.length, tl.length_north, tl.length_east, tl.C, tl.R, tl.gma, tl.Z0],
        )
    _cable_case("event1989", sections, ["X", "Y"], arrays, meta)


def case_current(arrays: dict[str, np.ndarray], meta: dict[str, Any]) -> None:
    """The ionospheric line-current source over the Quebec profile."""
    from scubas.datasets import PROFILES
    from scubas.sources.current_sources import Current, Direction

    cur = Current(
        h=110.0,
        a=5.0,
        I=1.0,
        site=PROFILES.Quebec,
        current_direction=Direction.GEO_EAST_WEST,
    )
    x = np.linspace(-500.0, 500.0, len(cur.omega))
    cur.compute_B_field(x=x.copy(), length_unit="km")
    put(arrays, f"current{SEP}p", cur.p)
    for i, b in enumerate(cur.B_fields):
        for key in ("total", "external", "internal"):
            put(arrays, f"current{SEP}B{i}{SEP}{key}", b[key])
        meta[f"current_B{i}_direction"] = b["direction"].name
    for i, e in enumerate(cur.E_fields):
        put(arrays, f"current{SEP}E{i}{SEP}total", e["total"])
        meta[f"current_E{i}_direction"] = e["direction"].name
    meta["current_units"] = [cur.B_field_units, cur.E_field_units]


def case_utils(arrays: dict[str, np.ndarray], meta: dict[str, Any]) -> None:
    """Spectral helpers, component maps, distances, the taper and detrend."""
    from scubas.models import Preprocess
    from scubas.utils import (
        GreatCircle,
        component_mappings,
        component_sign_mappings,
        fft,
        frexp102str,
        ifft,
    )

    t = np.arange(256) * 60.0
    sig = (
        5.0 * np.sin(2 * np.pi * t / 3600.0)
        + 2.0 * np.cos(2 * np.pi * t / 900.0)
        + 0.001 * t
    )
    spectrum, freqs = fft(sig, 60.0)
    put(arrays, f"utils{SEP}fft{SEP}spectrum", spectrum)
    put(arrays, f"utils{SEP}fft{SEP}freqs", freqs)
    put(arrays, f"utils{SEP}ifft", ifft(spectrum))
    spectrum0, freqs0 = fft(sig, 60.0, remove_zero_frequency=False)
    put(arrays, f"utils{SEP}fft0{SEP}spectrum", spectrum0)
    put(arrays, f"utils{SEP}fft0{SEP}freqs", freqs0)
    pairs = [
        ((0.0, 0.0), (0.0, 90.0)),
        ((39.6, -74.33), (50.79, -4.55)),
        ((-33.9, 151.2), (51.5, -0.1)),
    ]
    for i, (a, b) in enumerate(pairs):
        gc = GreatCircle(
            SimpleNamespace(lat=a[0], lon=a[1]), SimpleNamespace(lat=b[0], lon=b[1])
        )
        put(arrays, f"utils{SEP}gc{SEP}{i}", [gc.great_circle(), gc.haversine()])
    meta["strings"] = [
        frexp102str(v, p) for v, p in ((0.0, 2), (1.234e5, 2), (-9.876e-4, 3), (5.0, 0))
    ]
    meta["maps"] = [
        component_mappings("B2E", "X"),
        component_mappings("B2E", "Y"),
        component_sign_mappings("BxEy"),
        component_sign_mappings("ByEx"),
    ]
    tt = np.arange(200.0)
    field = 10.0 + np.sin(tt / 7.0) + 0.02 * tt
    for p in (0.0, 0.1, 0.3):
        proc = Preprocess(tt, field, p)
        put(arrays, f"preprocess{SEP}taper{SEP}{p}", proc.get_tapering_function())
        put(arrays, f"preprocess{SEP}detrend{SEP}{p}", proc.detrend_magnetic_field())


def case_readers(arrays: dict[str, np.ndarray], meta: dict[str, Any]) -> None:
    """IAGA-2002 and CSV B-field readers, with data gaps, and the B-to-E transform of the result."""
    from scubas.datasets import PROFILES
    from scubas.models import OceanModel

    model = OceanModel(PROFILES.CS, flim=(1e-4, 1e-2))
    for name in ("tst20200101vmin.min", "tsthd20200101vmin.min"):
        frame, header = model.read_iaga(DATA / name, return_header=True)
        put_frame(arrays, f"iaga{SEP}{name}", frame)
        meta[f"iaga_{name}_header"] = {k: str(v) for k, v in header.items()}
        meta[f"iaga_{name}_columns"] = list(frame.columns)
    for name in ("tst20200101vmin.min", "gap.csv"):
        model = OceanModel(PROFILES.CS, flim=(1e-4, 1e-2))
        b = model.read_Bfield_data([DATA / name])
        put_frame(arrays, f"bfield{SEP}{name}", b)
        meta[f"bfield_{name}_columns"] = list(b.columns)
        model.to_Efields(p=0.1)
        put_frame(arrays, f"efield{SEP}{name}", model.Efield)
        meta[f"efield_{name}_columns"] = list(model.Efield.columns)


def case_conductivity(
    arrays: dict[str, np.ndarray], meta: dict[str, Any], workdir: Path
) -> None:
    """ConductivityProfile on a synthetic LITHO1.0-shaped grid (the cache path is relative to the working directory)."""
    from scubas.conductivity import ConductivityProfile

    litho = workdir / ".scubas_config" / "LITHO1.0.nc"
    write_synthetic_litho(litho)
    with chdir(workdir):
        cp = ConductivityProfile()
        meta["conductivity_defaults"] = dict(vars(cp.cprop).items())
        meta["conductivity_method"] = cp.grid_interpolation_method
        pts = cp.get_interpolation_points([30.0, -70.0], [35.0, -40.0])
        put(arrays, f"litho{SEP}points", pts)
        frame = cp._compile_profile_(np.array([30.0, -70.0]))
        put(arrays, f"litho{SEP}point_thickness", frame["thickness"].to_numpy())
        put(
            arrays,
            f"litho{SEP}point_resistivity",
            frame["resistivity"].to_numpy(dtype=float),
        )
        meta["litho_names"] = list(frame["name"])
        frame2 = cp._compile_profile_(pts)
        put(arrays, f"litho{SEP}line_thickness", frame2["thickness"].to_numpy())
        for kind in ("rounded", "exact"):
            site = ConductivityProfile.compile_profile(
                [30.4, -69.6], kind=kind, site_name="p", site_description="d"
            )
            put(
                arrays,
                f"litho{SEP}site{SEP}{kind}{SEP}thickness",
                site.get_thicknesses(),
            )
            put(
                arrays,
                f"litho{SEP}site{SEP}{kind}{SEP}conductivity",
                site.get_conductivities(),
            )
        raw = ConductivityProfile.compile_profile([30.0, -70.0], to_site=False)
        put(arrays, f"litho{SEP}raw{SEP}thickness", raw["thickness"].to_numpy())
        sites = ConductivityProfile.compile_profiles(
            [[30.0, -70.0], [35.0, -40.0]], site_names=["a", "b"]
        )
        put(
            arrays,
            f"litho{SEP}profiles{SEP}thickness",
            [s.get_thicknesses() for s in sites],
        )
        meta["litho_profile_names"] = [s.name for s in sites]
        binned = ConductivityProfile.compile_bined_profile(
            [[30.0, -70.0], [35.0, -40.0]]
        )
        put(arrays, f"litho{SEP}binned{SEP}thickness", binned.get_thicknesses())
        many = ConductivityProfile.compile_bined_profiles(
            np.array([[30.0, -70.0], [32.0, -55.0], [35.0, -40.0]])
        )
        put(
            arrays,
            f"litho{SEP}binned_many{SEP}thickness",
            [s.get_thicknesses() for s in many],
        )
        mc = ConductivityProfile.compile_mcmc_bined_profiles(
            [[30.0, -70.0], [32.0, -55.0], [35.0, -40.0]], n=3, random_seed=7
        )
        put(
            arrays,
            f"litho{SEP}mcmc{SEP}thickness",
            [[s.get_thicknesses() for s in draw] for draw in mc],
        )


BIG = 600  # arrays longer than this are stored thinned (every STRIDE-th value) plus summary statistics
STRIDE = 37


def compress(arrays: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Keep the baseline small: a long array becomes every 37th value plus [size, sum, sum of imaginary parts, sum of |x|, max |x|, NaN count]."""
    out: dict[str, np.ndarray] = {}
    for name, arr in arrays.items():
        if arr.size <= BIG or arr.dtype.kind not in "fc":
            out[name] = arr
            continue
        flat = arr.reshape(-1)
        out[name] = flat[::STRIDE]
        mag = np.abs(flat)
        out[f"{name}#stats"] = np.array(
            [
                flat.size,
                np.nansum(flat.real),
                np.nansum(flat.imag),
                np.nansum(mag),
                np.nanmax(mag) if np.isfinite(mag).any() else 0.0,
                np.isnan(mag).sum(),
            ]
        )
    return out


def build_cases(workdir: Path) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Run every case; ``workdir`` is an empty scratch directory (used for the LITHO1.0 cache)."""
    arrays: dict[str, np.ndarray] = {}
    meta: dict[str, Any] = {}
    with quiet():
        case_sites(arrays, meta)
        case_sections(arrays, meta)
        case_tutorial_cables(arrays, meta)
        case_event_1989(arrays, meta)
        case_current(arrays, meta)
        case_utils(arrays, meta)
        case_readers(arrays, meta)
        case_conductivity(arrays, meta, workdir)
    return compress(arrays), meta
