"""
Transfer-function utilities for electric and magnetic fields over layered oceans.
"""

from __future__ import annotations

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "chakras4@erau.edu"
__status__ = "Research"

import datetime as dt
from dataclasses import dataclass
from pathlib import Path
from typing import (
    Iterable,
    List,
    Mapping,
    MutableMapping,
    Optional,
    Sequence,
    Tuple,
    Union,
)

import numpy as np
import pandas as pd
from loguru import logger
from scipy import constants as C

from scubas.datasets import Site
from scubas.utils import component_mappings, component_sign_mappings, fft, ifft


def _validate_frequency_limits(flim: Sequence[float]) -> Tuple[float, float]:
    """
    Normalise and validate frequency bounds.
    """
    if len(flim) != 2:
        raise ValueError("Frequency limits must contain exactly two entries.")
    fmin, fmax = float(flim[0]), float(flim[1])
    if fmin <= 0 or fmax <= 0:
        raise ValueError("Frequency limits must be positive.")
    if fmin >= fmax:
        raise ValueError("Frequency limits must be strictly increasing.")
    return fmin, fmax


class OceanModel:
    """
    Emulates electric and magnetic fields for a layered ocean profile.
    """

    def __init__(
        self,
        site: Site,
        flim: Sequence[float] = (1e-6, 1e0),
    ) -> None:
        fmin, fmax = _validate_frequency_limits(flim)
        logger.info(f"Compile Ocean-model: {site.name} E- and B-Fields")
        self.site = site
        self.flim = (fmin, fmax)
        self.freqs = np.linspace(
            fmin,
            fmax,
            int(round(fmax / fmin)) + 1,
        )
        self.functions: Mapping[str, callable] = {
            "E2B": lambda Z, Zd, kd: Zd / (np.cosh(kd) + (Zd * np.sinh(kd) / Z)),
            "B2B": lambda Z, Zd, kd: 1.0 / (np.cosh(kd) + (Zd * np.sinh(kd) / Z)),
        }
        self.Z: MutableMapping[str, np.ndarray] = {}

    def calcZ(self, freqs: Optional[Iterable[float]] = None) -> None:
        """
        Compute ocean and seabed impedances for the provided frequencies.
        """
        freqs_arr = np.asarray(freqs, dtype=float) if freqs is not None else self.freqs
        omega = 2 * np.pi * freqs_arr
        seawater_resistivity = self.site.layers[0].resistivity
        sigma_s = 1.0 / seawater_resistivity
        k = np.sqrt(1.0j * omega * C.mu_0 * sigma_s)
        Zo = (1.0j * omega * C.mu_0 / k) / (C.mu_0 / 1.0e-3)
        Kf = self.site.calcZ(freqs_arr, layer=1)[1]
        self.Z = {"ocean": Zo, "floor": Kf}

    def calcTF(
        self,
        kinds: Sequence[str] = ("E2B", "B2B"),
        freqs: Optional[Iterable[float]] = None,
    ) -> Mapping[str, np.ndarray]:
        """
        Calculate transfer functions for the supplied kind(s) and frequency grid.
        """
        freqs_arr = np.asarray(freqs, dtype=float) if freqs is not None else self.freqs
        self.calcZ(freqs_arr)
        omega = 2 * np.pi * freqs_arr
        Zd = self.Z["floor"]
        Z = self.Z["ocean"]
        sigma_s = 1.0 / self.site.layers[0].resistivity
        k = np.sqrt(1.0j * omega * C.mu_0 * sigma_s)
        kd = k * self.site.layers[0].thickness
        transfer_functions = {}
        for kind in kinds:
            if kind not in self.functions:
                raise KeyError(f"Unknown transfer function kind '{kind}'.")
            transfer_functions[kind] = self.functions[kind](Z, Zd, kd)
        return transfer_functions

    def get_TFs(
        self,
        key: str = "E2B",
        freqs: Optional[Iterable[float]] = None,
    ) -> pd.DataFrame:
        """
        Retrieve the specified transfer function as a dataframe.
        """
        tf_values = self.calcTF(freqs=freqs)
        freqs_arr = np.asarray(freqs, dtype=float) if freqs is not None else self.freqs
        if key not in tf_values:
            raise KeyError(f"Transfer function '{key}' not available.")
        frame = pd.DataFrame({"freq": freqs_arr, key: tf_values[key]})
        return frame

    def read_iaga(
        self,
        file: Union[str, Path],
        return_xyzf: bool = True,
        return_header: bool = False,
    ) -> Union[pd.DataFrame, Tuple[pd.DataFrame, Mapping[str, str]]]:
        """
        Parse an IAGA-2002 geomagnetic file into a dataframe.
        """
        header_records: MutableMapping[str, Union[str, int]] = {"header_length": 0}

        with open(file, "r", encoding="utf-8") as openfile:
            for line in openfile:
                if line and line[0] == " ":
                    header_records["header_length"] += 1
                    label = line[1:24].strip()
                    description = line[24:-2].strip()
                    header_records[label.lower()] = description

        reported = header_records.get("reported")
        if not reported:
            raise ValueError("IAGA header missing 'reported' column descriptors.")
        if len(reported) % 4 != 0:
            raise ValueError(f"The header record does not contain 4 values: {reported}")
        record_length = len(reported) // 4
        column_names = list(reported[record_length - 1 :: record_length])
        seen_count: MutableMapping[str, int] = {}
        for idx, column in enumerate(column_names):
            if column in seen_count:
                column_names[idx] += str(seen_count[column])
                seen_count[column] += 1
            else:
                seen_count[column] = 1

        df = pd.read_csv(
            file,
            header=header_records["header_length"],
            delim_whitespace=True,
            parse_dates=[[0, 1]],
            infer_datetime_format=True,
            index_col=0,
            usecols=[0, 1, 3, 4, 5, 6],
            na_values=[99999.90, 99999.0, 88888.80, 88888.00],
            names=["Date", "Time"] + column_names,
        )
        df.index.name = "Time"

        if return_xyzf and "X" not in column_names and "Y" not in column_names:
            if "H" not in column_names or "D" not in column_names:
                raise ValueError(
                    "Only HDZF to XYZF conversion is supported, but H or D is missing."
                )
            df["X"] = df["H"] * np.cos(np.deg2rad(df["D"] / 60.0))
            df["Y"] = df["H"] * np.sin(np.deg2rad(df["D"] / 60.0))
            df = df.drop(columns=["H", "D"])

        return (df, header_records) if return_header else df

    def read_Bfield_data(
        self,
        files: Sequence[Union[str, Path]],
        return_xyzf: bool = True,
        csv_file_date_name: str = "Date",
        check_for_gaps: bool = True,
    ) -> pd.DataFrame:
        """
        Aggregate B-field data from IAGA or CSV files.
        """
        data_frames: List[pd.DataFrame] = []
        for file in files:
            logger.info(f"Read file: {file}")
            file_suffix = Path(file).suffix.lower()
            if file_suffix in [".txt", ".min"]:
                df = self.read_iaga(file, return_xyzf, return_header=False)
            elif file_suffix == ".csv":
                df = pd.read_csv(file, parse_dates=[csv_file_date_name])
                df = df.rename(columns={csv_file_date_name: "Time"}).set_index("Time")
                df.index.name = "Time"
            else:
                raise ValueError(f"Unsupported B-field file extension '{file_suffix}'.")
            data_frames.append(df)
        self.Bfield = pd.concat(data_frames, axis=0).sort_index()
        if check_for_gaps:
            for component in ["X", "Y", "Z", "F"]:
                if component not in self.Bfield.columns:
                    raise ValueError(
                        f"B-field component '{component}' missing in data."
                    )
                else:
                    all_gaps = self.Bfield[component].isnull().sum()
                    logger.info(f"B-field[{component}] has {all_gaps} data gaps.")
                    if all_gaps > 0:
                        self.Bfield[component] = self.Bfield[component].interpolate(
                            method="spline", order=2
                        )
        logger.info(f"Compiled B-field data with {len(self.Bfield)} records.")

        return self.Bfield

    def to_Efields(
        self,
        Bfield: Optional[pd.DataFrame] = None,
        components: Sequence[str] = ("X", "Y"),
        p: Optional[float] = None,
    ) -> None:
        """
        Transform B-field measurements into E-field estimates.
        """
        self.Bfield = Bfield if Bfield is not None else getattr(self, "Bfield", None)
        if self.Bfield is None:
            raise RuntimeError(
                "B-field data not available; load data before converting."
            )
        self.components = list(components)

        self.Efield = pd.DataFrame({"Time": self.Bfield.index})
        first_time = self.Efield.Time.iloc[0]
        if isinstance(first_time, (dt.datetime, pd.Timestamp)):
            dt_seconds = self.Efield.Time.apply(
                lambda x: (x - first_time).total_seconds()
            )
            self.Efield["dTime"] = dt_seconds.tolist()
            self.Bfield["dTime"] = dt_seconds.tolist()
        else:
            dt_value = np.asarray(self.Efield.Time, dtype=float)
            delta = dt_value[1] - dt_value[0] if len(dt_value) > 1 else 0.0
            self.Efield["dTime"] = delta
            self.Bfield["dTime"] = delta
        for component in self.components:
            Bt = np.asarray(self.Bfield[component])
            proc = Preprocess(np.asarray(self.Bfield.dTime), np.asarray(Bt))
            Bt_proc = proc.detrend_magnetic_field(p)
            delta_t = (
                self.Bfield.dTime.iloc[1]
                if "dTime" in self.Bfield.columns and len(self.Bfield) > 1
                else proc.t[1] - proc.t[0]
            )
            Bf, freqs = fft(Bt_proc, delta_t)
            tf_values = np.asarray(self.get_TFs(freqs=freqs).E2B)
            mapped_component = component_mappings("B2E", component)
            sign = component_sign_mappings(
                f"B{component.lower()}E{mapped_component.lower()}"
            )
            Et = ifft(sign * tf_values * Bf)
            self.Efield[mapped_component] = Et

        self.Efield = self.Efield.set_index("Time")
        self.components = [
            component_mappings("B2E", component) for component in self.components
        ]


@dataclass
class Preprocess:
    """
    Preprocessing helper for tapering and detrending field data.
    """

    t: np.ndarray
    field: np.ndarray
    p: float = 0.1

    def get_tapering_function(self, p: Optional[float] = None) -> np.ndarray:
        """
        Construct the tapering window applied during detrending.
        """
        p = self.p if p is None else p
        if not 0 <= p <= 1:
            raise ValueError("Taper value 'p' must lie between 0 and 1.")
        T = len(self.t)
        P = int(T * p)
        if P == 0:
            return np.ones_like(self.t)
        P2 = max(int(P / 2), 1)
        window = np.zeros_like(self.t, dtype=float)
        window[:P2] = 0.5 * (1 - np.cos(2 * np.pi * self.t[:P2] / P))
        window[P2 : T - P2] = 1.0
        window[T - P2 :] = 0.5 * (
            1 - np.cos(2 * np.pi * (self.t[-1] - self.t[T - P2 :]) / P)
        )
        return window

    def detrend_magnetic_field(self, p: Optional[float] = None) -> np.ndarray:
        """
        Remove bias and taper magnetic field data to reduce spectral leakage.
        """
        p = self.p if p is None else p
        if len(self.field) == 0:
            raise ValueError("Field data is empty; cannot detrend.")
        fmed = np.median(self.field[: min(120, len(self.field))])
        window = self.get_tapering_function(p)
        return (self.field - fmed) * window
