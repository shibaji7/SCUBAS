"""
    model.py: Module is used to implement Ocean Model with Layered-Earth structures.
    Class:
    ------
    OceanModel: Dedicated for ocean - Earth B and E field transfer functions
"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"


import datetime as dt

import numpy as np
import pandas as pd
from loguru import logger
from scipy import constants as C

from scubas.utils import component_mappings, component_sign_mappings, fft, ifft


class OceanModel(object):
    """
    This class is to emulate electric and magnetic fields
    in a 1D ocean. We provide the 1D Earth conductivity model and
    one layer of ocean parameters (depth and resistivity) to emulate
    1D ocean.
    """

    def __init__(
        self,
        site,
        flim=[1e-6, 1e0],
    ):
        """
        Initialize the ocean model
        """
        logger.info(f"Compile Ocean-model: {site.name} E- and B-Fields")
        self.site = site
        self.flim = flim
        self.freqs = np.linspace(flim[0], flim[1], int(flim[1] / flim[0]) + 1)
        self.functions = {
            "E2B": lambda Z, Zd, kd: Zd / (np.cosh(kd) + (Zd * np.sinh(kd) / Z)),
            "B2B": lambda Z, Zd, kd: 1.0 / (np.cosh(kd) + (Zd * np.sinh(kd) / Z)),
        }
        return

    def calcZ(self, freqs=None):
        """
        Compute the Km(f) for different layers
        """
        freqs = np.copy(self.freqs) if freqs is None else freqs
        omega = 2 * C.pi * freqs
        sigma_s = 1 / self.site.layers[0].resistivity
        k2 = 1.0j * omega * C.mu_0 * sigma_s
        k = np.sqrt(k2)
        Zo = (1.0j * omega * C.mu_0 / k) / (C.mu_0 / 1.0e-3)
        Kf = self.site.calcZ(freqs, layer=1)[1]
        self.Z = {"ocean": Zo, "floor": Kf}
        return

    def calcTF(self, kinds=["E2B", "B2B"], freqs=None):
        """
        Calculate the transfer functions.
        """
        freqs = np.copy(self.freqs) if freqs is None else freqs
        self.calcZ(freqs)
        omega = 2 * C.pi * freqs
        Zd = self.Z["floor"]
        Z = self.Z["ocean"]
        TFs = {}

        omega = 2 * C.pi * freqs
        sigma_s = 1 / self.site.layers[0].resistivity
        k2 = 1.0j * omega * C.mu_0 * sigma_s
        k = np.sqrt(k2)
        kd = k * self.site.layers[0].thickness

        for kind in kinds:
            TFs[kind] = self.functions[kind](Z, Zd, kd)
        return TFs

    def get_TFs(self, key="E2B", freqs=None):
        """
        Fetch the transfer function based on key
        and frequency.
        """
        TFs = self.calcTF(freqs=freqs)
        tf = pd.DataFrame()
        tf["freq"], tf[key] = np.copy(self.freqs) if freqs is None else freqs, TFs[key]
        return tf

    def read_iaga(self, file, return_xyzf=True, return_header=False):
        """
        Read IAGA profiles
        """

        # Read Headers
        header_records = {"header_length": 0}

        with open(file, "r") as openfile:
            for newline in openfile:
                if newline[0] == " ":
                    header_records["header_length"] += 1
                    label = newline[1:24].strip()
                    description = newline[24:-2].strip()
                    header_records[label.lower()] = description

        if len(header_records["reported"]) % 4 != 0:
            raise ValueError(
                "The header record does not contain 4 values: {0}".format(
                    header_records["reported"]
                )
            )
        record_length = len(header_records["reported"]) // 4
        column_names = [
            x for x in header_records["reported"][record_length - 1 :: record_length]
        ]
        seen_count = {}
        for i, col in enumerate(column_names):
            if col in seen_count:
                column_names[i] += str(seen_count[col])
                seen_count[col] += 1
            else:
                seen_count[col] = 1
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
            # Convert the data to XYZF format
            # Only convert HD
            if "H" not in column_names or "D" not in column_names:
                raise ValueError(
                    "Only have a converter for HDZF->XYZF\n"
                    + "Input file is: "
                    + header_records["reported"]
                )

            # IAGA-2002 D is reported in minutes of arc.
            df["X"] = df["H"] * np.cos(np.deg2rad(df["D"] / 60.0))
            df["Y"] = df["H"] * np.sin(np.deg2rad(df["D"] / 60.0))
            del df["H"], df["D"]
        if return_header:
            return df, header_records
        else:
            return df

    def read_Bfield_data(self, files, return_xyzf=True, csv_file_date_name="Date"):
        """
        Read B-Files
        """
        self.Bfield = pd.DataFrame()
        for file in files:
            file_type = file.split(".")[-1]
            if file_type == "txt":
                o = self.read_iaga(file, return_xyzf, return_header=False)
            elif file_type == "csv":
                o = pd.read_csv(file, parse_dates=[csv_file_date_name])
                o = o.rename(columns={csv_file_date_name: "Time"})
                o = o.set_index("Time")
                o.index.name = "Time"
            self.Bfield = pd.concat([self.Bfield, o])
        return self.Bfield

    def to_Efields(self, Bfield=None, components=["X", "Y"], p=None):
        """
        Compute Et using numerical FFT and IFFT block
        Bfield: Dataframe containing B-field
        components: [X and Y] for B-fields
        """
        self.Bfield = Bfield if Bfield else self.Bfield
        self.components = components
        self.Efield = pd.DataFrame()
        self.Efield["Time"] = self.Bfield.index.tolist()
        tx = self.Efield.Time.tolist()[0]
        if isinstance(tx, dt.datetime) or isinstance(tx, pd.Timestamp):
            t = np.array(self.Efield.Time.apply(lambda x: (x - tx).total_seconds()))
            self.Efield["dTime"], self.Bfield["dTime"] = t, t
        else:
            t = np.array(self.Efield.Time)
            self.Efield["dTime"], self.Bfield["dTime"] = (t[1] - t[0], t[1] - t[0])
        for a in self.components:
            Bt = np.array(self.Bfield[a])
            proc = Preprocess(np.array(self.Bfield.dTime), np.array(Bt))
            Bt = proc.detrend_magnetic_field(p)
            dT = (
                self.Bfield.dTime.tolist()[1]
                if "dTime" in self.Bfield.columns
                else t[1] - t[0]
            )
            Bf, f = fft(Bt, dT)
            E2B = np.array(self.get_TFs(freqs=f).E2B)
            Et = ifft(
                component_sign_mappings(
                    "B%sE%s" % (a.lower(), component_mappings("B2E", a).lower())
                )
                * E2B
                * Bf
            )
            self.Efield[component_mappings("B2E", a)] = Et
        self.Efield = self.Efield.set_index("Time")
        self.components = [component_mappings("B2E", a) for a in self.components]
        return


class Preprocess(object):
    """
    This class is used to process the B-field and E-field data.
    """

    def __init__(self, t, field, p=0.1):
        """
        Parameters:
        -----------
        t: Array of time seconds
        field: Array containing the data
        p: Tapering value (0-1)
        """
        self.t = t
        self.field = field
        self.p = p
        return

    def get_tapering_function(self, p=None):
        """
        This method is resposible for generateing
        tapering function based on time sequence t
        and tapering coefficient p
        """
        p = self.p if p is None else p
        T = len(self.t)
        P, P2 = int(T * p), int(T * p / 2)
        w = np.zeros_like(self.t)
        w[:P2] = 0.5 * (1 - np.cos(2 * np.pi * self.t[:P2] / P))
        w[P2 : T - P2] = 1.0
        w[T - P2 :] = 0.5 * (
            1 - np.cos(2 * np.pi * (self.t[-1] - self.t[T - P2 :]) / P)
        )
        return w

    def detrend_magnetic_field(self, p=None):
        """
        This method is resposible for detrend
        magnetic field data and taper it to reduce
        spurious frequency components.
        """
        p = self.p if p is None else p
        fmed = np.median(self.field[:120])
        w = self.get_tapering_function(p)
        f = (self.field - fmed) * w
        return f
