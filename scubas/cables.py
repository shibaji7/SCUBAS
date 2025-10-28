"""
Utilities for cable section analysis, including nodal studies and
equivalent transmission-line representations.
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

import json
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from loguru import logger

from scubas.datasets import PROFILES
from scubas.models import OceanModel
from scubas.utils import RecursiveNamespace, frexp102str


def _update_recursive_namespace(
    target: RecursiveNamespace, values: Mapping[str, Any]
) -> None:
    """
    Recursively update ``target`` namespace with values from ``values``.
    """
    for key, value in values.items():
        if isinstance(value, Mapping):
            current = getattr(target, key, None)
            if isinstance(current, RecursiveNamespace):
                _update_recursive_namespace(current, value)
            else:
                setattr(target, key, RecursiveNamespace(**value))
        else:
            setattr(target, key, value)


class CableSection:
    """
    Representation of a single cable section with geometric metadata.
    """

    def __init__(
        self,
        sec_id: str,
        directed_length: Optional[Mapping[str, Any]] = None,
    ) -> None:
        """
        Parameters
        ----------
        sec_id :
            Identifier for the cable section.
        directed_length :
            Optional mapping describing the segment geometry. Supported keys
            include ``length``, ``length_north``, ``length_east``, and the
            nested ``edge_locations`` structure with ``initial`` and ``final``
            latitude/longitude entries.
        """
        default_directed_length = {
            "length": None,
            "length_north": None,
            "length_east": None,
            "edge_locations": {
                "initial": {"lat": 0.0, "lon": 0.0},
                "final": {"lat": 0.0, "lon": 0.0},
            },
        }
        self.sec_id = sec_id
        self.directed_length = RecursiveNamespace(**default_directed_length)
        if directed_length:
            _update_recursive_namespace(self.directed_length, directed_length)

        self.length: float = 0.0
        self.length_north: float = 0.0
        self.length_east: float = 0.0
        self.components: List[str] = ["X", "Y"]
        self.cable_lengths: Dict[str, float] = {}
        self.compute_lengths()

    @staticmethod
    def check_location(loc: Any) -> bool:
        """
        Return ``True`` if a location-like object provides ``lat`` and ``lon``.
        """
        return bool(loc is not None and hasattr(loc, "lat") and hasattr(loc, "lon"))

    def compute_lengths(
        self, length_method: str = "great_circle", **legacy_kwargs: Any
    ) -> None:
        """
        Populate length attributes for the cable section.

        Parameters
        ----------
        length_method :
            Currently unused placeholder for future interpolation choices.
        **legacy_kwargs :
            Includes support for legacy ``type`` keyword to maintain API
            compatibility.
        """
        if "type" in legacy_kwargs and not length_method:
            length_method = legacy_kwargs["type"]  # noqa: F841

        dl = self.directed_length
        length_total = 0.0
        length_north = 0.0
        length_east = 0.0

        if getattr(dl, "length", None):
            length_total = float(dl.length)
            diag = length_total / np.sqrt(2.0)
            length_north = diag
            length_east = diag
        elif (
            getattr(dl, "length_north", None) is not None
            or getattr(dl, "length_east", None) is not None
        ):
            logger.info("Cable length derived from provided north/east components.")
            length_north = float(getattr(dl, "length_north", 0.0) or 0.0)
            length_east = float(getattr(dl, "length_east", 0.0) or 0.0)
            length_total = np.hypot(length_east, length_north)
        elif (
            hasattr(dl, "edge_locations")
            and self.check_location(dl.edge_locations.initial)
            and self.check_location(dl.edge_locations.final)
        ):
            lat0 = float(dl.edge_locations.initial.lat)
            lon0 = float(dl.edge_locations.initial.lon)
            lat1 = float(dl.edge_locations.final.lat)
            lon1 = float(dl.edge_locations.final.lon)
            lamb = 0.5 * (lat0 + lat1)
            length_north = (111.133 - 0.56 * np.cos(np.deg2rad(2 * lamb))) * abs(
                lat1 - lat0
            )
            length_east = (
                (111.5065 - 0.1872 * np.cos(np.deg2rad(2 * lamb)))
                * np.cos(np.deg2rad(lamb))
                * abs(lon0 - lon1)
            )
            length_total = np.hypot(length_east, length_north)
        else:
            logger.warning(
                "No cable edge information available for section '%s'; "
                "defaulting lengths to zero.",
                self.sec_id,
            )

        self.length = length_total
        self.length_north = length_north
        self.length_east = length_east
        self.cable_lengths = {"X": self.length_north, "Y": self.length_east}

    def _pot_alongCS_(
        self,
        Vi: Optional[np.ndarray] = None,
        Vk: Optional[np.ndarray] = None,
        ln: int = 1000,
        idx: Optional[int] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculate potentials along the cable section using a distributed model.

        Parameters
        ----------
        Vi, Vk :
            Optional arrays describing the potential at the initial/final node.
        ln :
            Number of subdivisions along the cable length.
        idx :
            Optional time index used when ``Vi``/``Vk`` are time series arrays.

        Returns
        -------
        tuple
            Potential profile (in volts) and distance samples (km).
        """
        if not hasattr(self, "end_pot"):
            raise RuntimeError("End potentials unavailable; run nodal analysis first.")

        logger.info("Potential along cable section '%s' at index %s", self.sec_id, idx)
        Vi = Vi if Vi is not None else self.end_pot.Vi
        Vk = Vk if Vk is not None else self.end_pot.Vk
        if idx is not None:
            Vi = self.end_pot.Vi[idx]
            Vk = self.end_pot.Vk[idx]

        L = self.length * 1e3
        if L == 0:
            raise ValueError("Cable section length is zero; cannot compute profile.")

        x = np.linspace(0, L, ln + 1)
        denom = np.exp(self.gma * L) - np.exp(-self.gma * L)
        if np.isclose(denom, 0.0):
            raise ValueError("Singular solution encountered for cable potentials.")

        V = ((Vk * np.exp(self.gma * L) - Vi) * np.exp(-self.gma * (L - x)) / denom) + (
            (Vi * np.exp(self.gma * L) - Vk) * np.exp(-self.gma * x) / denom
        )
        return V, x / 1.0e3


class TransmissionLine(CableSection):
    """
    Cable section with electrical properties based on a stratified Earth model.
    """

    def __init__(
        self,
        sec_id: str,
        directed_length: Optional[Mapping[str, Any]] = None,
        elec_params: Optional[Mapping[str, Any]] = None,
        active_termination: Optional[Mapping[str, Any]] = None,
    ) -> None:
        """
        Parameters
        ----------
        sec_id :
            Identifier for the cable section.
        directed_length :
            Optional geometric overrides (see :class:`CableSection`).
        elec_params :
            Optional mapping with `site`, `width`, and `flim` entries.
        active_termination :
            Optional mapping with `left` and `right` termination definitions.
        """
        super().__init__(sec_id, directed_length=directed_length)

        default_elec_params = {
            "site": PROFILES.CS,
            "width": 1.0,
            "flim": [1e-6, 1.0],
        }
        default_active_termination = {"right": None, "left": None}

        self.elec_params = RecursiveNamespace(**default_elec_params)
        if elec_params:
            _update_recursive_namespace(self.elec_params, elec_params)

        self.active_termination = RecursiveNamespace(**default_active_termination)
        if active_termination:
            _update_recursive_namespace(self.active_termination, active_termination)

        (
            self.C,
            self.R,
            self.Z,
            self.Y,
            self.gma,
            self.Z0,
        ) = self.calc_trasmission_line_parameters()
        self.end_pot = RecursiveNamespace()

    def to_str(self) -> str:
        """
        Return a formatted string summarising the transmission-line properties.
        """
        lines = [
            f"Z: {frexp102str(self.Z * 1e3)} (Ohm/km)",
            f"Y: {frexp102str(self.Y * 1e3)} (S/km)",
            f"Z0: {frexp102str(self.Z0)} (Ohm)",
            f"gma: {frexp102str(self.gma * 1e3)} (/km)",
            f"Ad: {frexp102str(1e-3 / self.gma)} (km)",
        ]
        return "\n".join(lines)

    def compile_oml(
        self,
        bfield_data_files: Optional[Sequence[Union[str, Path]]] = None,
        p: Optional[Sequence[float]] = None,
        csv_file_date_name: str = "Date",
    ) -> "TransmissionLine":
        """
        Instantiate and populate an :class:`OceanModel` for this section.

        Parameters
        ----------
        bfield_data_files :
            Optional sequence of geomagnetic data files.
        p :
            Optional frequency-domain smoothing parameters propagated to the
            ocean model.
        csv_file_date_name :
            Column name used for timestamps when ingesting CSV data.

        Returns
        -------
        TransmissionLine
            Self reference to enable fluent usage.

        Raises
        ------
        RuntimeError
            When B-field ingestion fails.
        """
        self.model = OceanModel(
            self.elec_params.site,
            flim=self.elec_params.flim,
        )
        if bfield_data_files:
            try:
                self.model.read_Bfield_data(
                    bfield_data_files, csv_file_date_name=csv_file_date_name
                )
            except Exception as exc:  # pragma: no cover - propagating model errors
                raise RuntimeError("Failed to ingest B-field data.") from exc
            self.model.to_Efields(p=p)
            self.compute_eqv_pi_circuit()
        return self

    def add_active_termination(self) -> None:
        """
        Update active termination definitions with derived admittance values.
        """
        if not hasattr(self, "Efield"):
            raise RuntimeError(
                "E-field values are unavailable; run 'compute_eqv_pi_circuit' first."
            )

        terminators = [
            self.active_termination.right,
            self.active_termination.left,
        ]
        for terminator in terminators:
            if terminator:
                C, R, Z, Y, gma, Z0 = self.calc_trasmission_line_parameters(
                    site=getattr(terminator, "site", None),
                    width=getattr(terminator, "width", None),
                )
                Jn: Dict[str, np.ndarray] = {}
                for component in self.components:
                    E = np.asarray(self.Efield[component]) * 1.0e-6
                    Jn[component] = E / Z
                terminator.Yn = 1.0 / Z0
                terminator.Jn = Jn
                terminator.Z0 = Z0
                terminator.R = R
                terminator.C = C
                terminator.Z = Z
                terminator.Y = Y
                terminator.gma = gma

    def calc_trasmission_line_parameters(
        self,
        site: Optional[Any] = None,
        width: Optional[float] = None,
    ) -> Tuple[float, float, float, float, float, float]:
        """
        Compute primary transmission-line parameters for the cable section.

        Parameters
        ----------
        site :
            Optional SCUBAS site description. Defaults to this section's
            configured site.
        width :
            Cable width (metres). Defaults to the configured width.

        Returns
        -------
        tuple
            Capacitance (m/Ohm), resistance (m*Ohm), series impedance (Ohm*m),
            shunt admittance (S/m), propagation constant (1/m), and characteristic
            impedance (Ohm).

        Raises
        ------
        ValueError
            If width is not positive.
        AttributeError
            When the supplied site object does not expose the expected API.
        """
        width = width if width is not None else self.elec_params.width
        site = site if site is not None else self.elec_params.site
        if width is None or width <= 0:
            raise ValueError("Cable width must be positive.")

        logger.info("Cable width %s: %s", self.sec_id, width)

        try:
            if getattr(site, "name", "") == "Land":
                C = width * ((site.get_thicknesses(0) / site.get_resistivities(0)))
                R = (
                    (site.get_thicknesses(1) * site.get_resistivities(1))
                    + (site.get_thicknesses(2) * site.get_resistivities(2))
                ) / width
            else:
                C = width * (
                    (site.get_thicknesses(1) / site.get_resistivities(1))
                    + (site.get_thicknesses(0) / site.get_resistivities(0))
                )
                R = (
                    (site.get_thicknesses(2) * site.get_resistivities(2))
                    + (site.get_thicknesses(3) * site.get_resistivities(3))
                ) / width
        except AttributeError as exc:
            raise AttributeError("Site object lacks required methods.") from exc

        Z = 1.0 / C
        Y = 1.0 / R
        gma = np.sqrt(Z * Y)
        Z0 = np.sqrt(Z / Y)
        return C, R, Z, Y, gma, Z0

    def compute_eqv_pi_circuit(
        self,
        Efield: Optional[pd.DataFrame] = None,
        components: Optional[Sequence[str]] = None,
    ) -> None:
        """
        Calculate equivalent pi-circuit parameters for the section.

        Parameters
        ----------
        Efield :
            Optional dataframe of electric fields in mV/km.
        components :
            Optional list of components; defaults to the model components.
        """
        Efield = Efield if Efield is not None else self.model.Efield
        components = components if components is not None else self.model.components

        self.Ye: Dict[str, np.ndarray] = {}
        self.Yp2: Dict[str, np.ndarray] = {}
        self.Ie: Dict[str, np.ndarray] = {}

        for component in components:
            L = self.cable_lengths.get(component)
            if L is None:
                raise KeyError(f"Component '{component}' missing cable length.")
            L_m = L * 1e3
            E = np.asarray(Efield[component]) * 1.0e-6
            sinh_term = np.sinh(self.gma * L_m)
            if np.isclose(sinh_term, 0.0):
                raise ValueError("Degenerate propagation constant encountered.")
            self.Ye[component] = 1.0 / (self.Z0 * sinh_term)
            self.Yp2[component] = (np.cosh(self.gma * L_m) - 1) * self.Ye[component]
            self.Ie[component] = E / self.Z

        self.Efield = Efield
        self.components = list(components)
        self.add_active_termination()
        self.compute_Vj(Efield.index.tolist())

    def compute_Vj(self, time: Sequence[Any]) -> None:
        """
        Compute the induced electric potential along the cable section.

        Parameters
        ----------
        time :
            Sequence of timestamps matching the electric-field samples.
        """
        self.V = pd.DataFrame({"Time": time, "Vj": 0.0})
        for component in self.components:
            length = self.cable_lengths.get(component, 0.0)
            self.V["Vj"] += np.asarray(self.Efield[component]) * length
        self.V = self.V.set_index("Time")


class Cable:
    """
    Aggregate representation of a multi-section cable system.
    """

    def __init__(
        self,
        cable_sections: Sequence[TransmissionLine],
        components: Sequence[str],
    ) -> None:
        if not cable_sections:
            raise ValueError("At least one cable section is required.")

        self.cable_sections = list(cable_sections)
        self.components = list(components)
        if not self.components:
            raise ValueError("At least one field component is required.")

        self.node_count = len(self.cable_sections) + 1
        self.node_ids = np.arange(self.node_count)
        self.left_edge, self.right_edge = 0, self.node_ids[-1]
        self.nodes: Dict[int, Dict[str, RecursiveNamespace]] = {}

        self.compile()

    def compile(self) -> None:
        """
        Execute nodal analysis and consolidate derived quantities.
        """
        self.run_nodal_analysis()
        self.solve_admitance_matrix()
        self.consolidate_final_result()
        U0, U1 = self._pot_end_cable_()

        self.tot_params = pd.DataFrame(
            {"Time": self.cable_sections[0].Efield.index.tolist()}
        )
        self.tot_params["V(v)"] = 0.0

        for component in self.components:
            self.tot_params[f"E.{component}"] = 0.0
            for idx, section in enumerate(self.cable_sections):
                column = f"E.{component}.{idx:02d}"
                self.tot_params[column] = np.asarray(section.Efield[component])
                self.tot_params[f"E.{component}"] += np.asarray(
                    section.Efield[component]
                )

        for idx, section in enumerate(self.cable_sections):
            column = f"V(v).{idx:02d}"
            self.tot_params[column] = np.asarray(section.V.Vj) / 1e3
            self.tot_params["V(v)"] += self.tot_params[column]

        self.tot_params["Vt(v)"] = U0 - U1 + np.asarray(self.tot_params["V(v)"])
        self.tot_params["U0"], self.tot_params["U1"] = U0, U1
        self.tot_params = self.tot_params.set_index("Time")

    def run_nodal_analysis(self) -> None:
        """
        Populate nodal admittance and current injections for each section.
        """
        logger.info("Running equivalent nodal analysis.")
        sections = self.cable_sections
        for nid in self.node_ids:
            self.nodes[nid] = {}
            logger.info("Node: %s", nid)
            for component in self.components:
                node = RecursiveNamespace()
                Yii = np.zeros_like(self.node_ids, dtype=float)

                if nid == self.left_edge:
                    Ji = -1.0 * sections[nid].Ie[component]
                    Yii[nid : nid + 2] = np.array(
                        [
                            sections[nid].Ye[component] + sections[nid].Yp2[component],
                            -sections[nid].Ye[component],
                        ]
                    )
                    if sections[nid].active_termination.left:
                        logger.info("Adding active termination: left")
                        Yii[nid] += sections[nid].active_termination.left.Yn
                        Ji = (
                            sections[nid].active_termination.left.Jn[component]
                            - sections[nid].Ie[component]
                        )
                elif nid == self.right_edge:
                    Ji = sections[-1].Ie[component]
                    Yii[nid - 1 : nid + 1] = np.array(
                        [
                            -sections[-1].Ye[component],
                            sections[-1].Yp2[component] + sections[-1].Ye[component],
                        ]
                    )
                    if sections[-1].active_termination.right:
                        logger.info("Adding active termination: right")
                        Yii[nid] += sections[-1].active_termination.right.Yn
                        Ji = Ji - sections[-1].active_termination.right.Jn[component]
                else:
                    Ji = sections[nid - 1].Ie[component] - sections[nid].Ie[component]
                    Yii[nid - 1 : nid + 2] = np.array(
                        [
                            -sections[nid - 1].Ye[component],
                            sections[nid - 1].Ye[component]
                            + sections[nid].Ye[component]
                            + sections[nid - 1].Yp2[component]
                            + sections[nid].Yp2[component],
                            -sections[nid].Ye[component],
                        ]
                    )

                node.Ji = Ji
                node.Yii = Yii
                self.nodes[nid][component] = node

    def solve_admitance_matrix(self) -> None:
        """
        Solve the nodal admittance system to recover potentials.
        """
        self.V: Dict[str, np.ndarray] = {}
        logger.info("Solving admittance matrix.")
        for component in self.components:
            logger.info("Solving for component %s.", component)
            J, Y = [], []
            for nid in self.node_ids:
                node = self.nodes[nid][component]
                J.append(node.Ji)
                Y.append(node.Yii)
            J_arr, Y_arr = np.array(J), np.array(Y)
            logger.info("Shapes -> J: %s, Y: %s", J_arr.shape, Y_arr.shape)
            try:
                iY = np.linalg.inv(Y_arr)
            except np.linalg.LinAlgError as exc:
                raise RuntimeError("Admittance matrix inversion failed.") from exc
            self.V[component] = np.matmul(iY, J_arr)

        for idx, section in enumerate(self.cable_sections):
            Vi = self.V[self.components[0]][idx, :]
            Vk = self.V[self.components[0]][idx + 1, :]
            if len(self.components) == 2:
                Vi = Vi + self.V[self.components[1]][idx, :]
                Vk = Vk + self.V[self.components[1]][idx + 1, :]
            section.end_pot.Vi = Vi
            section.end_pot.Vk = Vk

    def consolidate_final_result(self) -> None:
        """
        Collect per-section and per-node quantities into serialisable structures.
        """
        result = {"nodes": {}, "cables": {}}
        logger.info("Consolidating results.")
        for bid, section in enumerate(self.cable_sections, start=1):
            result["cables"][bid] = {
                "R": section.R,
                "C": section.C,
                "gma": section.gma,
                "Z0": section.Z0,
                "ln": section.length_north,
                "le": section.length_east,
                "len_km": section.length,
                "Ye": {},
                "Yp2": {},
                "Ie": {},
            }
            for component in self.components:
                result["cables"][bid]["Ye"][component] = section.Ye[component]
                result["cables"][bid]["Yp2"][component] = section.Yp2[component]
                result["cables"][bid]["Ie"][component] = section.Ie[component].tolist()

        for nid in self.node_ids:
            node_entry: Dict[str, Dict[str, List[float]]] = {}
            for component in self.components:
                node = self.nodes[nid][component]
                node_entry[component] = {
                    "Ji": node.Ji.tolist(),
                    "Yii": node.Yii.tolist(),
                }
            result["nodes"][str(nid)] = node_entry

        self.result = result

    def save(self, folder: Union[str, Path]) -> None:
        """
        Persist the nodal analysis artefacts to disk.

        Parameters
        ----------
        folder :
            Destination directory for outputs.

        Raises
        ------
        RuntimeError
            When writing output files fails.
        """
        output_dir = Path(folder)
        output_dir.mkdir(parents=True, exist_ok=True)
        try:
            (output_dir / "est_cable_props.json").write_text(
                json.dumps(self.result, sort_keys=True, indent=4)
            )
            self.tot_params.to_csv(output_dir / "sim-params.csv", float_format="%g")
        except OSError as exc:
            raise RuntimeError("Failed to persist cable analysis outputs.") from exc

    def _pot_endCS_byComp_(
        self,
        cable_section_id: int,
        comp: str,
        unit: str = "V",
        timestamp: Optional[int] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Voltage at both ends of a cable section for a given component.
        """
        factor = 1.0 if unit == "V" else 1000.0
        U0 = np.round(self.V[comp][cable_section_id, :] * factor, 2)
        U1 = np.round(self.V[comp][cable_section_id + 1, :] * factor, 2)
        logger.info(
            "Max(V) at the end of Section-%s (Component-%s): %s %s",
            cable_section_id,
            comp,
            np.max(U0),
            np.max(U1),
        )
        if timestamp is not None:
            U0, U1 = U0[timestamp], U1[timestamp]
        return U0, U1

    def _pot_endCS_(
        self, cable_section_id: int, unit: str = "V", timestamp: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Combined voltage at both ends of a cable section across components.
        """
        U0, U1 = self._pot_endCS_byComp_(
            cable_section_id, self.components[0], unit, timestamp
        )
        if len(self.components) == 2:
            u0, u1 = self._pot_endCS_byComp_(
                cable_section_id, self.components[1], unit, timestamp
            )
            U0 += u0
            U1 += u1
        return U0, U1

    def _pot_end_cable_byComp_(
        self, comp: str = "X", unit: str = "V", timestamp: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Voltage at both cable ends for a single component.
        """
        factor = 1.0 if unit == "V" else 1000.0
        U0 = np.round(self.V[comp][0, :] * factor, 2)
        U1 = np.round(self.V[comp][-1, :] * factor, 2)
        logger.info(
            "Max(V) at the cable ends (Component-%s): %s %s",
            comp,
            np.max(U0),
            np.max(U1),
        )
        if timestamp is not None:
            U0, U1 = U0[timestamp], U1[timestamp]
        return U0, U1

    def _pot_end_cable_(
        self, unit: str = "V", timestamp: Optional[int] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Combined voltage at both cable ends across components.
        """
        U0, U1 = self._pot_end_cable_byComp_(self.components[0], unit, timestamp)
        if len(self.components) == 2:
            u0, u1 = self._pot_end_cable_byComp_(self.components[1], unit, timestamp)
            U0 += u0
            U1 += u1
        return U0, U1

    def _pot_along_cable_(
        self, timestamp: Optional[int] = None, unit: str = "V"
    ) -> Tuple[List[float], List[float]]:
        """
        Voltage profile along the entire cable.
        """
        Vcable: List[float] = []
        Lcable: List[float] = []
        for idx, section in enumerate(self.cable_sections):
            V, Lx = section._pot_alongCS_(idx=timestamp)
            V_converted = V if unit == "V" else V * 1e3
            Vcable.extend(V_converted.tolist())
            if idx == 0:
                Lcable = Lx.tolist()
            else:
                Lcable.extend((Lx + Lcable[-1]).tolist())
        return Vcable, Lcable
