"""
    cables.py: Module is used to implement cable section analysis and an event study
"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"


import numpy as np
import pandas as pd
from loguru import logger

from scubas.datasets import PROFILES
from scubas.models import OceanModel
from scubas.utils import RecursiveNamespace, frexp102str


class CableSection(object):
    """
    This class holds a cable section of a big cable
    Parameters:
    -----------
    """

    def __init__(
        self,
        sec_id,
        directed_length=dict(
            length=None,
            length_north=None,
            length_east=None,
            edge_locations=dict(
                initial=dict(lat=0.0, lon=0.0), final=dict(lat=0.0, lon=0.0)
            ),
        ),
    ):
        self.sec_id = sec_id
        self.directed_length = RecursiveNamespace(**directed_length)
        self.components = []
        self.compute_lengths()
        return

    def check_location(self, loc):
        """
        Check lat/lon exists in file or not
        """
        tag = True if (hasattr(loc, "lat") and hasattr(loc, "lon")) else False
        return tag

    def compute_lengths(self):
        """
        Compute all the length of the Cable Section
        """
        dl = self.directed_length
        if hasattr(dl, "length") and dl.length:
            self.length_north, self.length_east = (
                dl.length / np.sqrt(2),
                dl.length / np.sqrt(2),
            )
        elif (hasattr(dl, "length_north") and dl.length_north is not None) or (
            hasattr(dl, "length_east") and dl.length_east is not None
        ):
            logger.info("Cable length from directed length")
            self.length_north = dl.length_north if hasattr(dl, "length_north") else 0.0
            self.length_east = dl.length_east if hasattr(dl, "length_east") else 0.0
            self.length = np.sqrt(self.length_east**2 + self.length_north**2)
        elif (
            hasattr(dl, "edge_locations")
            and self.check_location(dl.edge_locations.initial)
            and self.check_location(dl.edge_locations.final)
        ):
            # PARAMETERS OF THE WGS84 EARTH MODEL
            a = 6378.137  # Equatorial radius
            b = 6356.752  # Polar Radius
            e = np.sqrt(0.00669437999014)  # Eccentricity
            lamb = 0.5 * (
                self.directed_length.edge_locations.initial.lat
                + self.directed_length.edge_locations.final.lat
            )
            self.length_north = (
                111.133 - 0.56 * np.cos(np.deg2rad(2 * lamb))
            ) * np.abs(
                self.directed_length.edge_locations.final.lat
                - self.directed_length.edge_locations.initial.lat
            )
            self.length_east = (
                (111.5065 - 0.1872 * np.cos(np.deg2rad(2 * lamb)))
                * np.cos(np.deg2rad(lamb))
                * np.abs(
                    self.directed_length.edge_locations.initial.lon
                    - self.directed_length.edge_locations.final.lon
                )
            )
            self.length = np.sqrt(self.length_north**2 + self.length_east**2)
        else:
            logger.warning("No cable edge information available")
        self.components = ["X", "Y"]
        self.cable_lengths = {"X": self.length_north, "Y": self.length_east}
        return


class TransmissionLine(CableSection):
    """
    This class is dedicated for DSTL.
    Parameters:
    -----------
    """

    def __init__(
        self,
        sec_id,
        directed_length=dict(
            length=None,
            length_north=None,
            length_east=None,
            edge_locations=dict(
                initial=dict(lat=0.0, lon=0.0), final=dict(lat=0.0, lon=0.0)
            ),
        ),
        elec_params=dict(
            site=PROFILES.CS,
            width=1.0,
            flim=[1e-6, 1e0],
        ),
        active_termination=dict(
            right=None,
            left=None,
        ),
    ):
        """
        Properties:
        -----------
        """
        # Compute phyiscal properties of the cable section
        super().__init__(sec_id, directed_length=directed_length)
        self.elec_params = RecursiveNamespace(**elec_params)
        self.active_termination = RecursiveNamespace(**active_termination)
        # Extract electrical properties of the cable
        (
            self.C,
            self.R,
            self.Z,
            self.Y,
            self.gma,
            self.Z0,
        ) = self.calc_trasmission_line_parameters()
        self.end_pot = RecursiveNamespace(**dict())
        return

    def to_str(self):
        """
        Create a string of properties for display
        """
        o = "Z: %s (Ohm/km)\n" % (frexp102str(self.Z * 1e3))
        o += "Y: %s (S/km)\n" % (frexp102str(self.Y * 1e3))
        o += "Z0: %s (Ohm)\n" % (frexp102str(self.Z0))
        o += "gma: %s (/km)\n" % (frexp102str(self.gma * 1e3))
        o += "Ad: %s (km)" % (frexp102str(1e-3 / self.gma))
        return o

    def compile_oml(self, bfield_data_files=[], p=None, csv_file_date_name="Date"):
        """
        Create ocean model
        """
        self.model = OceanModel(
            self.elec_params.site,
            flim=self.elec_params.flim,
        )
        if bfield_data_files and len(bfield_data_files) > 0:
            self.model.read_Bfield_data(
                bfield_data_files, csv_file_date_name=csv_file_date_name
            )
            self.model.to_Efields(p=p)
            self.compute_eqv_pi_circuit()
        return self

    def add_active_termination(self):
        """
        Adding active terminations
        """
        terminators = [
            self.active_termination.right,
            self.active_termination.left,
        ]
        for at in terminators:
            if at:
                C, R, Z, Y, gma, Z0 = self.calc_trasmission_line_parameters(at)
                Jn = dict()
                for a in self.components:
                    E = np.array(self.Efield[a]) * 1.0e-3 / 1.0e3
                    Jn[a] = E / Z  # Assuming input mV/km convert to V/m
                setattr(at, "Yn", 1.0 / Z0)
                setattr(at, "Jn", Jn)
                setattr(at, "Z0", Z0)
                setattr(at, "R", R)
                setattr(at, "C", C)
                setattr(at, "Z", Z)
                setattr(at, "Y", Y)
                setattr(at, "gma", gma)
        return

    def calc_trasmission_line_parameters(self, site=None, width=None):
        """
        Compute the transmission line parameters
        """
        width = width if width else self.elec_params.width
        site = site if site else self.elec_params.site
        logger.info(f"Cable width: {width}")
        if site.name == "Land":
            C = width * ((site.get_thicknesses(0) / site.get_resistivities(0)))
            R = (
                (site.get_thicknesses(1) * site.get_resistivities(1))
                + (site.get_thicknesses(2) * site.get_resistivities(2))
            ) / width
        else:
            C = width * (
                (site.get_thicknesses(1) / site.get_resistivities(1))
                + (site.get_thicknesses(0) / site.get_resistivities(0))
            )  # in m/ohm
            R = (
                (site.get_thicknesses(2) * site.get_resistivities(2))
                + (site.get_thicknesses(3) * site.get_resistivities(3))
            ) / width  # in m*ohm
        Z, Y = 1.0 / C, 1.0 / R  # in Ohm-m and S/m
        gma, Z0 = np.sqrt(Z * Y), np.sqrt(Z / Y)  # in /m and Ohm
        return C, R, Z, Y, gma, Z0

    def compute_eqv_pi_circuit(self, Efield=None, components=None):
        """
        Calculate equivalent pi circuit model.
        X component is Nort (n), Y component East (e)
        dE: Dataframe containing E-field
        components: [X and Y] for E-fields
        """
        Efield = Efield if Efield is not None else self.model.Efield
        components = components if components is not None else self.model.components
        self.Ye, self.Yp2, self.Ie = {}, {}, {}
        for a in components:
            L = self.cable_lengths[a]
            L *= 1000.0  # Convert km to m
            E = (
                np.array(Efield[a]) * 1.0e-3 / 1.0e3
            )  # Assuming input mV/km convert to V/m
            self.Ye[a] = 1.0 / (self.Z0 * np.sinh(self.gma * L))
            self.Yp2[a] = (np.cosh(self.gma * L) - 1) * self.Ye[a]
            self.Ie[a] = E / self.Z
        self.Efield = Efield
        self.components = components
        self.add_active_termination()
        self.compute_Vj(Efield.index.tolist())
        return

    def compute_Vj(self, time):
        """
        Calculate total electric potential induced along the cable segment.
        Vj = Ej_n(t)Lj_n + Ej_e(t)Lj_e
        """
        self.V = pd.DataFrame()
        self.V["Time"] = time
        self.V["Vj"] = 0.0
        for a in self.components:
            lx = self.cable_lengths[a]
            self.V["Vj"] += (
                np.array(self.Efield[a]) * lx
            )  # Potential in mV: E(mV/km) length: km
        self.V = self.V.set_index("Time")
        return

    def _pot_alongCS_(self, Vi=None, Vk=None, ln=1000):
        """
        Caclulate potentials along the cable section
        """
        Vi = Vi if Vi is not None else self.end_pot.Vi
        Vk = Vk if Vk is not None else self.end_pot.Vk
        L = self.length * 1e3
        x = np.linspace(0, L, ln + 1)
        V = (
            (Vk * np.exp(self.gma * L) - Vi)
            * np.exp(-self.gma * (L - x))
            / (np.exp(self.gma * L) - np.exp(-self.gma * L))
        ) + (
            (Vi * np.exp(self.gma * L) - Vk)
            * np.exp(-self.gma * x)
            / (np.exp(self.gma * L) - np.exp(-self.gma * L))
        )
        return V, x / 1.0e3


class Cable(object):
    """
    This class holds a cable
    Parameters:
    -----------
    cable: Cable parameters
    Efields: Dataframe of E field
    Bfields: Dataframe for B Field
    components: Components for B or E fields
    """

    def __init__(
        self,
        cable_sections,
        components,
    ):
        self.cable_sections = cable_sections
        self.components = components
        self.nodes = len(self.cable_sections) + 1
        self.node_ids = np.arange(self.nodes)
        self.left_edge, self.right_edge = 0, self.node_ids[-1]
        self.nodes = {}
        self.compile()
        return

    def compile(self):
        """
        Run nodal analysis for the cable
        """
        self.run_nodal_analysis()
        self.solve_admitance_matrix()
        self.consolidate_final_result()
        U0, U1 = self._pot_end_cable_()
        # Total parameter calculations
        self.tot_params = pd.DataFrame()
        self.tot_params["Time"] = self.cable_sections[0].Efield.index.tolist()
        self.tot_params["V(v)"] = 0.0
        for a in self.components:
            self.tot_params["E." + a] = 0.0
            for i, tl in enumerate(self.cable_sections):
                self.tot_params["E.%s.%02d" % (a, i)] = np.array(tl.Efield[a])
                self.tot_params["E." + a] += np.array(tl.Efield[a])
        for i, tl in enumerate(self.cable_sections):
            self.tot_params["V(v).%02d" % (i)] = np.array(tl.V.Vj) / 1e3
            self.tot_params["V(v)"] += np.array(tl.V.Vj) / 1e3

        self.tot_params["Vt(v)"] = U0 - U1 + np.array(self.tot_params["V(v)"])
        self.tot_params["U0"], self.tot_params["U1"] = U0, U1
        self.tot_params = self.tot_params.set_index("Time")
        return

    def run_nodal_analysis(self):
        logger.info(f"Eq. nodal analysis.")
        sections = self.cable_sections
        for nid in self.node_ids:
            self.nodes[nid] = {}
            logger.info(f"Node:{nid}")
            for a in self.components:
                node = RecursiveNamespace(**dict())
                Yii = np.zeros_like(self.node_ids, dtype=float)
                if nid == self.left_edge:
                    Ji = -1.0 * sections[nid].Ie[a]
                    Yii[nid : nid + 2] = np.array(
                        [
                            sections[nid].Ye[a] + sections[nid].Yp2[a],
                            -sections[nid].Ye[a],
                        ]
                    )
                    if sections[nid].active_termination.left:
                        Yii[nid] = Yii[nid] + sections[nid].active_termination.left.Yn
                        Ji = (
                            sections[nid].active_termination.left.Jn[a]
                            - sections[nid].Ie[a]
                        )
                elif nid == self.right_edge:
                    Ji = sections[-1].Ie[a]
                    Yii[nid - 1 : nid + 1] = np.array(
                        [
                            -sections[-1].Ye[a],
                            sections[-1].Yp2[a] + sections[-1].Ye[a],
                        ]
                    )
                    if sections[-1].active_termination.right:
                        Yii[nid] = Yii[nid] + sections[-1].active_termination.right.Yn
                        Ji = Ji - sections[-1].active_termination.right.Jn[a]
                else:
                    Ji = sections[nid - 1].Ie[a] - sections[nid].Ie[a]
                    Yii[nid - 1 : nid + 2] = np.array(
                        [
                            -sections[nid - 1].Ye[a],
                            sections[nid - 1].Ye[a]
                            + sections[nid].Ye[a]
                            + sections[nid - 1].Yp2[a]
                            + sections[nid].Yp2[a],
                            -sections[nid].Ye[a],
                        ]
                    )
                setattr(node, "Ji", Ji)
                setattr(node, "Yii", Yii)
                self.nodes[nid][a] = node
        return

    def solve_admitance_matrix(self):
        """
        Solve: [V] = inv([Y]).[J]
        """
        self.V = {}
        logger.info(f"Solving admitance matrix.")
        for a in self.components:
            logger.info(f"Solving for component {a}.")
            J, Y = [], []
            for nid in self.node_ids:
                n = self.nodes[nid][a]
                J.append(n.Ji)
                Y.append(n.Yii)
            J, Y = np.array(J), np.array(Y)
            logger.info(f"Sh(J):{J.shape}, Sh(Y):{Y.shape}")
            iY = np.linalg.inv(Y)
            self.V[a] = np.matmul(iY, J)
            logger.info(f"Sh(V):{self.V[a].shape}")
            logger.info(f"Set V[a] in each cable sections")
        for k, cs in enumerate(self.cable_sections):
            Vi, Vk = (
                self.V[self.components[0]][k, :],
                self.V[self.components[0]][k + 1, :],
            )
            if len(self.components) == 2:
                Vi += self.V[self.components[1]][k, :]
                Vk += self.V[self.components[1]][k + 1, :]
            setattr(cs.end_pot, "Vi", Vi)
            setattr(cs.end_pot, "Vk", Vk)
        return

    def consolidate_final_result(self):
        """
        Estimated Vi,k are the voltages at the end
        of the cable sections. Here we store estimated
        voltages in csv format, and line parameters
        (R, C, gma, L, Z0, Ln, Le, Ye[n,e], Yp2[n,e], Ie[n,e])
        in json format.
        """
        o = {"nodes": {}, "cables": {}}
        logger.info(f"Consolidate all results.")
        for bid, tx in enumerate(self.cable_sections):
            bid += 1
            o["cables"][bid] = {
                "R": tx.R,
                "C": tx.C,
                "gma": tx.gma,
                "Z0": tx.Z0,
                "ln": tx.length_north,
                "le": tx.length_east,
                "len_km": tx.length,
                "Ye": {},
                "Yp2": {},
                "Ie": {},
            }
            for a in self.components:
                o["cables"][bid]["Ye"][a] = tx.Ye[a]
                o["cables"][bid]["Yp2"][a] = tx.Yp2[a]
                o["cables"][bid]["Ie"][a] = tx.Ie[a].tolist()
        for nid in self.node_ids:
            nid = str(nid)
            for a in self.components:
                n = self.nodes[int(nid)][a]
                o["nodes"][nid] = {a: {}}
                o["nodes"][nid][a]["Ji"] = n.Ji.tolist()
                o["nodes"][nid][a]["Yii"] = n.Yii.tolist()
        self.result = o
        return

    def save(self, folder):
        """
        Save all analyzed data including
        Nodal Analysis
        """
        with open(folder + "est_cable_props.json", "w") as f:
            f.write(json.dumps(self.result, sort_keys=True, indent=4))

        self.tot_params.to_csv(folder + "sim-params.csv", float_format="%g")
        return

    def _pot_endCS_byComp_(self, cable_section_id, comp, unit="V", timestamp=None):
        """
        Provide the voltage at the ends of the
        cable to calculate total voltage by each components
        """
        unit = 1.0 if unit == "V" else 1000.0
        U0, U1 = (
            np.round(self.V[comp][cable_section_id, :] * u, 2),
            np.round(self.V[comp][cable_section_id + 1, :] * u, 2),
        )
        logger.info(
            f"Max(V) at the end of Section-{b}(Component-{comp}), {np.max(U0)} {np.max(U1)}"
        )
        if timestamp:
            U0, U1 = U0[timestamp], U1[timestamp]
        return U0, U1

    def _pot_endCS_(self, cable_section_id, unit="V", timestamp=None):
        """
        Provide the voltage at the ends of the
        cable to calculate total voltage
        """
        U0, U1 = self._pot_alongCS_byComp_(
            cable_section_id, self.components[0], unit, timestamp
        )
        if len(self.components) == 2:
            u0, u1 = self._pot_alongCS_byComp_(
                cable_section_id, self.components[1], unit, timestamp
            )
            U0 += u0
            U1 += u1
        return U0, U1

    def _pot_end_cable_byComp_(self, comp="X", unit="V", timestamp=None):
        """
        Provide the voltage at the ends of the
        cable to calculate total voltage by each components
        """
        u = 1.0 if unit == "V" else 1000.0
        U0, U1 = np.round(self.V[comp][0, :] * u, 2), np.round(
            self.V[comp][-1, :] * u, 2
        )
        logger.info(f"Max(V) at the end (Component-{comp}), {np.max(U0)} {np.max(U1)}")
        if timestamp:
            U0, U1 = U0[timestamp], U1[timestamp]
        return U0, U1

    def _pot_end_cable_(self, unit="V", timestamp=None):
        """
        Provide the voltage at the ends of the
        cable to calculate total voltage
        """
        U0, U1 = self._pot_end_cable_byComp_(self.components[0], unit, timestamp)
        if len(self.components) == 2:
            u0, u1 = self._pot_end_cable_byComp_(self.components[1], unit, timestamp)
            U0 += u0
            U1 += u1
        return U0, U1

    def _pot_along_cable_(self, timestamp, unit="V"):
        """ """
        Vcable, Lcable = [], []
        for cid, csec in enumerate(self.cable_sections):
            V, Lx = csec._pot_alongCS_()
            Vcable.extend(V.tolist())
            if cid == 0:
                Lcable = Lx.tolist()
            else:
                Lcable.extend((Lx + Lcable[-1]).tolist())
        return Vcable, Lcable
