from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from scubas.cables import Cable, CableSection, TransmissionLine
from scubas.datasets import PROFILES


def build_section_with_length(length: float = 2.0) -> CableSection:
    directed = {"length": length}
    section = CableSection("S1", directed_length=directed)
    return section


def test_cable_section_length_from_total():
    section = build_section_with_length(4.0)
    expected_component = 4.0 / np.sqrt(2.0)
    assert pytest.approx(section.length) == 4.0
    assert pytest.approx(section.length_north) == expected_component
    assert pytest.approx(section.length_east) == expected_component


def test_cable_section_length_from_edge_locations():
    section = CableSection(
        "S2",
        directed_length={
            "edge_locations": {
                "initial": {"lat": 0.0, "lon": 0.0},
                "final": {"lat": 1.0, "lon": 1.0},
            }
        },
    )
    assert section.length > 0.0
    assert section.length_north > 0.0
    assert section.length_east > 0.0


def test_cable_section_potential_profile_round_trip():
    section = build_section_with_length(1.0)
    section.gma = 0.01
    section.end_pot = SimpleNamespace(
        Vi=np.array([1.0, 0.8, 0.6]), Vk=np.array([0.5, 0.4, 0.3])
    )
    profile, distances = section._pot_alongCS_(idx=1, ln=4)
    assert profile.shape[0] == distances.shape[0]
    assert distances[-1] == pytest.approx(section.length)


def test_cable_section_potential_profile_requires_length():
    section = build_section_with_length(0.0)
    section.length = 0.0
    section.end_pot = SimpleNamespace(Vi=np.array([1.0]), Vk=np.array([0.5]))
    with pytest.raises(ValueError):
        section._pot_alongCS_()


def make_transmission_line(
    sec_id: str, with_termination: bool = False
) -> TransmissionLine:
    tl = TransmissionLine(sec_id, directed_length={"length": 1.0})
    if with_termination:
        tl.active_termination.left = SimpleNamespace(site=PROFILES.CS, width=1.0)
        tl.active_termination.right = None
    efield = pd.DataFrame(
        {"X": np.array([1.0, 0.5, 0.25])},
        index=pd.RangeIndex(3, name="Time"),
    )
    tl.compute_eqv_pi_circuit(Efield=efield, components=["X"])
    return tl


def test_transmission_line_compute_eqv_pi_circuit_populates_fields():
    tl = make_transmission_line("TL1")
    assert "X" in tl.Ye
    assert "X" in tl.Ie
    assert "Vj" in tl.V.columns
    assert tl.V.index.name == "Time"


def test_transmission_line_active_termination_populated():
    tl = make_transmission_line("TL2", with_termination=True)
    left = tl.active_termination.left
    assert hasattr(left, "Yn")
    assert "X" in left.Jn
    assert left.Z0 > 0


def test_cable_nodal_solution_generates_totals(monkeypatch):
    section_a = make_transmission_line("TL_A")
    section_b = make_transmission_line("TL_B")

    def fake_run(self):
        for section in self.cable_sections:
            length = len(section.Efield.index)
            section.end_pot = SimpleNamespace(Vi=np.zeros(length), Vk=np.zeros(length))

    def fake_solve(self):
        for section in self.cable_sections:
            length = len(section.Efield.index)
            section.end_pot = SimpleNamespace(Vi=np.zeros(length), Vk=np.zeros(length))

    def fake_consolidate(self):
        self.result = {"cables": {}, "nodes": {}}

    def fake_end(self, unit="V", timestamp=None):
        zeros = np.zeros(len(self.cable_sections[0].Efield.index))
        return zeros, zeros

    monkeypatch.setattr(Cable, "run_nodal_analysis", fake_run, raising=False)
    monkeypatch.setattr(Cable, "solve_admitance_matrix", fake_solve, raising=False)
    monkeypatch.setattr(
        Cable, "consolidate_final_result", fake_consolidate, raising=False
    )
    monkeypatch.setattr(Cable, "_pot_end_cable_", fake_end, raising=False)

    cable = Cable(cable_sections=[section_a, section_b], components=["X"])
    assert "V(v)" in cable.tot_params.columns
    assert "E.X.00" in cable.tot_params.columns
    assert isinstance(cable.result, dict)
    assert "cables" in cable.result
