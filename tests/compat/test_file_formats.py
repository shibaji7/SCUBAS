# SPDX-License-Identifier: MIT
# Copyright (c) 2023 Shibaji Chakraborty
"""File formats SCUBAS reads and writes: IAGA-2002, CSV B-field input, and the two files ``Cable.save`` writes."""

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from loguru import logger

from scubas.cables import Cable, TransmissionLine
from scubas.datasets import PROFILES
from scubas.models import OceanModel

DATA = Path(__file__).parent / "data"


@pytest.fixture(autouse=True)
def _quiet():
    logger.disable("scubas")
    yield
    logger.enable("scubas")


@pytest.fixture
def model():
    return OceanModel(PROFILES.CS, flim=(1e-4, 1e-2))


def make_cable(components=("X", "Y"), sections=2):
    index = pd.date_range("2020-01-01", periods=4, freq="60s", name="Time")
    lines = []
    for i in range(sections):
        tl = TransmissionLine(
            f"S{i}", directed_length={"length_north": 100.0 + i, "length_east": 50.0}
        )
        ef = pd.DataFrame(
            {"X": [1.0, 2.0, 3.0, 4.0], "Y": [0.5, 0.25, 0.125, 0.0]}, index=index
        )
        tl.compute_eqv_pi_circuit(Efield=ef, components=list(components))
        lines.append(tl)
    return Cable(lines, components=list(components))


# ---- IAGA-2002 -------------------------------------------------------------------------------------------------------
def test_iaga_xyzf_file_is_read_with_a_time_index_and_gaps_as_nan(model):
    frame, header = model.read_iaga(DATA / "tst20200101vmin.min", return_header=True)
    assert list(frame.columns) == ["X", "Y", "Z", "F"]
    assert frame.index.name == "Time" and len(frame) == 24
    assert frame.index[1] - frame.index[0] == pd.Timedelta(minutes=1)
    assert frame["X"].isna().sum() == 2  # 99999.00 is the missing-value marker
    assert header["reported"] == "XYZF" and header["header_length"] == 12


def test_iaga_hdzf_is_converted_to_xy_with_the_documented_formula(model):
    """H and D (minutes of arc) become X = H cos(D/60 deg), Y = H sin(D/60 deg), and H, D are dropped."""
    frame = model.read_iaga(DATA / "tsthd20200101vmin.min")
    assert set(frame.columns) == {"Z", "F", "X", "Y"}
    i = 5
    h, d = 20000 + 50 * math.sin(i / 3), 120 + 4 * math.cos(i / 4)
    assert frame["X"].iloc[i] == pytest.approx(
        h * math.cos(math.radians(d / 60.0)), abs=0.01
    )
    assert frame["Y"].iloc[i] == pytest.approx(
        h * math.sin(math.radians(d / 60.0)), abs=0.01
    )


def test_iaga_hdzf_can_be_left_unconverted(model):
    frame = model.read_iaga(DATA / "tsthd20200101vmin.min", return_xyzf=False)
    assert list(frame.columns) == ["H", "D", "Z", "F"]


def test_iaga_without_reported_header_is_refused(model, tmp_path):
    bad = tmp_path / "bad.min"
    bad.write_text(
        " Format                 IAGA-2002                                   |\n"
    )
    with pytest.raises(ValueError, match="reported"):
        model.read_iaga(bad)


# ---- B-field aggregation -----------------------------------------------------------------------------------------------
def test_gaps_are_filled_by_second_order_spline_interpolation(model):
    b = model.read_Bfield_data([DATA / "tst20200101vmin.min"])
    assert not b.isna().any().any() and len(b) == 24
    csv = OceanModel(PROFILES.CS, flim=(1e-4, 1e-2)).read_Bfield_data(
        [DATA / "gap.csv"]
    )
    assert not csv.isna().any().any()
    assert 19970.0 < csv["X"].iloc[10] < 20000.0  # between its neighbours, not zero


def test_gaps_are_left_alone_when_the_check_is_off(model):
    b = model.read_Bfield_data([DATA / "gap.csv"], check_for_gaps=False)
    assert b["X"].isna().sum() == 1


def test_csv_needs_all_four_components(model, tmp_path):
    path = tmp_path / "xy.csv"
    path.write_text("Date,X,Y\n2020-01-01 00:00:00,1,2\n2020-01-01 00:01:00,1,2\n")
    with pytest.raises(ValueError, match="component"):
        model.read_Bfield_data([path])


def test_unknown_extension_is_refused(model, tmp_path):
    path = tmp_path / "b.dat"
    path.write_text("x")
    with pytest.raises(ValueError, match="Unsupported"):
        model.read_Bfield_data([path])


def test_files_are_concatenated_and_sorted_by_time(model, tmp_path):
    a = tmp_path / "a.csv"
    b = tmp_path / "b.csv"
    a.write_text(
        "Date,X,Y,Z,F\n2020-01-01 00:02:00,1,2,3,4\n2020-01-01 00:03:00,1,2,3,4\n"
    )
    b.write_text(
        "Date,X,Y,Z,F\n2020-01-01 00:00:00,5,6,7,8\n2020-01-01 00:01:00,5,6,7,8\n"
    )
    frame = model.read_Bfield_data([a, b])
    assert list(frame.index) == sorted(frame.index) and len(frame) == 4


def test_csv_date_column_name_is_configurable(model, tmp_path):
    path = tmp_path / "t.csv"
    path.write_text(
        "When,X,Y,Z,F\n2020-01-01 00:00:00,1,2,3,4\n2020-01-01 00:01:00,1,2,3,4\n"
    )
    frame = model.read_Bfield_data([path], csv_file_date_name="When")
    assert frame.index.name == "Time" and len(frame) == 2


# ---- what Cable.save writes ----------------------------------------------------------------------------------------------
def test_save_writes_two_files_with_the_documented_layout(tmp_path):
    cable = make_cable()
    out = tmp_path / "run"
    cable.save(out)
    assert sorted(p.name for p in out.iterdir()) == [
        "est_cable_props.json",
        "sim-params.csv",
    ]
    result = json.loads((out / "est_cable_props.json").read_text())
    assert sorted(result) == ["cables", "nodes"]
    assert sorted(result["cables"]) == ["1", "2"]
    assert set(result["cables"]["1"]) == {
        "C",
        "R",
        "Ie",
        "Yp2",
        "Ye",
        "Z0",
        "gma",
        "le",
        "ln",
        "len_km",
    }
    assert sorted(result["cables"]["1"]["Ye"]) == ["X", "Y"]
    assert sorted(result["nodes"]) == ["0", "1", "2"]
    assert sorted(result["nodes"]["0"]) == ["X", "Y"] and sorted(
        result["nodes"]["0"]["X"]
    ) == ["Ji", "Yii"]
    assert len(result["nodes"]["1"]["X"]["Yii"]) == 3  # one entry per node
    assert len(result["cables"]["1"]["Ie"]["X"]) == 4  # one per time step


def test_sim_params_csv_header_and_format(tmp_path):
    cable = make_cable()
    cable.save(tmp_path)
    header = (tmp_path / "sim-params.csv").read_text().splitlines()[0].split(",")
    assert header[0] == "Time"
    expected = [
        "V(v)",
        "E.X",
        "E.X.00",
        "E.X.01",
        "E.Y",
        "E.Y.00",
        "E.Y.01",
        "V(v).00",
        "V(v).01",
        "Vt(v)",
        "U0",
        "U1",
    ]
    assert header[1:] == expected
    frame = pd.read_csv(
        tmp_path / "sim-params.csv", index_col="Time", parse_dates=["Time"]
    )
    assert len(frame) == 4
    np.testing.assert_allclose(
        frame["Vt(v)"], cable.tot_params["Vt(v)"], rtol=1e-5
    )  # written with %g


def test_save_creates_missing_folders(tmp_path):
    make_cable(components=("X",), sections=1).save(tmp_path / "a" / "b")
    assert (tmp_path / "a" / "b" / "sim-params.csv").exists()


def test_cable_validates_its_arguments():
    with pytest.raises(ValueError, match="cable section"):
        Cable([], components=["X"])
    tl = TransmissionLine("S", directed_length={"length": 1.0})
    with pytest.raises(ValueError, match="component"):
        Cable([tl], components=[])


def test_width_must_be_positive():
    tl = TransmissionLine("S", directed_length={"length": 1.0})
    with pytest.raises(ValueError, match="width"):
        tl.calc_trasmission_line_parameters(width=0.0)


def test_frequency_limits_are_validated():
    with pytest.raises(ValueError, match="increasing"):
        OceanModel(PROFILES.CS, flim=(1.0, 0.1))
    with pytest.raises(ValueError, match="positive"):
        OceanModel(PROFILES.CS, flim=(0.0, 1.0))
    with pytest.raises(ValueError, match="two entries"):
        OceanModel(PROFILES.CS, flim=(1.0,))
