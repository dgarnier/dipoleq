# Tests for the IMAS interface
# for now just to include coverage
# should add more tests

import os
from pathlib import Path

import fixtures.typeguard_fix  # noqa: F401  # pylint: disable=unused-import
import pytest
from dipoleq import Machine
from dipoleq.imas import imas_input_params, load_imas_data_structure

data_dir = Path(os.path.realpath(__file__)).parent / "data"


def test_create_ods() -> None:  # pylint: disable=missing-function-docstring
    m = Machine.from_file(data_dir / "beta1.yaml")
    m.solve()
    ods = m.to_omas()
    assert ods["equilibrium.code.name"] == "DipolEQ"


@pytest.fixture(scope="session")
def test_imas_save(tmp_path_factory: pytest.TempPathFactory) -> Path:
    # check that one can solve and save to h5 from yaml input
    m1 = Machine.from_yaml(data_dir / "beta1.yaml")
    m1.solve()
    ods = m1.to_omas()
    fpath = tmp_path_factory.mktemp("data") / "beta1_imas.h5"
    ods.save(str(fpath))
    return fpath


def test_imas_params(test_imas_save: Path) -> None:
    """IMAS XML input parameter re-read test"""
    ods = load_imas_data_structure(test_imas_save)
    assert ods["equilibrium.code.name"] == "DipolEQ"
    input_data = imas_input_params(ods)
    assert input_data
    assert input_data["Oname"] == "beta1"
    m1 = Machine.from_dict(input_data)
    m0 = Machine.from_file(data_dir / "beta1.yaml")
    assert m0 == m1
    # when solving, they results should be the same
    # but somehow not.. why?
    # m1.solve()
    # m0.solve()
