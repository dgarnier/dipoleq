# Tests for the IMAS interface
# for now just to include coverage
# should add more tests

import os
from pathlib import Path

from .fixtures import typeguard_fix  # noqa: F401  # pylint: disable=unused-import
import pytest
from dipoleq import Machine
from dipoleq.omas import omas_input_params, load_omas_data_structure, ODS
import warnings

data_dir = Path(os.path.realpath(__file__)).parent / "data"


def test_create_ods() -> None:  # pylint: disable=missing-function-docstring
    m = Machine.from_file(data_dir / "beta1.yaml")
    m.solve()
    ods = m.to_omas()
    assert ods["equilibrium.code.name"] == "DipolEQ"


@pytest.fixture(scope="session")
def omas_h5_file(tmp_path_factory: pytest.TempPathFactory) -> Path:
    # check that one can solve and save to h5 from yaml input
    m1 = Machine.from_yaml(data_dir / "beta1.yaml")
    m1.solve()
    ods = m1.to_omas()
    fpath = tmp_path_factory.mktemp("data") / "beta1_omas.h5"
    ods.save(str(fpath))
    return fpath


def test_omas_params(omas_h5_file: Path) -> None:
    """IMAS XML input parameter re-read test"""
    ods = load_omas_data_structure(omas_h5_file)
    assert ods["equilibrium.code.name"] == "DipolEQ"
    input_data = omas_input_params(ods)
    assert input_data
    assert input_data["Oname"] == "beta1"
    m1 = Machine.from_dict(input_data)
    m0 = Machine.from_file(data_dir / "beta1.yaml")
    assert m0 == m1
    # when solving, they results should be the same
    # but somehow not.. why?
    # m1.solve()
    # m0.solve()


@pytest.fixture(scope="session")
def omas_nc_file(tmp_path_factory: pytest.TempPathFactory) -> Path:
    # check that one can solve and save to h5 from yaml input
    dir_path = tmp_path_factory.mktemp("data")
    return make_omas_nc_file(dir_path)


def make_omas_nc_file(dir_path: Path) -> Path:
    m1 = Machine.from_yaml(data_dir / "beta1.yaml")
    m1.solve()
    ods = m1.to_omas()
    warnings.filterwarnings('ignore')
    file_path = dir_path / "beta1_omas.nc"
    ods.save(str(file_path))
    return file_path


def test_omas_nc(omas_nc_file: Path) -> None:
    ods = load_omas_data_structure(omas_nc_file)
    reference_ods = load_omas_data_structure(data_dir / "reference/beta1_omas.nc")

    diff = ods.diff(reference_ods)
    assert not diff, diff
