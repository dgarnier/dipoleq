import os
import warnings
from pathlib import Path

from .fixtures import typeguard_fix  # noqa: F401  # pylint: disable=unused-import  # isort: skip
import numpy as np
import pytest
from dipoleq import Machine

try:
    from dipoleq.imas import imas_input_params
    from imas import DBEntry
except ImportError:
    pass

data_dir = Path(os.path.realpath(__file__)).parent / "data"


@pytest.fixture(scope="session")
def test_imas_save(tmp_path_factory: pytest.TempPathFactory) -> Path:
    fpath = tmp_path_factory.mktemp("data") / "beta1_imas.nc"
    m = Machine.from_file(data_dir / "beta1.yaml")
    m.solve()
    warnings.filterwarnings("ignore")
    print(data_dir)
    with DBEntry(str(fpath), mode="w", dd_version="3.41.0+dipole") as db:
        m.to_imas(db)
    return fpath


@pytest.fixture(scope="session")
def test_imas_save_v4(tmp_path_factory: pytest.TempPathFactory) -> Path:
    fpath = tmp_path_factory.mktemp("data") / "beta1_imas.nc"
    m = Machine.from_file(data_dir / "beta1.yaml")
    m.solve()
    warnings.filterwarnings("ignore")
    print(data_dir)
    with DBEntry(str(fpath), mode="w", dd_version="4.0.0+dipole") as db:
        m.to_imas(db)
    return fpath


@pytest.mark.skipif(
    "DBEntry" not in globals(), reason="IMAS tests require IMAS-Python and dipmas"
)
def test_imas_params(test_imas_save: Path) -> None:
    with DBEntry(test_imas_save, "r", dd_version="3.41.0+dipole") as db:
        assert db.get("equilibrium")["code/name"] == "DipolEQ"

        input_data = imas_input_params(db.get("equilibrium"))
        assert input_data
        assert input_data["Oname"] == "beta1"
        m1 = Machine.from_dict(input_data)
        m0 = Machine.from_file(data_dir / "beta1.yaml")
        assert m0 == m1


@pytest.mark.skipif(
    "DBEntry" not in globals(), reason="IMAS tests require IMAS-Python and dipmas"
)
def test_imas_v4(test_imas_save: Path, test_imas_save_v4: Path) -> None:
    with DBEntry(test_imas_save, "r", dd_version="3.41.0+dipole") as db:
        psi_v3 = np.array(db.get("equilibrium")["time_slice[0]/profiles_1d/psi"])
    
    with DBEntry(test_imas_save_v4, "r", dd_version="4.0.0+dipole") as db:
        psi_v4 = np.array(db.get("equilibrium")["time_slice[0]/profiles_1d/psi"])
    
    np.testing.assert_allclose(psi_v3, -psi_v4)
