import os
from pathlib import Path

import fixtures.typeguard_fix  # noqa: F401  # pylint: disable=unused-import
import pytest
from dipoleq import Machine
from imas import DBEntry
import warnings

data_dir = Path(os.path.realpath(__file__)).parent / "data"


@pytest.fixture(scope="session")
def test_imas_save(tmp_path_factory: pytest.TempPathFactory) -> Path:
    fpath = tmp_path_factory.mktemp("data") / "beta1_imas.nc"
    m = Machine.from_file(data_dir / "beta1.yaml")
    m.solve()
    warnings.filterwarnings('ignore')
    print(data_dir)
    with DBEntry(str(fpath), mode='w', dd_version='3.41.0') as db:
        m.to_imas(db)
    return fpath


def test_imas_params(test_imas_save: Path) -> None:
    pass
