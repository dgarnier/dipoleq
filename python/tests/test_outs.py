# Tests weird outputs...
# GS2 for now, but there are others
# for now just to include coverage
# should add more tests

import os
from pathlib import Path

from dipoleq import Machine

data_dir = Path(os.path.realpath(__file__)).parent / "data"


def test_gs2_out() -> None:  # pylint: disable=missing-function-docstring
    m = Machine.from_file(data_dir / "beta1.yaml")
    m.solve()
    m.write_GS2_geo()
