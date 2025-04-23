# Tests weird outputs...
# GS2 for now, but there are others
# for now just to include coverage
# should add more tests

import os
from pathlib import Path
from tempfile import TemporaryDirectory

from .fixtures import typeguard_fix  # noqa: F401  # pylint: disable=unused-import  # isort: skip
from dipoleq import Machine

data_dir = Path(os.path.realpath(__file__)).parent / "data"


def test_gs2_out() -> None:  # pylint: disable=missing-function-docstring
    m = Machine.from_file(data_dir / "beta1.yaml")
    m.solve()
    with TemporaryDirectory() as tmpdir:
        curdir = Path.cwd()
        os.chdir(tmpdir)
        m.write_GS2_geo()
        os.chdir(curdir)
        tmpfile = Path(tmpdir) / f"{m.Oname}_gs2.out"
        assert Path(tmpfile).exists()


def test_geqdsk_out() -> None:  # pylint: disable=missing-function-docstring
    m = Machine.from_file(data_dir / "beta1.yaml")
    m.solve()
    with TemporaryDirectory() as tmpdir:
        tmpfile = tmpdir + "/test.geqdsk"
        m.to_geqdsk(tmpfile)
        assert Path(tmpfile).exists()
