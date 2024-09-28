"""Test plasma model = 7"""

import os
from pathlib import Path

import fixtures.typeguard_fix  # noqa: F401  # pylint: disable=unused-import
from dipoleq import Machine
from dipoleq.file_input import input_from_dotin
from dipoleq.input import CDipoleIntStableIn, CDipoleStablePsiNIn

data_dir = Path(os.path.realpath(__file__)).parent / "data"


def test_psinpeak_in() -> None:
    mIn = input_from_dotin(data_dir / "beta1_psinpeak.in")
    m = Machine.from_input_data(mIn)
    m.solve()


def test_psinpeak_py() -> None:
    mIn = input_from_dotin(data_dir / "beta1.in")
    pmOld = mIn.Plasma.Model
    assert isinstance(pmOld, CDipoleIntStableIn)
    mPM = CDipoleStablePsiNIn(
        Type=7,
        PsiNPeak=0.07,
        PEdge=pmOld.PEdge,
        PsiFlat=pmOld.PsiFlat,
        NSurf=pmOld.NSurf,
        fCrit=pmOld.fCrit,
    )
    mIn.Plasma.Model = mPM
    m = Machine.from_input_data(mIn)
    m.solve()
    m.plot_eq(show_peak=True)
