# Test dipoleq with basic LDX beta = 1 equilibrium
import os

# from os import PathLike
from pathlib import Path
from typing import Any

from .fixtures import typeguard_fix  # noqa: F401  # pylint: disable=unused-import  # isort: skip
import numpy as np
import pytest
from dipoleq import Machine
from dipoleq.core import Machine as coreMachine
from dipoleq.file_input import input_from_dotin
from dipoleq.solver import solve

data_dir = Path(os.path.realpath(__file__)).parent / "data"


def props(x: Any) -> dict[str, Any]:
    return {
        key: getattr(x, key)
        for key in dir(x)
        if not callable(getattr(x, key)) and not key.startswith("__")
    }


def test_FileInput() -> None:
    # see if we can read the file and get the right number of coils
    m = coreMachine(str(data_dir / "beta1.in"))
    assert m.NumCoils == 7
    i = 0
    for _ in m.Coils:
        i += 1
    assert i == 7


def test_solve_old() -> None:
    # see if we get the right current
    m = coreMachine(str(data_dir / "beta1.in"))
    solve(m)
    assert m.Plasma.Ip == pytest.approx(32984, rel=1e-4)


def test_solve_new() -> None:
    # see if we get the right current
    m = Machine.from_file(data_dir / "beta1.in")
    m.solve()
    assert m.Plasma.Ip == pytest.approx(32984, rel=1e-4)


def test_solve_yaml() -> None:
    # see if we get the right current
    m = Machine.from_file(data_dir / "beta1.yaml")
    m.solve()
    assert m.Plasma.Ip == pytest.approx(32984, rel=1e-4)
    assert m.is_diverted() is False


def test_read_dotin() -> None:
    # check that the input data is read
    mIn = input_from_dotin(data_dir / "beta1.in")
    assert mIn.NumCoils == 7


def test_old_new() -> None:
    # test that the C initialization and the Python initialization are the same
    m2 = coreMachine(str(data_dir / "beta1.in"))
    m1 = Machine.from_fileinput(data_dir / "beta1.in")
    m1.Iname = m2.Iname = "test"
    assert m1 == m2


def test_yaml_new() -> None:
    # test that the yaml and fileinput are the same
    m1 = Machine.from_yaml(data_dir / "beta1.yaml")
    m2 = Machine.from_fileinput(data_dir / "beta1.in")
    m1.Iname = m2.Iname = "test"
    m1.Name = m2.Name = "test"
    assert m1 == m2


@pytest.fixture(scope="session")
def test_yaml_save(tmp_path_factory: pytest.TempPathFactory) -> Path:
    # check that one can solve and save to h5 from yaml input
    m1 = Machine.from_yaml(data_dir / "beta1.yaml")
    m1.solve()
    fn = tmp_path_factory.mktemp("data") / "beta1.h5"
    m1.to_hdf5(fn)
    return fn


def test_h5togeqdsk(test_yaml_save: Path) -> None:
    # check that th old h5 reader can work
    from dipoleq.h5togeqdsk import h5togeqdsk

    gdata = h5togeqdsk(test_yaml_save)
    assert gdata["cpasma"] == pytest.approx(32984, rel=1e-4)


def test_plot_eq() -> None:
    m = Machine.from_file(data_dir / "beta1.yaml")
    m.solve()
    m.plot_eq()


def test_dipoleq_is_cocos_11() -> None:
    m = Machine.from_file(data_dir / "beta1.in")
    m.solve()

    # From Sauter et al, in "Computer Physics Communications"
    # Vol 184., p. 293, 2013.
    # Table I contains the expected values of various signs, which can be used to
    # distinguish between different COCOS versions.
    # It is difficult to check the difference between COCOS 11 and COCOS 12, because
    # almost all quantities are the same, so this test only checks that dipoleq
    # is COCOS 11 or 12.
    sign_delta_psi = np.sign(m.Plasma.PsiLCFS - m.Plasma.PsiFCFS)
    sign_q = np.sign(m.Plasma.q_pr)

    current_grid = np.array(m.PsiGrid.Current)
    pprime_grid = np.vectorize(m.Plasma.plasmaPp)(m.PsiGrid.Psi)

    # This value is the scale factor and sign in the Grad-Shafranov equation.
    # In COCOS terms, it should be equal to sigma_Bp * (2*pi)^e_Bp (which is +2pi for COCOS 11).
    with np.errstate(divide="ignore", invalid="ignore"):
        grad_shafranov_scale = -current_grid / (np.array(m.PsiGrid.R) * pprime_grid)

    assert sign_delta_psi == +1  # COCOS 1,11,2,12,5,15,6,16
    assert np.all(sign_q == +1)  # COCOS 1,11,2,12,7,17,8,18
    violations = (
        (~np.isclose(grad_shafranov_scale, +2 * np.pi, rtol=0.1))
        & np.isfinite(grad_shafranov_scale)
        & (grad_shafranov_scale != 0.0)
    )
    # Due to the residual not being precisely zero, there may be a few points
    # which are not equal to 2pi, but the vast majority (>99%) should be.
    assert (
        np.count_nonzero(violations) < 0.01 * grad_shafranov_scale.size
    )  # COCOS 11,12,15,16


if __name__ == "__main__":
    test_yaml_save(pytest.TempPathFactory())


# Expected UNITS attribute for every dataset written by saveh5.
# Keep this exhaustive: test_h5_units fails on any dataset not listed here,
# so a new dataset has to declare its units.
H5_UNITS = {
    "/Boundaries/FCFS": "m",
    "/Boundaries/LCFS": "m",
    "/Boundaries/ilim": "m",
    "/Boundaries/olim": "m",
    "/FluxFunctions/B2Ave": "T^2",
    "/FluxFunctions/BBetaMax": "T",
    "/FluxFunctions/BMax": "T",
    "/FluxFunctions/BetaAve": "",
    "/FluxFunctions/BetaMax": "",
    "/FluxFunctions/G2prime": "1/Wb",  # d(G^2)/dPsi, G dimensionless
    "/FluxFunctions/Gpsi": "1",  # G = F / (B0 R0)
    "/FluxFunctions/JAve": "A/m^2",
    "/FluxFunctions/PsiNorm": "1",
    "/FluxFunctions/RBMax": "m",
    "/FluxFunctions/RBetaMax": "m",
    "/FluxFunctions/Shear": "",  # d(ln q)/d(PsiNorm)
    "/FluxFunctions/Vprime": "m^3/Wb",
    "/FluxFunctions/Vpsi": "m^3",
    "/FluxFunctions/Well": "",  # differentiated wrt PsiNorm
    "/FluxFunctions/ZBMax": "m",
    "/FluxFunctions/ZBetaMax": "m",
    "/FluxFunctions/pprime": "Pa/Wb",
    "/FluxFunctions/ppsi": "Pa",
    "/FluxFunctions/psi": "Wb",
    "/FluxFunctions/qpsi": "",
    "/Grid/B2": "T^2",
    "/Grid/Beta": "",
    "/Grid/Bp_R": "T",  # Wb/m/m = T
    "/Grid/Bp_Z": "T",
    "/Grid/Current": "A/m^2",
    "/Grid/Pressure": "Pa",
    "/Grid/Psi": "Wb",
    "/Grid/R": "m",
    "/Grid/Residuals": "A/m^2",  # scaled by 1/mu0, same as Current
    "/Grid/ToroidalFlux": "1",  # G, dimensionless
    "/Grid/Z": "m",
    "/Scalars/B0": "T",
    "/Scalars/Ip": "A",
    "/Scalars/PsiFCFS": "Wb",
    "/Scalars/PsiLCFS": "Wb",
    "/Scalars/PsiMagX": "Wb",
    "/Scalars/R0": "m",
    "/Scalars/R0Z0": "T m",  # B0 * R0, the G -> F scale factor
    "/Scalars/RMagX": "m",
    "/Scalars/Z0": "m",
    "/Scalars/ZMagX": "m",
}


def test_h5_units(test_yaml_save: Path) -> None:
    # every dataset carries UNITS, and they are the expected ones
    import h5py

    found: dict[str, str] = {}

    def visit(name: str, obj: Any) -> None:
        if isinstance(obj, h5py.Dataset):
            units = obj.attrs["UNITS"]
            found["/" + name] = (
                units if isinstance(units, str) else units.decode("utf-8")
            )

    with h5py.File(test_yaml_save, "r") as h5f:
        h5f.visititems(visit)

    assert found == H5_UNITS


def test_h5_bpol_is_tesla(test_yaml_save: Path) -> None:
    # Bp_R and Bp_Z are fields, not gradients: Bp^2 + Bt^2 must equal B2 [T^2]
    import h5py

    with h5py.File(test_yaml_save, "r") as h5f:
        grid = h5f["/Grid"]
        R = grid["R"][()]
        br = grid["Bp_R"][()]
        bz = grid["Bp_Z"][()]
        bt = grid["ToroidalFlux"][()] * h5f["/Scalars/R0Z0"][()] / R
        b2 = grid["B2"][()]

    # skip the R -> 0 edge where Bt diverges
    sl = np.s_[5:-5, 5:-5]
    assert br[sl] ** 2 + bz[sl] ** 2 + bt[sl] ** 2 == pytest.approx(b2[sl], rel=1e-6)
