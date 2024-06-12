#!/usr/bin/env python
# Convert a dipoleq h5 file to a g-eqdsk file
# Usage: python h5togeqdsk.py <dipoleq.h5> ...

import typing as tp
from collections.abc import Mapping
from os import PathLike
from pathlib import Path

import h5py
import numpy as np
from freeqdsk import geqdsk
from numpy.typing import ArrayLike, NDArray


def dipoleq_lim_to_eqdsk(lim: NDArray[np.float64]) -> NDArray[np.float64]:
    """lim is a 3D array with shape (nlim, 2, 2)
    with a start and end point for each limiter
    in this we will assume that the start is the end of the last"""
    newlim = np.zeros((lim.shape[0] + 1, 2))
    newlim[:-1] = lim[:, 0]
    newlim[-1] = lim[-1, 1]
    return newlim


# h5py is not typed, and the stubs is horrible
@tp.no_type_check
def dipoleq_h5f_to_freeqdsk(
    h5f: h5py.Group, COCOS: int = 3, NormalizeAtAxis: bool = True
) -> tuple[dict[str, int | float | NDArray[np.float64]], str]:
    """Extract geqdsk data from a dipoleq h5 file"""

    # future version of geqdsk will have type hints
    # from typing import Dict, Union
    # from numpy.typing import ArrayLike
    # gdata = Dict[str, Union[int, float, ArrayLike]]

    # HDF5 and DipolEq is COCOS=11 (and HDF5 is in COCOS=11 for ITER)
    # this converts to COCOS=1 (for EFIT and g-eqdsk)
    # EFIT is COCOS=3 according to COCOS paper

    gdata = {}

    scale_psi = 0.5 / np.pi if COCOS < 10 else 1.0
    if COCOS % 10 == 3:
        scale_psi = -scale_psi

    # commputational domain
    Grid = h5f["/Grid"]
    Flux = h5f["/FluxFunctions"]
    eq0d = h5f["/Scalars"]

    R = Grid["R"][()]
    Z = Grid["Z"][()]

    # 0D values
    gdata["rdim"] = max(R) - min(R)
    gdata["rleft"] = min(R)
    gdata["zdim"] = max(Z) - min(Z)
    gdata["zmid"] = (max(Z) + min(Z)) / 2
    # reference values
    gdata["rcentr"] = eq0d["R0"][()]
    gdata["zcentr"] = eq0d["Z0"][()]
    gdata["bcentr"] = eq0d["B0"][()]
    Fscale = eq0d["R0"][()] * eq0d["B0"][()]

    # plasma current
    gdata["cpasma"] = eq0d["Ip"][()]
    gdata["rmagx"] = eq0d["RMagX"][()]
    gdata["zmagx"] = eq0d["ZMagX"][()]
    # psi values
    PsiFCFS: float = eq0d["PsiFCFS"][()] * scale_psi
    PsiLCFS: float = eq0d["PsiLCFS"][()] * scale_psi
    PsiMagX: float = eq0d["PsiMagX"][()] * scale_psi

    # 1D values
    # geqdsk assumes that the radial resolution is the same
    # as the number of X gridpoints
    # also.. some code requires that psi normalized STARTS at
    # the magnetic axis

    psi1D = Flux["psi"][()] * scale_psi
    if NormalizeAtAxis:
        psi = np.linspace(PsiMagX, PsiLCFS, len(R))
    else:
        psi = np.linspace(PsiFCFS, PsiLCFS, len(R))

    def regrid(y: ArrayLike) -> NDArray[np.float64]:
        return np.interp(psi, psi1D, y)

    gdata["simagx"] = psi[0]
    gdata["sibdry"] = psi[-1]
    gdata["fpol"] = regrid(Flux["Gpsi"][()] * Fscale)
    gdata["pres"] = regrid(Flux["ppsi"][()])  # codespell:ignore pres
    gdata["ffprime"] = regrid(Flux["G2prime"][()] * Fscale / scale_psi)
    gdata["pprime"] = regrid(Flux["pprime"][()] / scale_psi)
    gdata["qpsi"] = regrid(Flux["qpsi"][()])

    # 2D values
    gdata["psi"] = Grid["Psi"][()].T * scale_psi
    # should be fortran order
    # should be in Wb/rad

    # Boundary values
    lcfs: ArrayLike = h5f["/Boundaries/LCFS"][()]
    fcfs: ArrayLike = h5f["/Boundaries/FCFS"][()]
    gdata["rbdry"] = lcfs[:, 0]
    gdata["zbdry"] = lcfs[:, 1]
    gdata["ribdry"] = fcfs[:, 0]
    gdata["zibdry"] = fcfs[:, 1]

    olimq = h5f["/Boundaries/olim"][()]
    ilimq = h5f["/Boundaries/ilim"][()]
    olim = dipoleq_lim_to_eqdsk(olimq)
    gdata["rlim"] = olim[:, 0]
    gdata["zlim"] = olim[:, 1]
    ilim = dipoleq_lim_to_eqdsk(ilimq)
    gdata["rlimi"] = ilim[:, 0]
    gdata["zlimi"] = ilim[:, 1]
    oname = h5f.attrs["ONAME"]
    oname = oname if isinstance(oname, str) else str(oname, "utf-8")

    return (gdata, oname)


@tp.no_type_check
def plot_h5eq(h5eq: h5py.Group) -> None:
    """Plot the equilibrium from a dipoleq h5 group"""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib is required for plotting.")
        return

    _, ax = plt.subplots()
    ax.contour(h5eq["Grid"]["R"], h5eq["Grid"]["Z"], h5eq["Grid"]["Psi"], 100)
    runits = h5eq["Grid"]["R"].attrs["UNITS"]
    zunits = h5eq["Grid"]["Z"].attrs["UNITS"]
    ax.set_xlabel(f"R [{runits if isinstance(runits, str) else runits.decode()}]")
    ax.set_ylabel(f"Z [{zunits if isinstance(zunits, str) else zunits.decode()}]")
    ax.set_aspect("equal")
    ax.plot(h5eq["Boundaries"]["LCFS"][:, 0], h5eq["Boundaries"]["LCFS"][:, 1], "b--")
    ax.plot(h5eq["Boundaries"]["FCFS"][:, 0], h5eq["Boundaries"]["FCFS"][:, 1], "b--")
    ilim = h5eq["Boundaries"]["ilim"]
    olim = h5eq["Boundaries"]["olim"]
    for lim in [ilim, olim]:
        for i in range(lim.shape[0]):
            ax.plot(lim[i, :, 0], lim[i, :, 1], "k-")
    plt.show()


def write_geqdsk(
    gdata: Mapping[str, int | float | ArrayLike],
    filename: str | PathLike[str],
    oname: str,
) -> None:
    """Write a g-eqdsk file from a dictionary of data"""
    path = Path(filename)
    with path.open("w") as fh:
        geqdsk.write(gdata, fh, label=f"DipEq:{oname}")


def write_fcfs_csv(
    h5path: Path, gdata: dict[str, int | float | NDArray[np.float64]]
) -> None:
    """Write the FCFS and FLIM to a CSV file"""
    stem = h5path.stem
    parent = h5path.parent
    with (parent / f"{stem}_fcfs.csv").open("w", encoding="utf-8") as fh:
        fcfs = np.column_stack((gdata["ribdry"], gdata["zibdry"]))
        np.savetxt(fh, fcfs, delimiter=",", header="r,z")

    with (parent / f"{stem}_flim.csv").open("w", encoding="utf-8") as fh:
        flim = np.column_stack((gdata["rlimi"], gdata["zlimi"]))
        np.savetxt(fh, flim, delimiter=",", header="r,z")


def h5togeqdsk(
    h5file: Path | PathLike[str],
    plot: bool = False,
    NormalizeAtAxis: bool = True,
    suffix: str = ".geqdsk",
) -> dict[str, int | float | NDArray[np.float64]]:
    """Save a dipoleq h5 file to a g-eqdsk file"""

    h5path = Path(h5file)
    with h5py.File(h5file) as h5f:  # type: ignore[arg-type]
        if plot:
            plot_h5eq(h5f)
        gdata, oname = dipoleq_h5f_to_freeqdsk(h5f, NormalizeAtAxis=NormalizeAtAxis)

    write_geqdsk(gdata, h5path.with_suffix(suffix), oname)
    write_fcfs_csv(h5path, gdata)
    return gdata  # type: ignore[no-any-return]


def main() -> None:
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Convert a dipoleq hdf5 file to a g-eqdsk file")
    parser.add_argument(
        "h5files", metavar="h5file", type=str, nargs="+", help="dipoleq hdf5 file(s)"
    )
    parser.add_argument(
        "--plot", "-p", action="store_true", default=False, help="Plot the g-eqdsk"
    )

    args = parser.parse_args()

    for h5file in args.h5files:
        h5togeqdsk(Path(h5file), plot=args.plot)


if __name__ == "__main__":
    main()
