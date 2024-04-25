#!/usr/bin/env python
# Convert a dipoleq h5 file to a g-eqdsk file
# Usage: python dipoleq_h5togeqdsk.py <dipoleq.h5> <g-eqdsk>

from typing import Dict, Union, Any
import h5py
import numpy as np
from freeqdsk import geqdsk


def dipoleq_lim_to_eqdsk(lim):
    '''lim is a 3D array with shape (nlim, 2, 2)
    with a start and end point for each limiter
    in this we will assume that the start is the end of the last'''
    newlim = np.zeros((lim.shape[0] + 1, 2))
    newlim[:-1] = lim[:, 0]
    newlim[-1] = lim[-1, 1]
    return newlim


def dipoleq_to_geqdsk(h5f):
    '''Extract geqdsk data from a dipoleq h5 file'''

    # future version of geqdsk will have type hints
    # from typing import Dict, Union
    # from numpy.typing import ArrayLike
    # gdata = Dict[str, Union[int, float, ArrayLike]]
    gdata = {}

    # 0D values
    # commputational domain
    Grid = h5f["/Grid"]
    Flux = h5f["/FluxFunctions"]
    eq0d = h5f["/Scalars"]

    R = Grid["R"][()]
    Z = Grid["Z"][()]
    gdata["rdim"] = max(R) - min(R)
    gdata["rleft"] = min(R)
    gdata["zdim"] = max(Z) - min(Z)
    gdata["zmid"] = (max(Z) + min(Z)) / 2
    # reference values
    gdata["rcentr"] = eq0d["rcentr"][()]
    gdata["bcentr"] = eq0d["bcentr"][()]
    # plasma current
    gdata["cpasma"] = eq0d["cpasma"][()]
    gdata["rmagx"] = eq0d["rmagx"][()]
    gdata["zmagx"] = eq0d["zmagx"][()]
    # psi values
    gdata["simagx"] = h5f["/Grid/Psi"].attrs["PsiAxis"][0]
    gdata["sibdry"] = h5f["/Grid/Psi"].attrs["PsiLim"][0]

    # 1D values
    # geqdsk assumes that the radial resolution is the same
    # as the number of X gridpoints
    psiX = Flux["PsiX"][()]
    if len(psiX) != len(R):
        # must interpolate
        psiXn = np.linspace(psiX[0], psiX[-1], len(R))

        def regrid(y):
            return np.interp(psiXn, psiX, y)

    else:
        def regrid(y):
            return y
    gdata["fpol"] = regrid(Flux["fpol"][()])
    gdata["pres"] = regrid(Flux["pres"][()])
    # gdata['ffprime']   = Flux['ffprime'][()]
    gdata["ffprime"] = regrid(Flux["dG2dPsi_1D"][()] / 2)
    gdata["pprime"] = regrid(Flux["pprime"][()])
    gdata["qpsi"] = regrid(Flux["qpsi"][()])

    # 2D values
    gdata["psi"] = Grid["Psi"][()].T # should be fortran order

    # Boundary values
    lcfs = h5f["/Boundaries/LCFS"][()]
    fcfs = h5f["/Boundaries/FCFS"][()]
    gdata["rbdry"] = lcfs[:, 0]
    gdata["zbdry"] = lcfs[:, 1]
    gdata["ribdry"] = fcfs[:, 0]
    gdata["zibdry"] = fcfs[:, 1]

    olimq = h5f["/Boundaries/lim"][()]
    ilimq = h5f["/Boundaries/ilim"][()]
    olim = dipoleq_lim_to_eqdsk(olimq)
    gdata["rlim"] = olim[:, 0]
    gdata["zlim"] = olim[:, 1]
    ilim = dipoleq_lim_to_eqdsk(ilimq)
    gdata["rlimi"] = ilim[:, 0]
    gdata["zlimi"] = ilim[:, 1]
    oname = str(h5f.attrs["ONAME"], "utf-8")

    return (gdata, oname)


def plot_h5eq(h5eq):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.contour(h5eq["Grid"]["R"], h5eq["Grid"]["Z"], h5eq["Grid"]["Psi"], 100)
    ax.set_xlabel(f'R [{h5eq["Grid"]["R"].attrs["UNITS"].decode()}]')
    ax.set_ylabel(f'Z [{h5eq["Grid"]["Z"].attrs["UNITS"].decode()}]')
    ax.set_aspect("equal")
    ax.plot(h5eq["Boundaries"]["LCFS"][:, 0], h5eq["Boundaries"]["LCFS"][:, 1], "b--")
    ax.plot(h5eq["Boundaries"]["FCFS"][:, 0], h5eq["Boundaries"]["FCFS"][:, 1], "b--")
    ilim = h5eq["Boundaries"]["ilim"]
    olim = h5eq["Boundaries"]["lim"]
    for lim in [ilim, olim]:
        for i in range(lim.shape[0]):
            ax.plot(lim[i, :, 0], lim[i, :, 1], "k-")
    plt.show()


if __name__ == "__main__":
    from argparse import ArgumentParser
    from pathlib import Path
    parser = ArgumentParser(description="Convert a dipoleq hdf5 file to a g-eqdsk file")
    parser.add_argument('h5files', metavar='h5file', type=str, nargs='+',
                        help="dipoleq hdf5 file(s)")
    parser.add_argument('--plot', '-p', action='store_true', 
                        help="Plot the g-eqdsk")
    
    args = parser.parse_args()

    for h5file in args.h5files:
        stem = Path(h5file).stem
        with h5py.File(h5file) as h5f:
            if args.plot:
                plot_h5eq(h5f)
            gdata, oname = dipoleq_to_geqdsk(h5f)

        with open(f"{stem}.geqdsk", "w") as fh:
            geqdsk.write(gdata, fh, label=f"DipQ:{oname}")

        with open(f"{stem}_fcfs.csv", "w", encoding='utf-8') as fh:
            fcfs = np.column_stack((gdata["ribdry"], gdata["zibdry"]))
            np.savetxt(fh, fcfs, delimiter=',', header='r,z', )

        with open(f"{stem}_flim.csv", "w", encoding='utf-8') as fh:
            flim = np.column_stack((gdata["rlimi"], gdata["zlimi"]))
            np.savetxt(fh, flim, delimiter=',', header='r,z', )
