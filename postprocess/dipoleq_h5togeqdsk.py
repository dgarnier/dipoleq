#!/usr/bin/env python
# Convert a dipoleq h5 file to a g-eqdsk file
# Usage: python dipoleq_h5togeqdsk.py <dipoleq.h5> <g-eqdsk>

import re
import sys
import h5py
import numpy as np
from freeqdsk import geqdsk


def dipoleq_lim_to_eqdsk(lim):
    import numpy as np

    # lim is a 3D array with shape (nlim, 2, 2)
    # with a start and end point for each limiter
    # in this we will assume that the start is the end of the last
    newlim = np.zeros((lim.shape[0] + 1, 2))
    newlim[:-1] = lim[:, 0]
    newlim[-1] = lim[-1, 1]
    return newlim


def dipoleq_to_geqdsk(h5f):
    import numpy as np
    from typing import Dict, Union
    from numpy.typing import ArrayLike

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
        regrid = lambda y: y
    gdata["fpol"] = regrid(Flux["fpol"][()])
    gdata["pres"] = regrid(Flux["pres"][()])
    # gdata['ffprime']   = Flux['ffprime'][()]
    gdata["ffprime"] = regrid(Flux["dG2dPsi_1D"][()] / 2)
    gdata["pprime"] = regrid(Flux["pprime"][()])
    gdata["qpsi"] = regrid(Flux["qpsi"][()])

    # 2D values
    gdata["psi"] = Grid["Psi"][()]

    # Boundary values
    lcfs = h5f["/Boundaries/LCFS"][()]
    fcfs = h5f["/Boundaries/FCFS"][()]
    gdata["rbdry"] = lcfs[:, 0]
    gdata["zbdry"] = lcfs[:, 1]

    olimq = h5f["/Boundaries/lim"][()]
    ilimq = h5f["/Boundaries/ilim"][()]
    olim = dipoleq_lim_to_eqdsk(olimq)
    gdata["rlim"] = olim[:, 0]
    gdata["zlim"] = olim[:, 1]
    oname = str(h5f.attrs["ONAME"], "utf-8")

    return (gdata, oname)

if __name__ == "__main__":
    from argparse import ArgumentParser
    from pathlib import Path
    parser = ArgumentParser(description="Convert a dipoleq hdf5 file to a g-eqdsk file")
    parser.add_argument('h5files', metavar='h5file', type=str, nargs='+',
                        help="dipoleq hdf5 file(s)")
    args = parser.parse_args()

    for h5file in args.h5files:
        with h5py.File(h5file, "r") as h5f:
            gdata, oname = dipoleq_to_geqdsk(h5f)
            ofile = f"{Path(h5file).stem}.geqdsk"
        with open(ofile, "w") as fh:
            geqdsk.write(gdata, fh, label=f"DipQ:{oname}")
