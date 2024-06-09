# Create density and temperature profiles
# Author: Darren Garnier, OpenStar Technologies.


import h5py
import numpy as np
import scipy as sp


def dipoleq_profiles(h5eq, npeak=1e20, eta=0.67, NormalizeAtAxis=True):
    FF = h5eq["FluxFunctions"]
    eq0d = h5eq["Scalars"]
    # L = FF["RBetaMax"] / eq0d["RMagX"][()]
    p = FF["ppsi"][:]
    p[0] = p[-1]  # don't have zero edge pressure
    psi = FF["psi"][:]
    # d2Vdpsi = np.gradient(FF["Vprime"][:], FF["psi"][:])
    dlnp = FF["pprime"][:] / p
    dlnp[0] = dlnp[1]
    # dlnV = d2Vdpsi / FF["Vprime"][:]
    # d = -dlnp / dlnV

    # print(psi.shape, dlnp.shape, d.shape)
    peaks, _ = sp.signal.find_peaks(p)
    # print(psi[peaks])
    peak = peaks[0]

    # d = dlnp/dlnV
    # eta = dlnT/dlnn
    # dlnp = dlnT + dlnn
    # dlnn = dlnp - dlnT = dlnp - eta * dlnn
    # dlnn = dlnp / (1 + eta)
    # dlnT = eta * dlnn = eta * dlnp / (1 + eta)
    eV = 1.602e-19
    dlnn = dlnp / (1 + eta)
    dlnT = eta * dlnn
    n = np.exp(sp.integrate.cumulative_trapezoid(dlnn, psi, initial=0))
    T = np.exp(sp.integrate.cumulative_trapezoid(dlnT, psi, initial=0))
    n = n * npeak / n[peak]
    T = T * p[peak] / (npeak * T[peak]) / eV  # in eV

    PsiFCFS = eq0d["PsiFCFS"][()]
    PsiLCFS = eq0d["PsiLCFS"][()]
    PsiMagX = eq0d["PsiMagX"][()]
    R = h5eq["Grid"]["R"][()]

    # 1D values
    # geqdsk assumes that the radial resolution is the same
    # as the number of X gridpoints
    # also.. some code requires that psi normalized STARTS at
    # the magnetic axis

    if NormalizeAtAxis:
        psin = np.linspace(PsiMagX, PsiLCFS, len(R))
    else:
        psin = np.linspace(PsiFCFS, PsiLCFS, len(R))

    def regrid(y):
        return np.interp(psin, psi, y)

    n = regrid(n)
    T = regrid(T)
    return psin, n, T


if __name__ == "__main__":
    from argparse import ArgumentParser
    from pathlib import Path

    parser = ArgumentParser(
        description="Create kinetic profiles from a dipoleq hdf5 file"
    )
    parser.add_argument(
        "h5files", metavar="h5file", type=str, nargs="+", help="dipoleq hdf5 file(s)"
    )
    parser.add_argument("--plot", "-p", action="store_true", help="Plot the profiles")
    parser.add_argument("--npeak", "-n", type=float, default=1e20, help="Peak density")
    parser.add_argument("--eta", "-e", type=float, default=0.67, help="eta = dlnT/dlnn")
    parser.add_argument(
        "--magX", "-x", action="store_true", default=True, help="Start at magnetic axis"
    )

    args = parser.parse_args()

    for h5file in args.h5files:
        pdir = Path(h5file).parent
        stem = Path(h5file).stem
        with h5py.File(h5file) as h5f:
            oname = str(h5f.attrs["ONAME"], "utf-8")
            psi, n, T = dipoleq_profiles(
                h5f, npeak=args.npeak, eta=args.eta, NormalizeAtAxis=args.magX
            )
        if args.plot:
            # plot_profiles(h5f)
            pass

        with (pdir / f"{stem}_profile.csv").open("w", encoding="utf-8") as fh:
            data = np.column_stack(psi, n, T)
            np.savetxt(
                fh,
                data,
                delimiter=",",
                header="psi, n, T",
            )
