"""Save the solved equilibrium data to an HDF5 file
this should be identical to what's in the c code
"""

from pathlib import Path

import numpy as np
import numpy.typing as npt
from h5py import Dataset, File, Group

from .core import IMatrixView, Machine, MatrixView, VectorView
from .core import ModelType as MT
from .ds_names import DS_NAME, GROUP

# github copilot translation from C++ to Python


MU0 = 4.0e-7 * 3.14159265358979323846


ArrayLike = npt.ArrayLike | MatrixView | VectorView | IMatrixView


def _save_0D(loc: Group, name: str, units: str, data: float) -> Dataset:
    ds = loc.create_dataset(name, data=data)  # type: ignore[arg-type]
    ds.attrs["UNITS"] = units
    return ds


def _save_1D(
    loc: Group, name: str, units: str, dim: Dataset, data: ArrayLike
) -> Dataset:
    arr = np.array(data)
    ds = loc.create_dataset(name, data=arr)
    ds.attrs["UNITS"] = units
    ds.dims[0].attach_scale(dim)  # type: ignore[attr-defined]
    return ds


def _save_2D(
    loc: Group, name: str, units: str, scl_r: Dataset, scl_z: Dataset, data: ArrayLike
) -> Dataset:
    arr = np.array(data)
    ds = loc.create_dataset(name, data=arr)
    ds.attrs["UNITS"] = units
    ds.dims[0].attach_scale(scl_r)  # type: ignore[attr-defined]
    ds.dims[1].attach_scale(scl_z)  # type: ignore[attr-defined]
    return ds


def save_flux_functions(flux: Group, m: Machine) -> None:
    # create flux dimensions
    pl = m.Plasma
    if pl.PsiX_pr:  # if this exists the rest should too
        dimp = flux.create_dataset(DS_NAME.DIMX_NAME, data=pl.PsiX_pr)
        dimp.attrs["UNITS"] = "1"
        dimp.make_scale("Normalized Magnetic Flux")  # type: ignore[attr-defined]

        _save_1D(flux, DS_NAME.PSI_1D, "Wb", dimp, pl.Psi_pr)
        _save_1D(flux, DS_NAME.PRESS_1D, "Pa", dimp, pl.P_pr)
        _save_1D(flux, DS_NAME.G_1D, "1", dimp, pl.G_pr)
        _save_1D(flux, DS_NAME.PP_1D, "Pa/Wb", dimp, pl.Pp_pr)
        _save_1D(flux, DS_NAME.G2P_1D, "1/Wb", dimp, pl.G2p_pr)
        _save_1D(flux, DS_NAME.Q_1D, "", dimp, pl.q_pr)
        _save_1D(flux, DS_NAME.V_1D, "m^3/Wb", dimp, pl.Volp_pr)
        _save_1D(flux, DS_NAME.VOL_1D, "m^3", dimp, pl.Vol_pr)
        _save_1D(flux, DS_NAME.SHEAR_1D, "", dimp, pl.S_pr)
        _save_1D(flux, DS_NAME.WELL_1D, "", dimp, pl.Well_pr)
        _save_1D(flux, DS_NAME.B2_1D, "T^2", dimp, pl.B2_pr)
        _save_1D(flux, DS_NAME.BETA_1D, "", dimp, pl.Beta_pr)
        _save_1D(flux, DS_NAME.J_1D, "A/m^2", dimp, pl.J_pr)
        _save_1D(flux, DS_NAME.BETAMAX_1D, "", dimp, pl.BetaMax_pr)
        _save_1D(flux, DS_NAME.R_BETAMAX_1D, "m", dimp, pl.RBetaMax_pr)
        _save_1D(flux, DS_NAME.Z_BETAMAX_1D, "m", dimp, pl.ZBetaMax_pr)
        _save_1D(flux, DS_NAME.B_BETAMAX_1D, "T", dimp, pl.BBetaMax_pr)
        _save_1D(flux, DS_NAME.BMAX_1D, "T", dimp, pl.BMax_pr)
        _save_1D(flux, DS_NAME.R_BMAX_1D, "m", dimp, pl.RBMax_pr)
        _save_1D(flux, DS_NAME.Z_BMAX_1D, "m", dimp, pl.ZBMax_pr)


def save_boundaries(loc: Group, m: Machine) -> None:
    """Create and save the LCFS and FCFS boundaries"""
    pg = m.PsiGrid
    lcfs_r, lcfs_z = pg.get_contour(1.0)
    if lcfs_r is not None:
        lcfs = np.vstack((lcfs_r, lcfs_z)).T
        lcfs_ds = loc.create_dataset(DS_NAME.LCFS_NAME, data=lcfs)
        lcfs_ds.attrs["UNITS"] = "m"
        lcfs_ds.attrs["DIMENSION"] = "cylindrical"
        lcfs_ds.attrs["FORMAT"] = "F7.4"

    if pg.PsiAxis == pg.PsiMagAxis:
        # no fcfs when the inner boundary is at the axis.
        return

    fcfs_r, fcfs_z = pg.get_contour(0.0)
    if fcfs_r is not None:
        fcfs = np.vstack((fcfs_r, fcfs_z)).T
        fcfs_ds = loc.create_dataset(DS_NAME.FCFS_NAME, data=fcfs)
        fcfs_ds.attrs["UNITS"] = "m"
        fcfs_ds.attrs["DIMENSION"] = "cylindrical"
        fcfs_ds.attrs["FORMAT"] = "F7.4"


def save_limiters(loc: Group, m: Machine) -> None:
    """Save the limiters"""

    olim = np.array(
        [[[lim.R1, lim.Z1], [lim.R2, lim.Z2]] for lim in m.Limiters if lim.Enabled > 0]
    )
    ilim = np.array(
        [[[lim.R1, lim.Z1], [lim.R2, lim.Z2]] for lim in m.Limiters if lim.Enabled < 0]
    )

    olim_ds = loc.create_dataset(DS_NAME.OLIM_NAME, data=olim)
    olim_ds.attrs["UNITS"] = "m"
    olim_ds.attrs["DIMENSION"] = "cylindrical"
    olim_ds.attrs["FORMAT"] = "F7.4"

    ilim_ds = loc.create_dataset(DS_NAME.ILIM_NAME, data=ilim)
    ilim_ds.attrs["UNITS"] = "m"
    ilim_ds.attrs["DIMENSION"] = "cylindrical"
    ilim_ds.attrs["FORMAT"] = "F7.4"


def save_to_hdf5(m: Machine, filename: str | Path | None = None) -> None:
    """Save the equilibrium data to an HDF5 file"""
    if filename is None:
        filename = m.Oname + ".h5"

    # create the file.. (truncate if exists)

    pg = m.PsiGrid
    pl = m.Plasma
    R = np.array(pg.R)
    Z = np.array(pg.Z)

    with File(filename, "w") as h5f:  # type: ignore[call-arg, arg-type]
        # put some info into the file
        h5f.attrs["TITLE"] = "Equilibrium data from dipoleq"
        h5f.attrs["VERSION"] = "0.1"
        h5f.attrs["ONAME"] = m.Oname.encode()

        # create the groups
        grid = h5f.create_group(GROUP.GRID)
        bound = h5f.create_group(GROUP.BOUND)
        flux = h5f.create_group(GROUP.FLUX)
        scal = h5f.create_group(GROUP.SCALAR)

        # create grid dimensions
        dimr = grid.create_dataset(DS_NAME.DIMX_NAME, data=R)
        dimr.attrs["UNITS"] = "m"
        dimr.make_scale("R")  # type: ignore[attr-defined]
        dimz = grid.create_dataset(DS_NAME.DIMZ_NAME, data=Z)
        dimz.attrs["UNITS"] = "m"
        dimz.make_scale("Z")  # type: ignore[attr-defined]

        # write the scalar grid data
        _save_0D(scal, DS_NAME.RMAGX_0D, "m", pg.RMagAxis)
        _save_0D(scal, DS_NAME.ZMAGX_0D, "m", pg.ZMagAxis)

        # write 2d psigrid data
        _save_2D(
            grid, DS_NAME.CUR_NAME, "A/m^2", dimr, dimz, np.array(pg.Current) / MU0
        )
        _save_2D(grid, DS_NAME.PSI_NAME, "Wb", dimr, dimz, pg.Psi)
        _save_2D(grid, DS_NAME.RES_NAME, "Wb", dimr, dimz, np.array(pg.Residual) / MU0)

        # write 0d plasma data
        _save_0D(scal, DS_NAME.R0_0D, "m", pl.R0)
        _save_0D(scal, DS_NAME.Z0_0D, "m", pl.Z0)
        _save_0D(scal, DS_NAME.R0Z0_0D, "m", pl.R0 * pl.Z0)
        _save_0D(scal, DS_NAME.BT_0D, "T", pl.B0)
        _save_0D(scal, DS_NAME.IP_0D, "A", pl.Ip)
        _save_0D(scal, DS_NAME.PSIAXIS_0D, "Wb", pl.PsiMagAxis)
        _save_0D(scal, DS_NAME.PSIFCFS_0D, "Wb", pl.PsiAxis)
        _save_0D(scal, DS_NAME.PSILCFS_0D, "Wb", pl.PsiLim)

        # write 2d plasma data
        _save_2D(grid, DS_NAME.MODB_NAME, "T^2", dimr, dimz, pl.B2)
        Br = np.array(pl.GradPsiZ) / (2 * np.pi * R)
        Bz = -np.array(pl.GradPsiR) / (2 * np.pi * R)
        _save_2D(grid, DS_NAME.BpX_NAME, "T/m", dimr, dimz, Br)
        _save_2D(grid, DS_NAME.BpZ_NAME, "T/m", dimr, dimz, Bz)
        _save_2D(grid, DS_NAME.TFLUX_NAME, "Wb/R0B0", dimr, dimz, pl.G)

        match pl.ModelType:
            case MT.Std | MT.DipoleIntStable | MT.DipoleStd | MT.DipoleStablePsiN:
                _save_2D(
                    grid, DS_NAME.PRESS_NAME, "Pa", dimr, dimz, np.array(pl.Piso) / MU0
                )
                beta = 2 * np.array(pl.Piso) / np.array(pl.B2)
                _save_2D(grid, DS_NAME.BETA_NAME, "", dimr, dimz, beta)
            case _:
                pass

        # write 1d flux functions

        save_flux_functions(flux, m)

        # do boundaries

        save_boundaries(bound, m)
        save_limiters(bound, m)
