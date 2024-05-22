"""Save the solved equilibrium data to an HDF5 file

"""

from enum import StrEnum
from re import A

import h5py
import numpy as np
import numpy.typing as npt

from .core import Machine, MatrixView, VectorView, IMatrixView, ModelType

# github copilot translation from C++ to Python


MU0 = 4.0e-7 * 3.14159265358979323846


class GROUP(StrEnum):
    GRID   = "/Grid"
    BOUND  = "/Boundaries"
    FLUX   = "/FluxFunctions"
    SCALAR = "/Scalars"


class DS_NAME(StrEnum):
    CUR_NAME = "Current"
    PSI_NAME = "Psi"
    RES_NAME = "Residuals"
    DIMX_NAME = "R"
    DIMZ_NAME = "Z"
    PSIX_NAME = "PsiNorm"
    MODB_NAME = "B2"
    BpX_NAME = "Bp_R"
    BpZ_NAME = "Bp_Z"
    TFLUX_NAME = "ToroidalFlux"
    PRESS_NAME = "Pressure"
    BETA_NAME = "Beta"
    LCFS_NAME = "LCFS"
    FCFS_NAME = "FCFS"
    PSI_1D = "psi"
    PRESS_1D = "ppsi"
    G_1D = "Gpsi"
    PP_1D = "pprime"
    G2P_1D = "G2prime"
    FFP_1D = "ffprime"
    Q_1D = "qpsi"
    V_1D = "Vprime"
    VOL_1D = "Vpsi"
    SHEAR_1D = "Shear"
    WELL_1D = "Well"
    B2_1D = "B2Ave"
    BETA_1D = "BetaAve"
    J_1D = "JAve"
    BETAMAX_1D = "BetaMax"
    X_BETAMAX_1D = "RBetaMax"
    Z_BETAMAX_1D = "ZBetaMax"
    B_BETAMAX_1D = "BBetaMax"
    BMAX_1D = "BMax"
    X_BMAX_1D = "RBMax"
    Z_BMAX_1D = "ZBMax"
    IP_0D = "Ip"
    BT_0D = "B0"
    R0_0D = "R0"
    Z0_0D = "Z0"
    R0Z0_0D = "R0Z0"
    PSIAXIS_0D = "PsiMagX"
    PSIFCFS_0D = "PsiFCFS"
    PSILCFS_0D = "PsiLCFS"
    RMAGX_0D = "RMagX"
    ZMAGX_0D = "ZMagX"
    OLIM_NAME = "olim"
    ILIM_NAME = "ilim"


ArrayLike = npt.ArrayLike | MatrixView | VectorView | IMatrixView


def _save_0D(loc: h5py.Group, name: str, units: str,
             data: float) -> h5py.Dataset:
    ds = loc.create_dataset(name, data=[data])
    ds.attrs["UNITS"] = units
    return ds


def _save_1D(loc: h5py.Group, name: str, units: str, dim: h5py.Dataset,
             data: ArrayLike):
    arr = np.array(data)
    ds = loc.create_dataset(name, data=arr)
    ds.attrs["UNITS"] = units
    ds.dims[0].attach_scale(dim)
    return ds


def _save_2D(loc: h5py.Group, name: str, units: str, dimr: h5py.Dataset,
             dimz: h5py.Dataset, data: ArrayLike):
    arr = np.array(data)
    ds = loc.create_dataset(name, data=arr)
    ds.attrs["UNITS"] = units
    ds.dims[0].attach_scale(dimr)
    ds.dims[1].attach_scale(dimz)
    return ds


def save_to_hdf5(m: Machine, filename: str = None):

    if filename is None:
        filename = m.Oname + ".h5"

    # create the file.. (truncae if exists)

    pg = m.PsiGrid
    pl = m.Plasma
    R = np.array(pg.R)
    Z = np.array(pg.Z)

    with h5py.File(filename, mode="w") as h5f:
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
        dimr = grid.create_dataset(DS_NAME.DIMX_NAME,
                                   data=R)
        dimr.attrs["UNITS"] = "m"
        dimr.make_scale('R')
        dimz = grid.create_dataset(DS_NAME.DIMZ_NAME,
                                   data=Z)
        dimz.attrs["UNITS"] = "m"
        dimz.make_scale('Z')

        # write the scalar data
        _save_0D(scal, DS_NAME.RMAGX_0D, "m", pg.RMagAxis)
        _save_0D(scal, DS_NAME.ZMAGX_0D, "m", pg.ZMagAxis)

        # write 2d psigrid data
        _save_2D(grid, DS_NAME.CUR_NAME, "A/m^2", dimr, dimz,
                 np.array(pg.Current)/MU0)
        _save_2D(grid, DS_NAME.PSI_NAME, "Wb", dimr, dimz, pg.Psi)
        _save_2D(grid, DS_NAME.RES_NAME, "Wb", dimr, dimz,
                 np.array(pg.Residual)/MU0)

        # write 0d plasma data
        _save_0D(scal, DS_NAME.R0_0D, "m", pl.R0)
        _save_0D(scal, DS_NAME.Z0_0D, "m", pl.Z0)
        _save_0D(scal, DS_NAME.BT_0D, "T", pl.B0)
        _save_0D(scal, DS_NAME.IP_0D, "A", pl.Ip)
        
        # write 2d plasma data
        _save_2D(grid, DS_NAME.MODB_NAME, "T^2", dimr, dimz, pl.B2)
        _save_2D(grid, DS_NAME.BpX_NAME, "T/m", dimr, dimz,
                 np.array(pl.GradPsiZ)/(2*np.pi*R)) 
        _save_2D(grid, DS_NAME.BpX_NAME, "T/m", dimr, dimz,
                 -np.array(pl.GradPsiR)/(2*np.pi*R))
        _save_2D(grid, DS_NAME.TFLUX_NAME, "Wb/R0B0", dimr, dimz, pl.G)
        
        match(p.ModelType):
            case (ModelType.Std
                | ModelType.DipoleIntStable
                | ModelType.DipoleStd):
                _save_2D(grid, DS_NAME.PRESS_NAME, "Pa", dimr, dimz, 
                         np.array(pl.Piso)/MU0)
                beta = 2*np.array(pl.Piso)/np.array(pl.B2)
                _save_2D(grid, DS_NAME.BETA_NAME, "", dimr, dimz, beta)
            case _:
                pass
        
        # write 1d flux functions
            
        
                
        
        
        
        
        
        
        