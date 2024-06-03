"""ds_names.py
The names are used to access the datasets in the HDF5 file.
They can also be used in NetCDF files /etc.
This might, likely will be changed to use IMAS definitions in the future.
"""

try:
    from enum import StrEnum
except ImportError:  # python < 3.11
    from enum import Enum

    class StrEnum(str, Enum):  # type: ignore[no-redef]
        pass


class GROUP(StrEnum):
    GRID = "/Grid"
    BOUND = "/Boundaries"
    FLUX = "/FluxFunctions"
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
    J_1D = "JAve"  # codespell:ignore
    BETAMAX_1D = "BetaMax"
    R_BETAMAX_1D = "RBetaMax"
    Z_BETAMAX_1D = "ZBetaMax"
    B_BETAMAX_1D = "BBetaMax"
    BMAX_1D = "BMax"
    R_BMAX_1D = "RBMax"
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
