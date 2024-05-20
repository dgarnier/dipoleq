"""Save the solved equilibrium data to an HDF5 file

"""
    
from pyparsing import C
from . import MACHINE
import hdf5
from enum import Enum

# github copilot translation from C++ to Python

class Groups(Enum):
    GRID_GROUP = "/Grid"
    BOUND_GROUP = "/Boundaries"
    FLUX_GROUP = "/FluxFunctions"
    SCALAR_GROUP = "/Scalars"

class DataSets(Enum):
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

definitions = {
    "Groups": Groups,
    "DataSets": DataSets
}



def save_to_hdf5(m: MACHINE, filename: str = None)
    
    if filename is None:
        filename = m.Oname + '.h5'
    
    # create the file.. (truncae if exists)
    
    with h5py.File(filename, 'w') as h5f:
        # put some info into the file
        h5f.attrs['Title'] = "Equilibrium data from dipoleq"
        h5f.attts['Version'] = "0.1"
        
        # create the groups
        grid = h5f.create_group(Groups.GRID_GROUP.value)
        bound = h5f.create_group(Groups.BOUND_GROUP.value)
        flux = h5f.create_group(Groups.FLUX_GROUP.value)
        scal = h5f.create_group(Groups.SCALAR_GROUP.value)
        
        
        
        
        
    
    
    
    
    