"""
Python bindings for DipolEq
"""
from __future__ import annotations
import numpy
import typing
__all__ = ['CPlasmaModel', 'CircleType', 'Coil', 'Coils', 'IMatrixView', 'Limiter', 'Limiters', 'Machine', 'MatrixView', 'MeasType', 'Measure', 'Measures', 'ModelType', 'Plasma', 'PsiGrid', 'Separatricies', 'Separatrix', 'Shell', 'Shells', 'SubCoil', 'SubCoils', 'SubShell', 'SubShells', 'VectorView']
class CPlasmaModel:
    def model_input(self, arg0: str, arg1: str, arg2: str) -> None:
        """
        Input model parameters
        """
    def update_model(self, arg0: Machine) -> None:
        """
        Update the plasma model
        """
class CircleType:
    """
    Members:
    
      btcos
    
      brsin
    
      brcos
    """
    __members__: typing.ClassVar[dict[str, CircleType]]  # value = {'btcos': <CircleType.btcos: 1>, 'brsin': <CircleType.brsin: 2>, 'brcos': <CircleType.brcos: 3>}
    brcos: typing.ClassVar[CircleType]  # value = <CircleType.brcos: 3>
    brsin: typing.ClassVar[CircleType]  # value = <CircleType.brsin: 2>
    btcos: typing.ClassVar[CircleType]  # value = <CircleType.btcos: 1>
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class Coil:
    CoilCurrent: float
    Enabled: int
    def __init__(self, arg0: int) -> None:
        """
        Create Coil
        """
    def compute_SubCoils(self, arg0: PsiGrid) -> None:
        """
        Compute subcoils
        """
    @property
    def Name(self) -> str:
        """
        Name of the coil
        """
    @Name.setter
    def Name(self, arg1: str) -> None:
        ...
    @property
    def NumSubCoils(self) -> int:
        """
        Number of subcoils, setting will erase old subcoils
        """
    @NumSubCoils.setter
    def NumSubCoils(self, arg1: int) -> None:
        ...
    @property
    def R(self) -> float:
        """
        Coil centroid R
        """
    @R.setter
    def R(self, arg0: float) -> None:
        ...
    @property
    def SubCoils(self) -> SubCoils:
        """
        Return vector of SubCoils
        """
    @property
    def Z(self) -> float:
        """
        Coil centroid Z
        """
    @Z.setter
    def Z(self, arg0: float) -> None:
        ...
    @property
    def dR(self) -> float:
        """
        Coil radial width
        """
    @dR.setter
    def dR(self, arg0: float) -> None:
        ...
    @property
    def dZ(self) -> float:
        """
        Coil vertical width
        """
    @dZ.setter
    def dZ(self, arg0: float) -> None:
        ...
class Coils:
    def __getitem__(self, arg0: int) -> Coil:
        ...
    def __len__(self) -> int:
        ...
    def __setitem__(self, arg0: int, arg1: Coil) -> None:
        ...
    def new_Coil(self, arg0: int) -> Coil:
        """
        Add a new coil
        """
class IMatrixView:
    pass
class Limiter:
    Enabled: int
    PsiMin: float
    R1: float
    R2: float
    Rmin: float
    Z1: float
    Z2: float
    Zmin: float
    def __init__(self) -> None:
        """
        Create Limiter
        """
    @property
    def Name(self) -> str:
        """
        Name of the limiter
        """
    @Name.setter
    def Name(self, arg1: str) -> None:
        ...
class Limiters:
    def __getitem__(self, arg0: int) -> Limiter:
        ...
    def __len__(self) -> int:
        ...
    def __setitem__(self, arg0: int, arg1: Limiter) -> None:
        ...
    def new_limiter(self) -> Limiter:
        """
        Add a new limiter
        """
class Machine:
    Confidence: float
    IterFixed: int
    IterFree: int
    LHGreenStatus: int
    MGreenStatus: int
    MaxIterFixed: int
    MaxIterFree: int
    MaxIterMCarlo: int
    NumEqualEq: int
    NumMCarloData: int
    NumMCarloEq: int
    RestartStatus: int
    RestartUnkns: int
    SGreenStatus: int
    SInductStatus: int
    VacuumOnly: int
    def __del__(self) -> None:
        """
        Free Machine and what is below it
        """
    @typing.overload
    def __init__(self) -> None:
        """
        Return new Machine struct
        """
    @typing.overload
    def __init__(self, arg0: str) -> None:
        """
        Initialize Machine with a file
        """
    def add_coil_J(self) -> None:
        """
        Add Coil currents
        """
    def add_shell_J(self) -> None:
        """
        Add Shell currents
        """
    def find_J(self) -> None:
        """
        Find J
        """
    def find_plasma_boundary(self) -> None:
        """
        Find plasma boundary
        """
    @typing.overload
    def find_shell_current(self) -> None:
        """
        Find shell current
        """
    @typing.overload
    def find_shell_current(self) -> None:
        """
        Find shell current
        """
    def free_bndry_greens(self) -> None:
        """
        Load boundary greens
        """
    def free_meas_greens(self) -> None:
        """
        Free measurement greens
        """
    def get_plasma_parameters(self) -> None:
        """
        Get plasma parameters
        """
    def init(self) -> None:
        """
        Initialize Machine
        """
    def least_squares(self, arg0: int) -> None:
        """
        Least squares
        """
    def load_bndry_greens(self) -> None:
        """
        Load boundary greens
        """
    def load_meas_greens(self) -> None:
        """
        Load measurement greens
        """
    def psi_boundary(self) -> None:
        """
        Calculate psi boundary
        """
    def read_restart(self) -> None:
        """
        Read the restart file
        """
    def set_NumSubShells(self, arg0: int, arg1: int) -> None:
        """
        Set the number of subshells of shell i to n
        """
    def set_NumSubcoils(self, arg0: int, arg1: int) -> None:
        """
        Set the number of subcoils for a coil
        """
    def set_coil_NumSubCoils(self, arg0: int, arg1: int) -> None:
        """
        Set the number of subcoils for a coil
        """
    def set_start_time(self) -> None:
        """
        Set the start time
        """
    def set_stop_time(self) -> None:
        """
        Set the end time
        """
    def write_restart(self) -> None:
        """
        Write the restart file
        """
    def zero_J(self) -> None:
        """
        Zero J
        """
    @property
    def Coils(self) -> Coils:
        """
        Array of COILS
        """
    @property
    def Iname(self) -> str:
        """
        Input name
        """
    @Iname.setter
    def Iname(self, arg1: str) -> None:
        ...
    @property
    def Info(self) -> str:
        """
        Info about the machine
        """
    @Info.setter
    def Info(self, arg1: str) -> None:
        ...
    @property
    def LHname(self) -> str:
        """
        LH name
        """
    @LHname.setter
    def LHname(self, arg1: str) -> None:
        ...
    @property
    def Limiters(self) -> Limiters:
        """
        Get the limiters
        """
    @property
    def MGname(self) -> str:
        """
        MG name
        """
    @MGname.setter
    def MGname(self, arg1: str) -> None:
        ...
    @property
    def Measures(self) -> Measures:
        """
        Get the measurements
        """
    @property
    def Name(self) -> str:
        """
        Set the name
        """
    @Name.setter
    def Name(self, arg1: str) -> None:
        ...
    @property
    def NumCoils(self) -> int:
        """
        Number of coils, setting will erase old coils
        """
    @NumCoils.setter
    def NumCoils(self, arg1: int) -> None:
        ...
    @property
    def NumLimiters(self) -> int:
        """
        Number of limiters, setting will erase old limiters
        """
    @NumLimiters.setter
    def NumLimiters(self, arg1: int) -> None:
        ...
    @property
    def NumMeasures(self) -> int:
        """
        Number of measurements, setting will erase old measurements
        """
    @NumMeasures.setter
    def NumMeasures(self, arg1: int) -> None:
        ...
    @property
    def NumSeps(self) -> int:
        """
        Number of separatrixes, setting will erase old separatrixes
        """
    @NumSeps.setter
    def NumSeps(self, arg1: int) -> None:
        ...
    @property
    def NumShells(self) -> int:
        """
        Number of shells, setting will erase old shells
        """
    @NumShells.setter
    def NumShells(self, arg1: int) -> None:
        ...
    @property
    def Oname(self) -> str:
        """
        Output name
        """
    @Oname.setter
    def Oname(self, arg1: str) -> None:
        ...
    @property
    def Plasma(self) -> Plasma:
        ...
    @property
    def PsiGrid(self) -> PsiGrid:
        ...
    @property
    def RSname(self) -> str:
        """
        RS name
        """
    @RSname.setter
    def RSname(self, arg1: str) -> None:
        ...
    @property
    def SGname(self) -> str:
        """
        SG name
        """
    @SGname.setter
    def SGname(self, arg1: str) -> None:
        ...
    @property
    def SMname(self) -> str:
        """
        SM name
        """
    @SMname.setter
    def SMname(self, arg1: str) -> None:
        ...
    @property
    def Seps(self) -> Separatricies:
        """
        Get the separatrixes
        """
    @property
    def Shells(self) -> Shells:
        """
        Get the shells
        """
    @property
    def Start(self) -> str:
        ...
    @property
    def Stop(self) -> str:
        ...
class MatrixView:
    pass
class MeasType:
    """
    Members:
    
      unk
    
      bp
    
      press
    
      pperp
    
      ppar
    
      flux
    
      saddle
    
      circle
    
      coilcur
    
      plasmacur
    
      bt
    
      diam
    
      bangle
    
      flowt
    
      flowp
    
      ne
    
      Te
    
      Ti
    
      Zeff
    
      rho
    
      rot
    
      ppsix
    
      bpangle
    
      pnorm
    
      J0
    """
    J0: typing.ClassVar[MeasType]  # value = <MeasType.J0: 24>
    Te: typing.ClassVar[MeasType]  # value = <MeasType.Te: 16>
    Ti: typing.ClassVar[MeasType]  # value = <MeasType.Ti: 18>
    Zeff: typing.ClassVar[MeasType]  # value = <MeasType.Zeff: 17>
    __members__: typing.ClassVar[dict[str, MeasType]]  # value = {'unk': <MeasType.unk: 0>, 'bp': <MeasType.bp: 1>, 'press': <MeasType.press: 2>, 'pperp': <MeasType.pperp: 3>, 'ppar': <MeasType.ppar: 4>, 'flux': <MeasType.flux: 5>, 'saddle': <MeasType.saddle: 6>, 'circle': <MeasType.circle: 7>, 'coilcur': <MeasType.coilcur: 8>, 'plasmacur': <MeasType.plasmacur: 9>, 'bt': <MeasType.bt: 10>, 'diam': <MeasType.diam: 11>, 'bangle': <MeasType.bangle: 12>, 'flowt': <MeasType.flowt: 13>, 'flowp': <MeasType.flowp: 14>, 'ne': <MeasType.ne: 15>, 'Te': <MeasType.Te: 16>, 'Ti': <MeasType.Ti: 18>, 'Zeff': <MeasType.Zeff: 17>, 'rho': <MeasType.rho: 19>, 'rot': <MeasType.rot: 20>, 'ppsix': <MeasType.ppsix: 21>, 'bpangle': <MeasType.bangle: 12>, 'pnorm': <MeasType.pnorm: 23>, 'J0': <MeasType.J0: 24>}
    bangle: typing.ClassVar[MeasType]  # value = <MeasType.bangle: 12>
    bp: typing.ClassVar[MeasType]  # value = <MeasType.bp: 1>
    bpangle: typing.ClassVar[MeasType]  # value = <MeasType.bangle: 12>
    bt: typing.ClassVar[MeasType]  # value = <MeasType.bt: 10>
    circle: typing.ClassVar[MeasType]  # value = <MeasType.circle: 7>
    coilcur: typing.ClassVar[MeasType]  # value = <MeasType.coilcur: 8>
    diam: typing.ClassVar[MeasType]  # value = <MeasType.diam: 11>
    flowp: typing.ClassVar[MeasType]  # value = <MeasType.flowp: 14>
    flowt: typing.ClassVar[MeasType]  # value = <MeasType.flowt: 13>
    flux: typing.ClassVar[MeasType]  # value = <MeasType.flux: 5>
    ne: typing.ClassVar[MeasType]  # value = <MeasType.ne: 15>
    plasmacur: typing.ClassVar[MeasType]  # value = <MeasType.plasmacur: 9>
    pnorm: typing.ClassVar[MeasType]  # value = <MeasType.pnorm: 23>
    ppar: typing.ClassVar[MeasType]  # value = <MeasType.ppar: 4>
    pperp: typing.ClassVar[MeasType]  # value = <MeasType.pperp: 3>
    ppsix: typing.ClassVar[MeasType]  # value = <MeasType.ppsix: 21>
    press: typing.ClassVar[MeasType]  # value = <MeasType.press: 2>
    rho: typing.ClassVar[MeasType]  # value = <MeasType.rho: 19>
    rot: typing.ClassVar[MeasType]  # value = <MeasType.rot: 20>
    saddle: typing.ClassVar[MeasType]  # value = <MeasType.saddle: 6>
    unk: typing.ClassVar[MeasType]  # value = <MeasType.unk: 0>
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class Measure:
    Fit: float
    Now: float
    R: float
    StdDev: float
    Value: float
    Z: float
    def __init__(self, arg0: int) -> None:
        """
        Create Measurement
        """
    @property
    def Angle(self) -> float:
        """
        Angle
        """
    @Angle.setter
    def Angle(self, arg1: float) -> None:
        ...
    @property
    def CircleType(self) -> CircleType:
        """
        CircleType
        """
    @CircleType.setter
    def CircleType(self, arg1: CircleType) -> None:
        ...
    @property
    def CoilNum(self) -> int:
        """
        CoilNum
        """
    @CoilNum.setter
    def CoilNum(self, arg1: int) -> None:
        ...
    @property
    def Name(self) -> str:
        """
        Name of the measurement
        """
    @Name.setter
    def Name(self, arg1: str) -> None:
        ...
    @property
    def Number(self) -> int:
        """
        Number
        """
    @Number.setter
    def Number(self, arg1: int) -> None:
        ...
    @property
    def R1(self) -> float:
        """
        Saddle Loop R1
        """
    @R1.setter
    def R1(self, arg1: float) -> None:
        ...
    @property
    def R2(self) -> float:
        """
        Saddle Loop R2
        """
    @R2.setter
    def R2(self, arg1: float) -> None:
        ...
    @property
    def Radius(self) -> float:
        """
        Radius
        """
    @Radius.setter
    def Radius(self, arg1: float) -> None:
        ...
    @property
    def Type(self) -> MeasType:
        """
        Type of the measurement
        """
    @Type.setter
    def Type(self, arg1: MeasType) -> None:
        ...
    @property
    def Z1(self) -> float:
        """
        Saddle Loop Z1
        """
    @Z1.setter
    def Z1(self, arg1: float) -> None:
        ...
    @property
    def Z2(self) -> float:
        """
        Saddle Loop Z2
        """
    @Z2.setter
    def Z2(self, arg1: float) -> None:
        ...
class Measures:
    def __getitem__(self, arg0: int) -> Measure:
        ...
    def __len__(self) -> int:
        ...
    def __setitem__(self, arg0: int, arg1: Measure) -> None:
        ...
    def new_meas(self, arg0: MeasType) -> Measure:
        """
        Add a new measurement of type mtype
        """
class ModelType:
    """
    Members:
    
      Std
    
      IsoNoFlow
    
      IsoFlow
    
      AnisoNoFlow
    
      AnisoFlow
    
      DipoleStd
    
      DipoleIntStable
    """
    AnisoFlow: typing.ClassVar[ModelType]  # value = <ModelType.AnisoFlow: 4>
    AnisoNoFlow: typing.ClassVar[ModelType]  # value = <ModelType.AnisoNoFlow: 3>
    DipoleIntStable: typing.ClassVar[ModelType]  # value = <ModelType.DipoleIntStable: 6>
    DipoleStd: typing.ClassVar[ModelType]  # value = <ModelType.DipoleStd: 5>
    IsoFlow: typing.ClassVar[ModelType]  # value = <ModelType.IsoFlow: 2>
    IsoNoFlow: typing.ClassVar[ModelType]  # value = <ModelType.IsoNoFlow: 1>
    Std: typing.ClassVar[ModelType]  # value = <ModelType.Std: 0>
    __members__: typing.ClassVar[dict[str, ModelType]]  # value = {'Std': <ModelType.Std: 0>, 'IsoNoFlow': <ModelType.IsoNoFlow: 1>, 'IsoFlow': <ModelType.IsoFlow: 2>, 'AnisoNoFlow': <ModelType.AnisoNoFlow: 3>, 'AnisoFlow': <ModelType.AnisoFlow: 4>, 'DipoleStd': <ModelType.DipoleStd: 5>, 'DipoleIntStable': <ModelType.DipoleIntStable: 6>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class Plasma:
    G2pTerms: int
    HTerms: int
    Nsize: int
    NumBndMomts: int
    NumPsiPts: int
    PpTerms: int
    RotTerms: int
    SisoTerms: int
    SparTerms: int
    SperTerms: int
    def __init__(self) -> None:
        """
        Create Plasma
        """
    def init(self) -> None:
        """
        Initialize Plasma
        """
    def plasmaG(self, arg0: float) -> float:
        """
        Calculate plasma G
        """
    def plasmaG2p(self, arg0: float) -> float:
        """
        Calculate plasma G2prime
        """
    def plasmaP(self, arg0: float) -> float:
        """
        Calculate plasma pressure
        """
    def plasmaPp(self, arg0: float) -> float:
        """
        Calculate plasma Pprime
        """
    @property
    def Alpha(self) -> MatrixView:
        ...
    @property
    def B0(self) -> float:
        """
        Vacuum magnetic field at R0, Z0
        """
    @B0.setter
    def B0(self, arg0: float) -> None:
        ...
    @property
    def B0R0(self) -> float:
        """
        B0 * R0
        """
    @B0R0.setter
    def B0R0(self, arg0: float) -> None:
        ...
    @property
    def B2(self) -> MatrixView:
        ...
    @property
    def Bt(self) -> MatrixView:
        ...
    @property
    def ChiSqr(self) -> float:
        """
        Chi squared
        """
    @property
    def CrossSection(self) -> float:
        """
        Cross section
        """
    @property
    def Diamag(self) -> float:
        """
        Diamagnetism
        """
    @property
    def Elongation(self) -> float:
        """
        Elongation
        """
    @property
    def G(self) -> MatrixView:
        ...
    @property
    def G2p(self) -> VectorView:
        ...
    @property
    def GradPsi2(self) -> MatrixView:
        ...
    @property
    def GradPsiX(self) -> MatrixView:
        ...
    @property
    def GradPsiZ(self) -> MatrixView:
        ...
    @property
    def H(self) -> VectorView:
        ...
    @property
    def HalfWidth(self) -> float:
        """
        Half width
        """
    @property
    def Ip(self) -> float:
        """
        Plasma current
        """
    @property
    def Ip0(self) -> float:
        """
        initial plasma current
        """
    @Ip0.setter
    def Ip0(self, arg0: float) -> None:
        ...
    @property
    def Jedge(self) -> float:
        """
        Edge current density
        """
    @Jedge.setter
    def Jedge(self, arg0: float) -> None:
        ...
    @property
    def Ltotal(self) -> float:
        """
        Total inductance
        """
    @property
    def Model(self) -> CPlasmaModel:
        ...
    @property
    def ModelType(self) -> ModelType:
        """
        Plasma model type, see ModelType enum
        """
    @ModelType.setter
    def ModelType(self, arg1: ModelType) -> None:
        ...
    @property
    def Perimeter(self) -> float:
        """
        Perimeter
        """
    @property
    def Piso(self) -> MatrixView:
        ...
    @property
    def Pp(self) -> VectorView:
        ...
    @property
    def Ppar(self) -> MatrixView:
        ...
    @property
    def Pper(self) -> MatrixView:
        ...
    @property
    def PsiAxis(self) -> float:
        """
        Psi at axis
        """
    @property
    def PsiLim(self) -> float:
        """
        Psi at plasma/vacuum boundary
        """
    @property
    def PsiMagAxis(self) -> float:
        """
        Psi at magnetic axis
        """
    @property
    def PsiXmax(self) -> float:
        """
        Outermost normalized Psi from 0.0 to 1.0
        """
    @PsiXmax.setter
    def PsiXmax(self, arg0: float) -> None:
        ...
    @property
    def R0(self) -> float:
        """
        Reference major radius
        """
    @R0.setter
    def R0(self, arg0: float) -> None:
        ...
    @property
    def RMagAxis(self) -> float:
        """
        Magnetic axis R
        """
    @property
    def Rho(self) -> MatrixView:
        ...
    @property
    def Rot(self) -> VectorView:
        ...
    @property
    def Siso(self) -> VectorView:
        ...
    @property
    def Spar(self) -> VectorView:
        ...
    @property
    def Sper(self) -> VectorView:
        ...
    @property
    def Volume(self) -> float:
        """
        Volume
        """
    @property
    def Z0(self) -> float:
        """
        Reference vertical position
        """
    @Z0.setter
    def Z0(self, arg0: float) -> None:
        ...
    @property
    def ZMagAxis(self) -> float:
        """
        Magnetic axis Z
        """
    @property
    def beta(self) -> float:
        """
        Average toroidal beta
        """
    @property
    def beta0(self) -> float:
        """
        Vacuum toroidal beta at R0
        """
    @property
    def betap(self) -> float:
        """
        Poloidal beta
        """
    @property
    def li(self) -> float:
        """
        Normalized internal inductance
        """
    @property
    def mu(self) -> float:
        """
        Normalized diamagnetism
        """
    @property
    def q0(self) -> float:
        """
        Central safety factor
        """
    @property
    def qCircular(self) -> float:
        """
        qCircular
        """
    @property
    def qStar(self) -> float:
        """
        qStar
        """
    @property
    def totKinEnergy(self) -> float:
        """
        Total kinetic energy
        """
    @property
    def totMagEnergy(self) -> float:
        """
        Total agnetic energy
        """
class PsiGrid:
    BoundError: float
    BoundThreshold: float
    DelPsi: float
    MaxRes: float
    Nsize: int
    PastMaxRes: float
    RMagAxis: float
    ResThreshold: float
    Rmax: float
    Rmin: float
    Symmetric: int
    UnderRelax1: float
    UnderRelax2: float
    ZMagAxis: float
    Zmax: float
    Zmin: float
    dr: float
    dz: float
    def __init__(self) -> None:
        """
        Create PsiGrid
        """
    def get_IsPlasma(self, arg0: float, arg1: float) -> float:
        """
        Get IsPlasma
        """
    def get_Psi(self, arg0: numpy.ndarray[numpy.float64], arg1: numpy.ndarray[numpy.float64]) -> typing.Any:
        """
        Get Psi
        """
    def get_new_residual(self) -> None:
        """
        Get new residual
        """
    def go_PDE(self) -> None:
        """
        Solve the PDE on the PsiGrid
        """
    def init(self) -> None:
        """
        Initialize PsiGrid
        """
    def init_J(self, arg0: Plasma) -> None:
        """
        Initialize J
        """
    def make_psi_symmetric(self) -> None:
        """
        Make Psi symmetric
        """
    def new_M_solution(self) -> None:
        """
        Get new M solution
        """
    def new_solution(self) -> None:
        """
        Get new solution
        """
    @property
    def Current(self) -> MatrixView:
        ...
    @property
    def IsPlasma(self) -> IMatrixView:
        ...
    @property
    def Psi(self) -> MatrixView:
        ...
    @property
    def PsiAxis(self) -> float:
        """
        Psi at FCFS or Magnetic Axis
        """
    @PsiAxis.setter
    def PsiAxis(self, arg0: float) -> None:
        ...
    @property
    def PsiLim(self) -> float:
        """
        Psi at plasma/vacuum boundary
        """
    @PsiLim.setter
    def PsiLim(self, arg0: float) -> None:
        ...
    @property
    def PsiMagAxis(self) -> float:
        """
        Psi at Magnetic Axis
        """
    @PsiMagAxis.setter
    def PsiMagAxis(self, arg0: float) -> None:
        ...
    @property
    def R(self) -> VectorView:
        ...
    @property
    def Residual(self) -> MatrixView:
        ...
    @property
    def Z(self) -> VectorView:
        ...
class Separatricies:
    def __getitem__(self, arg0: int) -> Separatrix:
        ...
    def __len__(self) -> int:
        ...
    def __setitem__(self, arg0: int, arg1: Separatrix) -> None:
        ...
    def new_separatrix(self) -> Separatrix:
        """
        Add a new separatrix
        """
class Separatrix:
    Enabled: int
    IsSeparatrix: int
    def __init__(self) -> None:
        """
        Create Separatrix
        """
    @property
    def Name(self) -> str:
        """
        Name of the separatrix
        """
    @Name.setter
    def Name(self, arg1: str) -> None:
        ...
    @property
    def Psi(self) -> float:
        """
        Value of Psi at separatrix
        """
    @property
    def R1(self) -> float:
        """
        Box to search for separatrix
        """
    @R1.setter
    def R1(self, arg0: float) -> None:
        ...
    @property
    def R2(self) -> float:
        """
        Box to search for separatrix
        """
    @R2.setter
    def R2(self, arg0: float) -> None:
        ...
    @property
    def RC(self) -> float:
        """
        Center of plasma from Sep
        """
    @RC.setter
    def RC(self, arg0: float) -> None:
        ...
    @property
    def Rs(self) -> float:
        """
        R of separatrix
        """
    @property
    def Z1(self) -> float:
        """
        Box to search for separatrix
        """
    @Z1.setter
    def Z1(self, arg0: float) -> None:
        ...
    @property
    def Z2(self) -> float:
        """
        Box to search for separatrix
        """
    @Z2.setter
    def Z2(self, arg0: float) -> None:
        ...
    @property
    def ZC(self) -> float:
        """
        Center of plasma from Sep
        """
    @ZC.setter
    def ZC(self, arg0: float) -> None:
        ...
    @property
    def Zs(self) -> float:
        """
        Z of separatrix
        """
class Shell:
    Enabled: int
    def __init__(self, arg0: int) -> None:
        """
        Create Shell
        """
    def set_NumSubShells(self, arg0: int, arg1: Machine) -> None:
        """
        Set the number of subshells
        """
    @property
    def Name(self) -> str:
        """
        Name of the shell
        """
    @Name.setter
    def Name(self, arg1: str) -> None:
        ...
    @property
    def NumSubShells(self) -> int:
        ...
    @property
    def SubShells(self) -> SubShells:
        """
        Return vector of SubShells
        """
class Shells:
    def __getitem__(self, arg0: int) -> Shell:
        ...
    def __len__(self) -> int:
        ...
    def __setitem__(self, arg0: int, arg1: Shell) -> None:
        ...
    def new_shell(self) -> Shell:
        """
        Add a new shell
        """
class SubCoil:
    R: float
    Z: float
    def __init__(self) -> None:
        """
        Create SubCoil
        """
    @property
    def Fraction(self) -> float:
        """
        Fraction of current
        """
    @Fraction.setter
    def Fraction(self, arg0: float) -> None:
        ...
    @property
    def Name(self) -> str:
        """
        Name of the subcoil
        """
    @Name.setter
    def Name(self, arg1: str) -> None:
        ...
class SubCoils:
    def __getitem__(self, arg0: int) -> SubCoil:
        ...
    def __len__(self) -> int:
        ...
    def __setitem__(self, arg0: int, arg1: SubCoil) -> None:
        ...
    def new_subcoil(self) -> SubCoil:
        """
        Add a new subcoil
        """
class SubShell:
    Current: float
    R: float
    Radius: float
    Z: float
    def __init__(self) -> None:
        """
        Create SubShell
        """
    @property
    def Name(self) -> str:
        """
        Name of the subshell
        """
    @Name.setter
    def Name(self, arg1: str) -> None:
        ...
class SubShells:
    def __getitem__(self, arg0: int) -> SubShell:
        ...
    def __len__(self) -> int:
        ...
class VectorView:
    pass
