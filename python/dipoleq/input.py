"""
dipoleq input file schema using pydantic for validation
"""

from re import M
from typing import Optional, List, Any, Union, Tuple
from unittest.mock import Base
from numpy import isin
from typing_extensions import Self, Annotated, Literal
import enum
from itertools import chain

from pydantic import BaseModel, Field, field_validator, model_validator

from ._c import Machine, PsiGrid, Plasma, CPlasmaModel, \
                Coil, Shell, SubCoil, SubShell, \
                Limiter, Separatrix, Measure, ModelType, MeasType, CircleType


MU0 = 4.0e-7*3.14159265358979323846


class PsiGridIn(BaseModel):
    Nsize: int = Field(default=256)
    Rmin: float = Field(alias="Xmin")
    Rmax: float = Field(alias="Xmax")
    Zmin: float
    Zmax: float
    Symmetric: Optional[bool] = False
    BoundThreshold: float
    ResThreshold: float
    UnderRelax1: float
    UnderRelax2: float
    
    def do_init(self, pg: PsiGrid) -> None:
        pg.Nsize = self.Nsize
        pg.Rmin = self.Rmin
        pg.Rmax = self.Rmax
        pg.Zmin = self.Zmin
        pg.Zmax = self.Zmax
        pg.Symmetric = 1 if self.Symmetric else 0
        pg.BoundThreshold = self.BoundThreshold
        pg.ResThreshold = self.ResThreshold
        pg.UnderRelax1 = self.UnderRelax1
        pg.UnderRelax2 = self.UnderRelax2


def model_type_ids(mts: List[ModelType]) -> List[Union[int, str]]:
    """ make a list of literals for the model types both name and value
    """
    lists_of_literals = [[mt.value, mt.name, str(mt.value)] for mt in mts]
    return list(chain.from_iterable(lists_of_literals))


OLD_MODELS = [ModelType.Std,
              ModelType.IsoNoFlow, ModelType.IsoFlow,
              ModelType.AnisoNoFlow, ModelType.AnisoFlow]

PLASMA_BASE_MODEL_TYPES = Literal[*model_type_ids(OLD_MODELS)]


class PlasmaModelBaseModel(BaseModel):
    Type: Any

    @field_validator('Type', mode='after')
    @classmethod
    def check_model_type(cls, v: Any) -> Any:
        if isinstance(v, str):
            if v.isnumeric():
                return ModelType(int(v))
            return ModelType.__members__[v]
        if isinstance(v, int):
            return ModelType(v)
        return v
        

class PlasmaModelOld(PlasmaModelBaseModel):
    Type: PLASMA_BASE_MODEL_TYPES = Field(alias="ModelType")
    G2pTerms: Optional[int]
    HTerms: Optional[int]
    PpTerms: Optional[int]
    HTerms: Optional[int]
    RotTerms: Optional[int]
    SisoTerms: Optional[int]
    SparTerms: Optional[int]
    SperTerms: Optional[int]    
    G2p: Optional[List[float]]
    H: Optional[List[float]]
    Pp: Optional[List[float]]
    Rot: Optional[List[float]]
    Siso: Optional[List[float]]
    Spar: Optional[List[float]]
    Sper: Optional[List[float]]
   
    @model_validator(mode='after')
    def check_old_plasma_model(self) -> Self:
        # should also check for consistency of the terms and models
        if self.Type not in OLD_MODELS:
            raise ValueError(f"ModelType {self.Type} is not a valid plasma model")
        return self
    
    def do_init(self, pm: Plasma) -> None:
        pm.G2pTerms = len(self.G2p) if self.G2p else 0
        if self.G2p:
            for i, v in enumerate(self.G2p):
                pm.G2p[i] = v
        pm.HTerms = len(self.H) if self.H else 0
        if self.H:
            for i, v in enumerate(self.H):
                pm.H[i] = v
        pm.PpTerms = len(self.Pp) if self.Pp else 0
        if self.Pp:
            for i, v in enumerate(self.Pp):
                pm.Pp[i] = v
        pm.RotTerms = len(self.Rot) if self.Rot else 0
        if self.Rot:
            for i, v in enumerate(self.Rot):
                pm.Rot[i] = v
        pm.SisoTerms = len(self.Siso) if self.Siso else 0
        if self.Siso:
            for i, v in enumerate(self.Siso):
                pm.Siso[i] = v
        pm.SparTerms = len(self.Spar) if self.Spar else 0
        if self.Spar:
            for i, v in enumerate(self.Spar):
                pm.Spar[i] = v
        pm.SperTerms = len(self.Sper) if self.Sper else 0
        if self.Sper:
            for i, v in enumerate(self.Sper):
                pm.Sper[i] = v


class CDipoleStdIn(PlasmaModelBaseModel):
    Type: Literal[*model_type_ids([ModelType.DipoleStd])] = \
        Field(ModelType.DipoleStd, alias="ModelType")
    RPeak: Optional[float]
    ZPeak: Optional[float]
    PsiPeak: Optional[float]
    PsiEdge: Optional[float]
    PsiDipole: Optional[float]
    PPeak: Optional[float]
    PrExp: Optional[float]
    def do_init(self, pl: Plasma) -> None:
        pm = pl.Model
        keys = ['RPeak', 'ZPeak', 'PsiPeak', 'PsiEdge',
                'PsiDipole', 'PPeak', 'PrExp']
        for key in keys:
            pm.model_input(key, str(getattr(self, key)), "")


DIPOLE_INTSTABLE_IDS = Annotated[Literal[*model_type_ids([ModelType.DipoleIntStable])],
                                 Field(ModelType.DipoleIntStable,
                                       alias="ModelType")]


class CDipoleIntStableIn(PlasmaModelBaseModel):
    Type: DIPOLE_INTSTABLE_IDS
    RPeak: Optional[float]
    ZPeak: Optional[float]
    PEdge: Optional[float]
    PsiFlat: Optional[float]
    NSurf: Optional[int]
    fCrit: Optional[float]

    def do_init(self, pl: Plasma) -> None:
        keys = ['RPeak', 'ZPeak', 'PEdge', 'PsiFlat', 'NSurf', 'fCrit']
        pm = pl.Model
        for key in keys:
            pm.model_input(key, str(getattr(self, key)), "")    


PlasmaModel = Annotated[Union[PlasmaModelOld, CDipoleStdIn,
                              CDipoleIntStableIn],
                        Field(discriminator='Type')]

class PlasmaIn(BaseModel):
    B0: float
    R0: float
    R0B0: Optional[float] = None
    Ip0: Optional[float]
    NumBndMomts: Optional[int] = None
    NumPsiPts: Optional[int]
    PsiXmax: Optional[float]
    Jedge: Optional[float]
    Model: PlasmaModel
    
    #@field_validator('ModelType')
    #@classmethod
    #def check_model_type(cls, v: Any) -> Any:
    #    if isinstance(v, int):
    #        return ModelType(v)
    #    return v

    def do_init(self, pl: Plasma) -> None:
        pl.B0 = self.B0
        pl.R0 = self.R0
        if self.R0B0:
            pl.R0B0 = self.R0B0
        if self.Ip0:
            pl.Ip0 = self.Ip0
        if self.NumBndMomts:
            pl.NumBndMomts = self.NumBndMomts
        if self.NumPsiPts:
            pl.NumPsiPts = self.NumPsiPts
        if self.PsiXmax:
            pl.PsiXmax = self.PsiXmax
        if self.Jedge:
            pl.Jedge = self.Jedge
        pl.ModelType = ModelType(self.Model.Type)

class LimiterIn(BaseModel):
    Name: Optional[str]
    R1: float = Field(alias="X1")
    Z1: float
    R2: float = Field(alias="X2")
    Z2: float
    Enabled: Optional[int] = True
    
    def do_init(self, lim: Limiter) -> None:
        if self.Name:
            lim.Name = self.Name
        lim.R1 = self.R1
        lim.Z1 = self.Z1
        lim.R2 = self.R2
        lim.Z2 = self.Z2
        lim.Enabled = self.Enabled


class SeparatrixIn(BaseModel):
    Name: Optional[str]
    R1: float = Field(alias="X1")
    Z1: float
    R2: float = Field(alias="X2")
    Z2: float
    RC: float = Field(alias="XC")
    ZC: float
    Enabled: Optional[bool] = True
    
    def do_init(self, sep: Separatrix) -> None:
        if self.Name:
            sep.Name = self.Name
        sep.R1 = self.R1
        sep.Z1 = self.Z1
        sep.R2 = self.R2
        sep.Z2 = self.Z2
        sep.RC = self.RC
        sep.ZC = self.ZC
        sep.Enabled = self.Enabled


class SubCoilIn(BaseModel):
    Name: Optional[str]
    Enabled: Optional[bool] = True
    Fraction: float
    R: float = Field(alias="X")
    Z: float

    def do_init(self, scoil: SubCoil) -> None:
        if self.Name:
            scoil.Name = self.Name
        scoil.Enabled = self.Enabled
        scoil.Fraction = self.Fraction
        scoil.R = self.R
        scoil.Z = self.Z


class CoilIn(BaseModel):
    Name: Optional[str]
    Enabled: Optional[bool] = True
    InitialCurrent: float
    R: Optional[float] = Field(alias="X")
    dR: Optional[float] = Field(alias="dX")
    Z: Optional[float]
    dZ: Optional[float]
    NumSubCoils: Optional[int] = None
    SubCoils: Optional[List[SubCoilIn]] = None

    @model_validator(mode='after')
    def check_shape_vs_subcoils(self) -> Self:
        if self.NumSubCoils:
            if self.NumSubCoils != (len(self.SubCoils) if self.SubCoils else 0):
                raise ValueError(f"NumSubCoils {self.NumSubCoils} does not match "
                                 f"the number of subcoils {len(self.SubCoils)}")
        else:
            self.NumSubCoils = len(self.SubCoils) if self.SubCoils else 0
        if not self.R or self.R < 0.0:
            if not self.SubCoils:
                raise ValueError("Coil must have subcoils, or defined limit points")
        else:
            if self.dR == 0.0 or self.dZ == 0.0:
                raise ValueError("For coils without subcoils, dR and dZ must be defined")
        return self

    def do_init(self, coil: Coil, psiGrid: PsiGrid) -> None:
        if self.Name:
            coil.Name = self.Name
        coil.Enabled = self.Enabled
        coil.CoilCurrent = self.InitialCurrent*MU0
        if self.R:
            coil.R = self.R
        if self.dR:
            coil.dR = self.dR
        if self.Z:
            coil.Z = self.Z
        if self.dZ:
            coil.dZ = self.dZ
        # these should be allocated earlier by machine init
        if self.SubCoils:
            self.NumSubCoils = len(self.SubCoils)
        if coil.NumSubCoils != self.NumSubCoils:
            if self.NumSubCoils:
                coil.NumSubCoils = self.NumSubCoils
        if self.SubCoils:
            for selfsc, subcoil in zip(self.SubCoils, coil.SubCoils):
                selfsc.do_init(subcoil)
        if self.dR > 0.0:
            coil.compute_SubCoils(psiGrid)


class SubShellIn(BaseModel):
    Name: Optional[str]
    Current: Optional[float] = 0.0
    R: float = Field(alias="X")
    Z: float

    def do_init(self, sshell: SubShell) -> None:
        if self.Name:
            sshell.Name = self.Name
        sshell.Current = self.Current
        sshell.R = self.R
        sshell.Z = self.Z


class ShellIn(BaseModel):
    Name: Optional[str]
    Enabled: Optional[bool] = True
    NumSubShells: Optional[int]
    SubShells: List[SubShellIn]

    def do_init(self, shell: Shell) -> None:
        if self.Name:
            shell.Name = self.Name
        if self.Enabled:
            shell.Enabled = self.Enabled
        if shell.NumSubShells != self.NumSubShells:
            raise ValueError(f"NumSubShells {self.NumSubShells} does not match"
                             f" the number of subshells {shell.NumSubShells}")
        for selfss, subshell in zip(self.SubShells, shell.SubShells):
            selfss.do_init(subshell)


class MeasureIn(BaseModel):
    Name: Optional[str]
    Enabled: Optional[bool] = True
    Type: str

    def do_init(self, meas: Measure) -> None:
        if self.Name:
            meas.Name = self.Name
        meas.Enabled = self.Enabled
        meas.Type = self.Type

class MachineIn(BaseModel):
    # Complete input file has these fields
    MaxIterFixed: Optional[int] = 0
    MaxIterFree: Optional[int] = 50
    Name: Optional[str] = "pyDipolEQ"
    Info: Optional[str] = "DipolEQ Equilibrium"
    Oname: Optional[str]
    MGname: Optional[str] = ""
    LHname: Optional[str] = ""
    RSname: Optional[str] = ""
    RestartStatus: Optional[bool] = False
    RestartUnkns: Optional[bool] = False
    LHGreenStatus: Optional[bool] = False
    MGreenStatus: Optional[bool] = False
    NumEqualEq: Optional[int] = 0
    VacuumOnly: Optional[bool] = False
    FitMeasurements: Optional[bool] = False
    PsiGrid: PsiGridIn
    Plasma: PlasmaIn
    Coils: List[CoilIn]
    Limiters: List[LimiterIn]
    Separatricies: Optional[List[SeparatrixIn]] = None
    Measures: Optional[List[MeasureIn]] = None
    Shells: Optional[List[ShellIn]] = None
    NumCoils: Optional[int] = None
    NumShells: Optional[int] = None
    NumLimiters: Optional[int] = None
    NumSeps: Optional[int] = None
    NumMeasures: Optional[int] = None

    @model_validator(mode='after')
    def check_numbers(self) -> Self:
        if self.NumCoils:
            if self.NumCoils != len(self.Coils):
                raise ValueError(f"NumCoils {self.NumCoils} does not match the"
                                 f" number of coils {len(self.Coils)}")
        else:
            self.NumCoils = len(self.Coils)
        if self.NumShells:
            numShells = len(self.Shells) if self.Shells else 0
            if self.NumShells != numShells:
                raise ValueError(f"NumShells {self.NumShells} does not match"
                                 f" the number of shells {numShells}")
        else:
            self.NumShells = len(self.Shells) if self.Shells else 0
        if self.NumLimiters:
            if self.NumLimiters != len(self.Limiters):
                raise ValueError(f"NumLimiters {self.NumLimiters} does not "
                                 f"match the number of limiters "
                                 f"{len(self.Limiters)}")
        else:
            self.NumLimiters = len(self.Limiters) if self.Limiters else 0
        if self.NumSeps:
            numSeps = len(self.Separatricies) if self.Separatricies else 0
            if self.NumSeps != numSeps:
                raise ValueError(f"NumSeps {self.NumSeps} does not match the"
                                 f" number of separatricies {numSeps}")
        else:
            self.NumSeps = len(self.Separatricies) if self.Separatricies else 0
        if self.NumMeasures:
            numMeasures = len(self.Measures) if self.Measures else 0
            if self.NumMeasures != numMeasures:
                raise ValueError(f"NumMeasures {self.NumMeasures} does not "
                                 f"match the number of measures {numMeasures}")
        else:
            self.NumMeasures = len(self.Measures) if self.Measures else 0

        if not self.Oname:
            self.Oname = self.Name
            
        return self

    def create_machine(self) -> Machine:
        """Creates a Machine object from the input data

        Returns:
            Machine: The machine object, ready to be solved hopefully
        """
        
        m = Machine()
        
        # names
        m.Name = self.Name
        m.Info = self.Info
        m.Oname = self.Oname
        m.MGname = self.MGname
        m.LHname = self.LHname
        m.RSname = self.RSname
        
        # control parameters
        m.RestartStatus = self.RestartStatus
        m.RestartUnkns = self.RestartUnkns
        m.LHGreenStatus = self.LHGreenStatus
        m.MGreenStatus = self.MGreenStatus
        m.NumEqualEq = self.NumEqualEq
        m.VacuumOnly = self.VacuumOnly
        
        m.MaxIterFixed = self.MaxIterFixed
        m.MaxIterFree = self.MaxIterFree
        
        # set numbers
        m.NumCoils = self.NumCoils
        m.NumShells = self.NumShells
        m.NumLimiters = self.NumLimiters
        m.NumSeps = self.NumSeps
        m.NumMeasures = self.NumMeasures
        
        # init PsiGrid and Plasma
        self.PsiGrid.do_init(m.PsiGrid)
        self.Plasma.do_init(m.Plasma)

        # init machine (do allocations)
        m.init()
        
        # now init the plasma model that was inited
        self.Plasma.Model.do_init(m.Plasma)
        
        # now do lists
        for i, coil in enumerate(self.Coils):
            coil.do_init(m.Coils[i], m.PsiGrid)
            
        for i, lim in enumerate(self.Limiters):
            if not m.Limiters[i]:
                m.Limiters[i] = m.Limiters.new_limiter()
            lim.do_init(m.Limiters[i])
            
        if self.Separatricies:
            for i, sep in enumerate(self.Separatricies):
                if not m.Seps[i]:
                    m.Seps[i] = m.Seps.new_separatrix()
                sep.do_init(m.Seps[i])
                
        if self.Measures:
            for i, meas in enumerate(self.Measures):
                if not m.Measures[i]:
                    m.Measures[i] = m.Measures.new_meas()
                meas.do_init(m.Measures[i])
                
        if self.Shells:
            for i, shell in enumerate(self.Shells):
                shell.do_init(m.Shells[i])
        
        return m
    
