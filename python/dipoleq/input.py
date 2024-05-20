"""
dipoleq input file schema using pydantic for validation
"""

from typing import Optional, List, Any, Union
from typing_extensions import Self

from pydantic import BaseModel, Field, field_validator, model_validator

from ._c import Machine, PsiGrid, Plasma, Coil, Shell, SubCoil, SubShell, \
                Limiter, Separatrix, Measure, ModelType, MeasType, CircleType

class PsiGridIn(BaseModel):
    Nsize: Optional[int] = 256
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
        pg.Symmetric = self.Symmetric
        pg.BoundThreshold = self.BoundThreshold
        pg.ResThreshold = self.ResThreshold
        pg.UnderRelax1 = self.UnderRelax1
        pg.UnderRelax2 = self.UnderRelax2


class PlasmaIn(BaseModel):
    B0: float
    R0: float
    Ip0: float
    NumPsiPts: int
    PsiXmax: float
    Jedge: float
    # ModelType: Union[ModelType, int]  # ModelType
    ModelType: int
    RPeak: Optional[float]
    ZPeak: Optional[float]
    PEdge: Optional[int]
    PsiFlat: Optional[float]
    NSurf: Optional[int]
    fCrit: Optional[int]

    @field_validator('ModelType')
    @classmethod
    def check_model_type(cls, v: Any) -> Any:
        if isinstance(v, int):
            return ModelType(v)
        return v

    def do_init(self, pl: Plasma) -> None:
        pl.B0 = self.B0
        pl.R0 = self.R0
        pl.Ip0 = self.Ip0
        pl.NumPsiPts = self.NumPsiPts
        pl.PsiXmax = self.PsiXmax
        pl.Jedge = self.Jedge
        pl.ModelType = ModelType(self.ModelType)
        pl.RPeak = self.RPeak
        pl.ZPeak = self.ZPeak
        pl.PEdge = self.PEdge
        pl.PsiFlat = self.PsiFlat
        pl.NSurf = self.NSurf
        pl.fCrit = self.fCrit


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
        if not self.R or self.R < 0.0:
            if not self.SubCoils:
                raise ValueError("Coil must have subcoils, or defined limit points")
        else:
            if self.dR == 0.0 or self.dZ == 0.0:
                raise ValueError("For coils without subcoils, dR and dZ must be defined")
        return self

    def do_init(self, coil: Coil) -> None:
        if self.Name:
            coil.Name = self.Name
        coil.Enabled = self.Enabled
        coil.InitialCurrent = self.InitialCurrent
        coil.R = self.R
        coil.dR = self.dR
        coil.Z = self.Z
        coil.dZ = self.dZ
        # these should be allocated earlier by machine init
        if coil.NumSubCoils != self.NumSubCoils:
            coil.NumSubCoils = self.NumSubCoils
        if self.SubCoils:
            for selfsc, subcoil in zip(self.SubCoils, coil.SubCoils):
                selfsc.do_init(scoil)


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
        shell.Enabled = self.Enabled
        if shell.NumSubShells != self.NumSubShells:
            shell.NumSubShells = self.NumSubShells
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

    def init_machine(self) -> Machine:
        
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
        
        # now do lists
        for i, coil in enumerate(self.Coils):
            coil.do_init(m.Coils[i])
            
        for i, lim in enumerate(self.Limiters):
            lim.do_init(m.Limiters[i])
            
        if self.Separatricies:
            for i, sep in enumerate(self.Separatricies):
                sep.do_init(m.Separatricies[i])
                
        if self.Measures:
            for i, meas in enumerate(self.Measures):
                meas.do_init(m.Measures[i])
                
        if self.Shells:
            for i, shell in enumerate(self.Shells):
                shell.do_init(m.Shells[i])
        
        return m
    
