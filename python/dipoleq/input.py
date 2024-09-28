"""
dipoleq input file schema using pydantic for validation
"""

from __future__ import annotations

from functools import reduce
from itertools import chain
from typing import Annotated, Any, Generic, Literal, TypeAlias, TypeVar, cast

from pydantic import (
    AliasChoices,
    BaseModel,
    Field,
    GetCoreSchemaHandler,
    field_validator,
    model_validator,
)
from pydantic_core import core_schema
from typing_extensions import Self

from .core import (
    CircleType,
    Coil,
    Limiter,
    Machine,
    MeasType,
    Measure,
    ModelType,
    Plasma,
    PsiGrid,
    Separatrix,
    Shell,
    SubCoil,
    SubShell,
)

MU0 = 4.0e-7 * 3.14159265358979323846


# deal with pybind11 enum classes
_T = TypeVar("_T", ModelType, MeasType, CircleType)


class PyBindEnumAnnotation(Generic[_T]):
    """Validation annotation for pybind11 enum classes."""

    @classmethod
    def __get_pydantic_core_schema__(
        cls, source: type[Any], handler: GetCoreSchemaHandler
    ) -> core_schema.CoreSchema:
        if source in [
            ModelType,
            MeasType,
            CircleType,
        ]:  # this line can't seem to not use specific choices here

            def validate(val: Any) -> _T:
                if isinstance(val, source):
                    return cast(_T, val)
                if isinstance(val, str):
                    if val.isnumeric():
                        return cast(_T, source(int(val)))
                    return cast(_T, source.__members__[val])
                if isinstance(val, int):
                    return cast(_T, source(val))
                raise ValueError(f"{val} is not convertible to {source}")

            return core_schema.no_info_after_validator_function(
                validate,
                core_schema.enum_schema(source, list(source.__members__.values())),
            )
        return handler.generate_schema(source)


MeasTypeA: TypeAlias = Annotated[MeasType, PyBindEnumAnnotation]
ModelTypeA: TypeAlias = Annotated[ModelType, PyBindEnumAnnotation]
CircleTypeA: TypeAlias = Annotated[CircleType, PyBindEnumAnnotation]


class PsiGridIn(BaseModel):
    Nsize: int = Field(default=256)
    Rmin: float = Field(validation_alias=AliasChoices("Rmin", "Xmin"))
    Rmax: float = Field(validation_alias=AliasChoices("Rmax", "Xmax"))
    Zmin: float
    Zmax: float
    Symmetric: bool | None = False
    BoundThreshold: float = 1.0e-6
    ResThreshold: float = 1.0e-6
    UnderRelax1: float | None = None
    UnderRelax2: float | None = None

    def do_init(self, pg: PsiGrid) -> None:
        pg.Nsize = self.Nsize
        pg.Rmin = self.Rmin
        pg.Rmax = self.Rmax
        pg.Zmin = self.Zmin
        pg.Zmax = self.Zmax
        pg.Symmetric = 1 if self.Symmetric else 0
        if self.BoundThreshold:
            pg.BoundThreshold = self.BoundThreshold
        if self.ResThreshold:
            pg.ResThreshold = self.ResThreshold
        if self.UnderRelax1:
            pg.UnderRelax1 = self.UnderRelax1
        if self.UnderRelax2:
            pg.UnderRelax2 = self.UnderRelax2


def model_type_literals(mts: list[ModelTypeA]) -> TypeAlias:
    """make a list of literals for the model types both name and value"""
    lists_of_vals: list[list[str | int]] = [
        [mt.value, mt.name, str(mt.value)] for mt in mts
    ]
    list_of_vals: list[str | int] = list(chain.from_iterable(lists_of_vals))
    list_of_literals: list[TypeAlias] = [Literal[v] for v in list_of_vals]
    # list_of_literals: list[TypeVar] = [Literal[v] for v in list_of_vals]

    def unionizer(a: TypeAlias, b: TypeAlias) -> TypeAlias:
        return a | b

    return reduce(unionizer, list_of_literals)
    # return Union[list_of_literals]


class PlasmaModelBaseModel(BaseModel):
    """Base class for plasma models. This is not a real model, but a base class"""

    Type: ModelTypeA = Field(
        validation_alias=AliasChoices("Type", "ModelType"), serialization_alias="int"
    )

    @field_validator("Type", mode="after")
    @classmethod
    def check_model_type(cls, v: Any) -> Any:
        if isinstance(v, str):
            if v.isnumeric():
                return ModelTypeA(int(v))
            return ModelTypeA.__members__[v]
        if isinstance(v, int):
            return ModelTypeA(v)
        return v


OLD_MODELS = [
    ModelType.Std,
    ModelType.IsoNoFlow,
    ModelType.IsoFlow,
    ModelType.AnisoNoFlow,
    ModelType.AnisoFlow,
]

PLASMA_BASE_MODEL_TYPES: TypeAlias = model_type_literals(OLD_MODELS)  # type: ignore[valid-type]


class PlasmaModelOld(PlasmaModelBaseModel):
    """Old plasma model data. Covers the standard polynomial models."""

    Type: PLASMA_BASE_MODEL_TYPES
    G2pTerms: int | None
    HTerms: int | None
    PpTerms: int | None
    RotTerms: int | None
    SisoTerms: int | None
    SparTerms: int | None
    SperTerms: int | None
    G2p: list[float] | None
    H: list[float] | None
    Pp: list[float] | None
    Rot: list[float] | None
    Siso: list[float] | None
    Spar: list[float] | None
    Sper: list[float] | None

    @model_validator(mode="after")
    def check_old_plasma_model(self) -> Self:
        """Check for validation of the old plasma model data."""
        # should also check for consistency of the terms and models
        if self.Type not in OLD_MODELS:
            raise ValueError(f"ModelType {self.Type} is not a valid plasma model")
        return self

    def do_init(self, pm: Plasma) -> None:
        """Initialize the plasma model with the input data."""
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


# python < 3.11 can't do argument unpacking in a literal
# so we have to do it manually
LTS: TypeAlias = model_type_literals([ModelType.DipoleStd])  # type: ignore[valid-type]


class CDipoleStdIn(PlasmaModelBaseModel):
    Type: LTS
    RPeak: float | None
    ZPeak: float | None
    PsiPeak: float | None
    PsiEdge: float | None
    PsiDipole: float | None
    PPeak: float | None
    PrExp: float | None

    def do_init(self, pl: Plasma) -> None:
        pm = pl.Model
        keys = ["RPeak", "ZPeak", "PsiPeak", "PsiEdge", "PsiDipole", "PPeak", "PrExp"]
        for key in keys:
            pm.model_input(key, str(getattr(self, key)), "")


LDS: TypeAlias = model_type_literals([ModelType.DipoleIntStable])  # type: ignore[valid-type]


class CDipoleIntStableIn(PlasmaModelBaseModel):
    Type: LDS
    RPeak: float | None
    ZPeak: float | None
    PEdge: float | None
    PsiFlat: float | None
    NSurf: int | None
    fCrit: float | None

    def do_init(self, pl: Plasma) -> None:
        """Initialize the plasma model with the input data."""
        keys = ["RPeak", "ZPeak", "PEdge", "PsiFlat", "NSurf", "fCrit"]
        pm = pl.Model
        for key in keys:
            pm.model_input(key, str(getattr(self, key)), "")


LDPN: TypeAlias = model_type_literals([ModelType.DipoleStablePsiN])  # type: ignore[valid-type]


class CDipoleStablePsiNIn(PlasmaModelBaseModel):
    Type: LDPN
    PsiNPeak: float | None
    PEdge: float | None
    PsiFlat: float | None
    NSurf: int | None
    fCrit: float | None

    def do_init(self, pl: Plasma) -> None:
        """Initialize the plasma model with the input data."""
        keys = ["PsiNPeak", "PEdge", "PsiFlat", "NSurf", "fCrit"]
        pm = pl.Model
        for key in keys:
            pm.model_input(key, str(getattr(self, key)), "")


PlasmaModel = Annotated[
    PlasmaModelOld | CDipoleStdIn | CDipoleIntStableIn | CDipoleStablePsiNIn,
    Field(discriminator="Type"),
]


class PlasmaIn(BaseModel):
    B0: float
    R0: float
    R0B0: float | None = None
    Ip0: float | None
    NumBndMomts: int | None = None
    NumPsiPts: int | None
    PsiXmax: float | None
    Jedge: float | None
    Model: PlasmaModel

    def do_init(self, pl: Plasma) -> None:
        pl.B0 = self.B0
        pl.R0 = self.R0
        if self.R0B0:
            pl.B0R0 = self.R0B0
        else:
            pl.B0R0 = self.R0 * self.B0
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
    Name: str | None
    R1: float = Field(validation_alias=AliasChoices("R1", "X1"))
    Z1: float
    R2: float = Field(validation_alias=AliasChoices("R2", "X2"))
    Z2: float
    Enabled: int | None = True

    def do_init(self, lim: Limiter) -> None:
        if self.Name:
            lim.Name = self.Name
        lim.R1 = self.R1
        lim.Z1 = self.Z1
        lim.R2 = self.R2
        lim.Z2 = self.Z2
        if self.Enabled is not None:
            lim.Enabled = self.Enabled


class SeparatrixIn(BaseModel):
    Name: str | None
    R1: float = Field(default=0.0, validation_alias=AliasChoices("R1", "X1"))
    Z1: float
    R2: float = Field(default=0.0, validation_alias=AliasChoices("R2", "X2"))
    Z2: float
    RC: float = Field(default=0.0, validation_alias=AliasChoices("RC", "XC"))
    ZC: float = 0.0
    Enabled: bool | None = True

    def do_init(self, sep: Separatrix) -> None:
        if self.Name:
            sep.Name = self.Name
        sep.R1 = self.R1
        sep.Z1 = self.Z1
        sep.R2 = self.R2
        sep.Z2 = self.Z2
        sep.RC = self.RC
        sep.ZC = self.ZC
        if self.Enabled is not None:
            sep.Enabled = self.Enabled


class SubCoilIn(BaseModel):
    """Subcoil input data."""

    Name: str | None = None
    Fraction: float
    R: float = Field(validation_alias=AliasChoices("R", "X"))
    Z: float

    def do_init(self, scoil: SubCoil) -> None:
        if self.Name:
            scoil.Name = self.Name
        scoil.Fraction = self.Fraction
        scoil.R = self.R
        scoil.Z = self.Z


class CoilIn(BaseModel):
    """Coil input data."""

    Name: str | None = None
    Enabled: bool | None = True
    InitialCurrent: float
    R: float | None = Field(default=None, validation_alias=AliasChoices("R", "X"))
    dR: float | None = Field(default=None, validation_alias=AliasChoices("dR", "dX"))
    Z: float | None = None
    dZ: float | None = None
    NumSubCoils: int | None = None
    SubCoils: list[SubCoilIn] | None = None

    @model_validator(mode="after")
    def check_shape_vs_subcoils(self) -> Self:
        if self.NumSubCoils:
            if self.NumSubCoils != (len(self.SubCoils) if self.SubCoils else 0):
                raise ValueError(
                    f"NumSubCoils {self.NumSubCoils} does not match "
                    f"the number of subcoils {len(self.SubCoils) if self.SubCoils else 0}"
                )
        else:
            self.NumSubCoils = len(self.SubCoils) if self.SubCoils else 0
        if not self.R or self.R < 0.0:
            if not self.SubCoils:
                raise ValueError("Coil must have subcoils, or defined limit points")
        elif self.dR == 0.0 or self.dZ == 0.0:
            raise ValueError("For coils without subcoils, dR and dZ must be defined")
        return self

    def do_init(self, coil: Coil, psiGrid: PsiGrid) -> None:
        if self.Name:
            coil.Name = self.Name
        if self.Enabled is not None:
            coil.Enabled = int(self.Enabled)
        coil.CoilCurrent = self.InitialCurrent * MU0
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
        if coil.NumSubCoils != self.NumSubCoils and self.NumSubCoils:
            coil.NumSubCoils = self.NumSubCoils
        if self.SubCoils:
            for selfsc, subcoil in zip(self.SubCoils, coil.SubCoils, strict=False):
                selfsc.do_init(subcoil)
        if self.dR and self.dR > 0.0:
            coil.compute_SubCoils(psiGrid)


class SubShellIn(BaseModel):
    "Subshell data."

    Name: str | None
    Current: float = 0.0
    R: float = Field(validation_alias=AliasChoices("R", "X"))
    Z: float

    def do_init(self, sshell: SubShell) -> None:
        if self.Name:
            sshell.Name = self.Name
        sshell.Current = self.Current
        sshell.R = self.R
        sshell.Z = self.Z


class ShellIn(BaseModel):
    """(Super) Conducting shell input."""

    Name: str | None
    Enabled: bool | None = True
    NumSubShells: int | None
    SubShells: list[SubShellIn]

    def do_init(self, shell: Shell) -> None:
        if self.Name:
            shell.Name = self.Name
        if self.Enabled:
            shell.Enabled = self.Enabled
        if shell.NumSubShells != self.NumSubShells:
            raise ValueError(
                f"NumSubShells {self.NumSubShells} does not match"
                f" the number of subshells {shell.NumSubShells}"
            )
        for selfss, subshell in zip(self.SubShells, shell.SubShells, strict=False):
            selfss.do_init(subshell)


class MeasureIn(BaseModel):
    """Measurement input data."""

    Name: str | None
    # Enabled: bool | None = True
    Type: MeasTypeA

    def do_init(self, meas: Measure) -> None:
        if self.Name:
            meas.Name = self.Name
        # meas.Enabled = self.Enabled  # no such field
        meas.Type = self.Type  # this is already converted during validation


class MachineIn(BaseModel):
    # Complete input file has these fields
    MaxIterFixed: int | None = 0
    MaxIterFree: int | None = 50
    Name: str | None = "pyDipolEQ"
    Info: str | None = "DipolEQ Equilibrium"
    Oname: str | None = None
    Iname: str | None = None
    MGname: str = ""
    LHname: str = ""
    RSname: str = ""
    RestartStatus: bool | None = False
    RestartUnkns: bool | None = False
    LHGreenStatus: bool | None = False
    MGreenStatus: bool | None = False
    NumEqualEq: int | None = 0
    VacuumOnly: bool | None = False
    FitMeasurements: bool | None = False
    PsiGrid: PsiGridIn
    Plasma: PlasmaIn
    Coils: list[CoilIn]
    Limiters: list[LimiterIn]
    Separatrices: list[SeparatrixIn] | None = None
    Measures: list[MeasureIn] | None = None
    Shells: list[ShellIn] | None = None
    NumCoils: int | None = None
    NumShells: int | None = None
    NumLimiters: int | None = None
    NumSeps: int | None = None
    NumMeasures: int | None = None

    @model_validator(mode="after")
    def check_numbers(self) -> Self:
        """check that Num* fields match the number of items in the lists
        Args:
            Self: the MachineIn object
        Returns:
            Self: the MachineIn object, validated
        """
        if self.NumCoils:
            if self.NumCoils != len(self.Coils):
                raise ValueError(
                    f"NumCoils {self.NumCoils} does not match the"
                    f" number of coils {len(self.Coils)}"
                )
        else:
            self.NumCoils = len(self.Coils)
        if self.NumShells:
            numShells = len(self.Shells) if self.Shells else 0
            if self.NumShells != numShells:
                raise ValueError(
                    f"NumShells {self.NumShells} does not match"
                    f" the number of shells {numShells}"
                )
        else:
            self.NumShells = len(self.Shells) if self.Shells else 0
        if self.NumLimiters:
            if self.NumLimiters != len(self.Limiters):
                raise ValueError(
                    f"NumLimiters {self.NumLimiters} does not "
                    f"match the number of limiters "
                    f"{len(self.Limiters)}"
                )
        else:
            self.NumLimiters = len(self.Limiters) if self.Limiters else 0
        if self.NumSeps:
            numSeps = len(self.Separatrices) if self.Separatrices else 0
            if self.NumSeps != numSeps:
                raise ValueError(
                    f"NumSeps {self.NumSeps} does not match the"
                    f" number of separatrices {numSeps}"
                )
        else:
            self.NumSeps = len(self.Separatrices) if self.Separatrices else 0
        if self.NumMeasures:
            numMeasures = len(self.Measures) if self.Measures else 0
            if self.NumMeasures != numMeasures:
                raise ValueError(
                    f"NumMeasures {self.NumMeasures} does not "
                    f"match the number of measures {numMeasures}"
                )
        else:
            self.NumMeasures = len(self.Measures) if self.Measures else 0

        if not self.Oname:
            self.Oname = self.Name

        if not self.Iname:
            self.Iname = self.Name

        return self

    def initalize_machine(self, m: Machine) -> None:
        """Initializes a Machine object from the input data"""

        # names
        if self.Name:
            m.Name = self.Name
        if self.Info:
            m.Info = self.Info
        if self.Oname:
            m.Oname = self.Oname
        if self.Iname:
            m.Iname = self.Iname
        m.MGname = self.MGname
        m.LHname = self.LHname
        m.RSname = self.RSname

        # control parameters
        if self.RestartStatus is not None:
            m.RestartStatus = self.RestartStatus
        if self.RestartUnkns is not None:
            m.RestartUnkns = self.RestartUnkns
        if self.LHGreenStatus is not None:
            m.LHGreenStatus = self.LHGreenStatus
        if self.MGreenStatus is not None:
            m.MGreenStatus = self.MGreenStatus
        if self.NumEqualEq is not None:
            m.NumEqualEq = self.NumEqualEq
        if self.VacuumOnly is not None:
            m.VacuumOnly = self.VacuumOnly

        if self.MaxIterFixed is not None:
            m.MaxIterFixed = self.MaxIterFixed
        if self.MaxIterFree is not None:
            m.MaxIterFree = self.MaxIterFree

        # set numbers
        if self.NumCoils is not None:
            m.NumCoils = self.NumCoils
        if self.NumShells is not None:
            m.NumShells = self.NumShells
        if self.NumLimiters is not None:
            m.NumLimiters = self.NumLimiters
        if self.NumSeps is not None:
            m.NumSeps = self.NumSeps
        if self.NumMeasures is not None:
            m.NumMeasures = self.NumMeasures

        # init PsiGrid and Plasma
        self.PsiGrid.do_init(m.PsiGrid)
        self.Plasma.do_init(m.Plasma)
        # these need to be equal
        m.Plasma.Nsize = m.PsiGrid.Nsize

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

        if self.Separatrices:
            for i, sep in enumerate(self.Separatrices):
                if not m.Seps[i]:
                    m.Seps[i] = m.Seps.new_separatrix()
                sep.do_init(m.Seps[i])

        if self.Measures:
            for i, meas in enumerate(self.Measures):
                if not m.Measures[i]:
                    m.Measures[i] = m.Measures.new_meas(meas.Type)
                meas.do_init(m.Measures[i])

        if self.Shells:
            for i, shell in enumerate(self.Shells):
                shell.do_init(m.Shells[i])
