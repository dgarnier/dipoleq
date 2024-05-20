"""DipolEq: Levitated Dipole Equilibrium Solver

DipolEq is a Python package for solving the equilibrium of a levitated dipole
magnetized plasma. It is a Python wrapper around the C code that was originally
written by Michael Mauel and his group at Columbia University in the early 1999s
called TokaMac. The C code was later modified by Darren Garnier at Columbia
for use in the LDX experiment. The Python wrapper was written by Darren Garnier
at OpenStar Technologies, LTD.
    
"""

from ._c import Machine, Plasma, PsiGrid, CPlasmaModel, \
                Measure, Coil, Shell, Separatrix, Limiter, \
                SubCoil, SubShell, ModelType, MeasType, CircleType,\
                VectorView, MatrixView, IMatrixView
#from ._c import *

from .solve import solve
from .input import MachineIn
from .read_dotin import read_dotin

__all__ = [
    'Machine', 'Plasma', 'PsiGrid', 'CPlasmaModel',
    'Measure', 'Coil', 'Shell', 'Separatrix', 'Limiter',
    'SubCoil', 'SubShell', 'ModelType', 'MeasType', 'CircleType',
    'VectorView', 'MatrixView', 'IMatrixView',
    'solve', 'MachineIn', 'read_dotin'
]
