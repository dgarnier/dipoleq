"""Utility functions for the dipoleq package"""

try:
    from enum import StrEnum
except ImportError:  # python < 3.11
    from enum import Enum

    class StrEnum(str, Enum):  # type: ignore[no-redef]
        pass


from typing import Any

import numpy as np

# from .core import (
#    Machine, Plasma, PsiGrid, CPlasmaModel,
#    Measure, Coil, Shell, Separatrix, Limiter, SubCoil, SubShell,
#    Measures, Coils, Shells, Separatrices, Limiters, SubCoils, SubShells,
#    ModelType, MeasType, CircleType, VectorView, MatrixView, IMatrixView
# )


def _props(x: Any) -> dict[str, Any]:
    """Gets just the properties of an object"""
    return {
        key: getattr(x, key)
        for key in dir(x)
        if not callable(getattr(x, key)) and not key.startswith("__")
    }


def machine_diff(obj1: Any, obj2: Any, root: str = "", verbose: bool = False) -> bool:
    """Print the differences between two objects
    Args:
        obj1 (object): The first object
        obj2 (object): The second object
        root (str): The root of the object
    """
    areDiff = False
    p1 = _props(obj1)  # filter out methods and class attributes
    if verbose:
        print(f"Comparing {root}.*")
    for key, value in p1.items():
        if verbose:
            print(f"Comparing {root}.{key}, class: {value.__class__.__name__}")
        match value.__class__.__name__:
            case (
                "Coils"
                | "Shells"
                | "Separatrices"
                | "Measures"
                | "Limiters"
                | "SubCoils"
                | "SubShells"
            ):
                for i, coil in enumerate(value):
                    areDiff |= machine_diff(
                        coil, getattr(obj2, key)[i], root=f"{root}.{key}[{i}]"
                    )
            case (
                "Coil"
                | "Shell"
                | "Separatrix"
                | "Measure"
                | "Limiter"
                | "SubCoil"
                | "SubShell"
                | "Plasma"
                | "PsiGrid"
                | "Machine"
            ):
                areDiff |= machine_diff(value, getattr(obj2, key), root=f"{root}.{key}")
            case "CPlasmaModel":  # not sure how to handle this yet
                pass
            case "VectorView" | "MatrixView" | "IMatrixView":
                s1, s2 = (np.array(value).shape, np.array(getattr(obj2, key)).shape)
                if s1 != s2:
                    print(f"{root}.{key}: Array shapes differ {s1}!={s2}")
                elif not np.array_equal(value, getattr(obj2, key)):
                    print(f"{root}.{key}: Arrays not equal")
            case "float":
                if not np.isclose(value, getattr(obj2, key)):
                    print(f"{root}.{key}: {value} != {getattr(obj2, key)}")
            case _:
                if value != getattr(obj2, key):
                    print(f"{root}.{key}: {value} != {getattr(obj2, key)}")
    return areDiff
