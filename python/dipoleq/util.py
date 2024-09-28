"""Utility functions for the dipoleq package"""

try:
    from enum import StrEnum
except ImportError:  # python < 3.11
    from enum import Enum

    class StrEnum(str, Enum):  # type: ignore[no-redef]
        pass


from collections.abc import Callable
from typing import Annotated, Any, Literal, ParamSpec, TypeVar

import numpy as np
from numpy.typing import NDArray

# lets pirate onto the c-class with extra functions
# just have to import this file... nice.
CT = TypeVar("CT")
T = TypeVar("T")
P = ParamSpec("P")


def add_method(cls: type[CT]) -> Callable[[Callable[P, T]], Callable[P, T]]:
    """Add a method to a class, used as a decorator to a function.
    The first argument to the method must be the class instance.
    """

    def decorator(func: Callable[P, T]) -> Callable[P, T]:
        setattr(cls, func.__name__, func)
        return func

    return decorator


# this typing business is a bit of a mess
# I think a package like "jaxtyping" would be better
# https://docs.kidger.site/jaxtyping/
ArrayN2 = Annotated[NDArray[np.generic], Literal["N", 2]]
ArrayN22 = Annotated[NDArray[np.generic], Literal["N", 2, 2]]


def segments_to_polygon(segments: ArrayN22) -> ArrayN2:
    """segments is a 3D array with shape (nsegs, 2, 2)
    with a start and end point for each segment
    """
    pts = [last_pt := np.empty(2) * np.nan] and [
        last_pt := pt
        for seg in segments  # noqa: F841, RUF100
        for pt in seg
        if not np.allclose(pt, last_pt)
    ]
    return np.array(pts)


def is_polygon_closed(polygon: ArrayN2) -> bool:
    """Check if a polygon is closed"""
    return np.allclose(polygon[0], polygon[-1])


def is_polygon_clockwise(polygon: ArrayN2) -> bool:
    """Check if a polygon is defined clockwise"""
    if polygon.shape[0] < 3:
        return True  # not enough points to determine
    # first find point on convex hull
    # find the min x, if there are more than one, take the min y
    if is_polygon_closed(polygon):  # remove last point if repeated
        polygon = polygon[:-1]
    # where is weird, puts the result in a tuple for some reason
    idxs = np.where(polygon[:, 0] == np.min(polygon[:, 0]))[0]
    # now find the lowest y index and convert to index in original
    idx = idxs[np.argmin(polygon[idxs, 1])]
    # make the convex hull point the second point
    rpoly = np.roll(polygon, -(idx - 1), axis=0)
    # now take the points around it and check the determinant
    # of the orientation matrix
    orientm = np.hstack([np.ones([3, 1]), rpoly[0:3]])
    return bool(np.linalg.det(orientm) < 0)


def _props(x: Any) -> dict[str, Any]:
    """Gets just the properties of an object"""
    return {
        key: getattr(x, key)
        for key in dir(x)
        if not key.startswith("__") and not callable(getattr(x, key))
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
            case "MachineIn":
                input2 = getattr(obj2, key, None)  # check if input_data exists
                if input2 is not None:
                    areDiff |= machine_diff(
                        value.model_dump(mode="json"),
                        getattr(obj2, key).model_dump(mode="json"),
                        root=f"{root}.{key}",
                    )

            case _:
                if value != getattr(obj2, key):
                    areDiff = True
                    print(f"{root}.{key}: {value} != {getattr(obj2, key)}")
    return areDiff
