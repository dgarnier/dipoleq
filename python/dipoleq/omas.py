"""
Map to OMAS data structures, a open-source implementation of IMAS.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from omas import ODS  # type: ignore[import-untyped]

from .mas import DS, fill_ds, mas_input_params
from .omas_dipole_extras import add_inner_boundary_to_omas
from .post_process import Machine

# export (or re-export) these functions
__all__ = [
    "to_omas",
    "ODS",
    "load_omas_data_structure",
    "omas_input_params",
]


@dataclass
class OmasDS(DS):
    _ods: ODS

    def _wrap_object(self, o: Any, force: bool) -> Any:
        if isinstance(o, ODS):
            return OmasDS(o)
        if force:
            raise TypeError(f"Unsupported type {type(o)} of {o} in force wrap mode")
        return o

    def _getitem(self, key: list[str | int], force: bool = False) -> Any:
        return self._wrap_object(self._ods[self._join_key(key)], force)

    def _setitem(self, key: list[str | int], value: Any) -> None:
        assert len(key) > 0
        if len(key) == 2 and isinstance(key[1], int):
            # Deal with time-dependent arrays
            orig_value = []
            if key[0] in self._ods:
                orig_value = np.atleast_1d(self._ods[key[0]]).tolist()
            if key[1] < len(orig_value):
                orig_value[key[1]] = value
            elif key[1] == len(orig_value):
                orig_value += [value]
            else:
                raise IndexError
            self._ods[key[0]] = np.atleast_1d(orig_value)
        elif len(key) == 1:
            self._ods[key[0]] = value
        else:
            self._getitem([key[0]], True)._setitem(key[1:], value)

    def _join_key(self, key: list[str | int]) -> str:
        assert all(isinstance(k, str | int) for k in key)
        return ".".join([str(k) for k in key])

    def __len__(self) -> int:
        return len(self._ods)

    @property
    def inner(self) -> ODS:
        return self._ods


def omas_input_params(ods: ODS) -> dict[str, Any] | None:
    """Create a DipolEQ input data from an OMAS equilibrium data structure"""
    return mas_input_params(OmasDS(ods["equilibrium"]))


def load_omas_data_structure(filename: str | Path) -> ODS:
    """Load an OMAS data structure from a file"""
    add_inner_boundary_to_omas()
    ods = ODS()
    ods.load(str(filename))
    return ods


def prepare_omas_ds(ods: ODS | None) -> tuple[ODS, DS, DS]:
    if ods is None:
        # need this for consistency check to pass
        # extend the data dictionary for the inner boundary
        add_inner_boundary_to_omas()
        # by default, consistency check is on, and cocos=11
        ods = ODS(cocos=11)
    eq = ods["equilibrium"]
    wall = ods["wall"]
    return ods, OmasDS(eq), OmasDS(wall)


def to_omas(
    m: Machine, ods: ODS | None = None, time_index: int | None = None, time: float = 0.0
) -> ODS:
    """Add the equilibrium data to an OMAS data structure."""
    ods, eq, wall = prepare_omas_ds(ods)

    fill_ds(m, eq, wall, time_index, time)

    return ods
