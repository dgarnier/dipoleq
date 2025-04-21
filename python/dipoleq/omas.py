"""
Map to IMAS data structures. Since IMAS is not yet available
under public license, use OMAS.
"""

import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Tuple

from omas import ODS  # type: ignore[import-untyped]

from ._version import __version__, __version_tuple__
from .omas_dipole_extras import add_inner_boundary_to_omas
from .post_process import Machine
from .mas import imas_input_params, fill_ds
from .ds import DS, OmasDS

# export (or re-export) these functions
__all__ = [
    "to_omas",
    "ODS",
    "load_omas_data_structure",
    "omas_input_params",
]


def omas_input_params(ods: ODS) -> dict[str, Any] | None:
    return imas_input_params(ods["equilibrium"])


def load_omas_data_structure(filename: str | Path) -> ODS:
    """Load an OMAS data structure from a file"""
    add_inner_boundary_to_omas()
    ods = ODS()
    ods.load(str(filename))
    return ods


def prepare_omas_ds(ods: ODS | None) -> Tuple[ODS, DS, DS]:
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
