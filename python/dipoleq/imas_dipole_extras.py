"""
In order to add to IMAS data structures, and make sense for dipoles, we need to
add information that mostly has to do with the inner wall.

This means we need an inner wall object, and inner boundary and the profiles
should be defined only between the walls.

IMAS has different walls.. so its possible to use a second wall object, but this
is not the case for the boundary.

"""

import json
import logging
from pathlib import Path
from typing import Any

from omas.omas_setup import IMAS_versions  # type: ignore[import-untyped]
from omas.omas_utils import (  # type: ignore[import-untyped]
    _extra_structures,
    structures_filenames,
)

logger = logging.getLogger(__name__)
# we won't be adding new limiter or wall structures, but will be adding new
# custom types of limiters of walls with negative type indices to indicate.
# since dipole limiter and walls are not in the IMAS schema, the type index
# for these will be negative.

# that said, there is no such type index for the equilibrium boundary
#  "equilibrium.time_slice[:].boundary" is used for lcfs.

# rather than create new structures by hand, we will use the existing
# structure and copy it for a new inner_boundary structure.
# so take the "active" structure (not obsolescent or alpha) and copy it
# changing the names and descriptions as needed.


def get_imas_latest_version() -> str:
    """Get the latest IMAS version"""
    imas_versions = IMAS_versions()
    if len(imas_versions.keys()):
        return str(list(imas_versions.keys())[-1])
    return "3.41.0"  # latest in OMAS as of 2024-09-24


def _get_default_equilibrium_imas_struct_defs() -> dict[str, dict[str, Any]]:
    fn = structures_filenames(get_imas_latest_version())["equilibrium"]
    with Path(fn).open(mode="r", encoding="utf-8") as f:
        return json.load(f)  # type: ignore[no-any-return]


def _create_inner_boundary_structure() -> dict[str, dict[str, Any]]:
    # get the default boundary definitions
    equil_defs = _get_default_equilibrium_imas_struct_defs()
    # get only the ones we need
    # boundary is for the 99.x% of the outer boundary (from fixed boundary codes?)
    #   - do we want to use this for the inner boundary?
    #   - do we want just one?
    # boundary_separatrix is for the actual separatrix
    # boundary_secondary_separatrix.. 2nd x-point?  : NO for now.
    bound_defs = {
        k: v
        for k, v in equil_defs.items()
        if k.startswith(
            (
                "equilibrium.time_slice[:].boundary.",
                "equilibrium.time_slice[:].boundary_separatrix.",
                "equilibrium.time_slice[:].global_quantities.psi_boundary",
            )
        )
        and v["lifecycle_status"] == "active"
    }

    # iterate and change the keys, paths, and descriptions
    inner_boundary_extras = {}
    for key, value in bound_defs.items():
        k = key.replace("boundary", "inner_boundary")
        v = value.copy()
        v["documentation"] = v["documentation"].replace("boundary", "inner boundary")
        v["full_path"] = v["full_path"].replace("boundary", "inner_boundary")
        logger.debug("Adding to IMAS: %s", k)
        inner_boundary_extras[k] = v
    return inner_boundary_extras


def add_inner_boundary_to_imas() -> None:
    """Add the inner boundary to the IMAS equilibrium schema"""
    extra_equilibrium_structures = _extra_structures.get("equilibrium", {})
    extra_equilibrium_structures.update(_create_inner_boundary_structure())
    _extra_structures["equilibrium"] = extra_equilibrium_structures
