"""
"""

import dipmas
from imas import DBEntry, ids_defs

from .post_process import Machine
from .mas import fill_ds
from .ds import ImasDS

# export (or re-export) these functions
__all__ = [
    "to_imas"
]
# TODO: Make sure we are inputting with COCOS 11

def prepare_imas_ds(db: DBEntry):
    eq = db.factory.equilibrium()
    eq.ids_properties.homogeneous_time = ids_defs.IDS_TIME_MODE_HOMOGENEOUS
    wall = db.factory.wall()
    wall.ids_properties.homogeneous_time = ids_defs.IDS_TIME_MODE_HOMOGENEOUS
    return ImasDS(eq), ImasDS(wall)


def to_imas(
    m: Machine, db: DBEntry, time_index: int | None = None, time: float = 0.0
) -> DBEntry:
    """Add the equilibrium data to an OMAS data structure."""
    eq, wall = prepare_imas_ds(db)

    fill_ds(m, eq, wall, time_index, time)

    eq.inner.validate()
    wall.inner.validate()

    db.put(eq.inner)
    db.put(wall.inner)

    return db
