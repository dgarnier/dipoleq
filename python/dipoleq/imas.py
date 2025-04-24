"""
Map to IMAS data structures, using the IMAS-Python library.
"""

from dataclasses import dataclass
from typing import Any

import dipmas  # noqa: F401  # pylint: disable=unused-import
from imas import DBEntry, IDSFactory, convert_ids, ids_defs
from imas.ids_base import IDSBase
from imas.ids_primitive import IDSPrimitive

from .mas import DS, fill_ds, mas_input_params
from .post_process import Machine

# export (or re-export) these functions
__all__ = [
    "to_imas",
    "imas_input_params",
    "DBEntry",
    "ImasDS",
]


@dataclass
class ImasDS(DS):
    _ids: IDSBase

    def _wrap_object(self, o: Any, force: bool) -> Any:
        if isinstance(o, IDSPrimitive) and not force:
            return o.value
        if isinstance(o, IDSBase):
            return ImasDS(o)
        if force:
            raise TypeError(f"Unsupported type {type(o)} of {o} in force wrap mode")
        return o

    def _getitem(self, key: list[str | int], force: bool = False) -> Any:
        if len(key) == 1:
            if isinstance(key[0], int) and key[0] == len(self._ids):
                self._ids.resize(len(self._ids) + 1)
            return self._wrap_object(self._ids[key[0]], force)
        return self._getitem([key[0]])._getitem(key[1:])

    def _setitem(self, key: list[str | int], value: Any) -> None:
        if len(key) == 1:
            if isinstance(key[0], int) and key[0] == len(self._ids):
                self._ids.resize(len(self._ids) + 1)
            self._ids[key[0]] = value
        else:
            self._getitem([key[0]], True)._setitem(key[1:], value)

    def _join_key(self, key: list[str | int]) -> str:
        result = ""
        for i, k in enumerate(key):
            if isinstance(k, int):
                result += f"[{k}]"
            elif isinstance(k, str):
                if i != 0:
                    result += "/"
                result += k
            else:
                raise ValueError("Invalid key type")
        return result

    def __len__(self) -> int:
        return len(self._ids)

    @property
    def inner(self) -> IDSBase:
        return self._ids


def prepare_imas_ds() -> tuple[ImasDS, ImasDS]:
    """Create the equilibrium and wall IDSs"""
    # OMAS uses 3.41.0 (COCOS 11), so the DipolEQ-IMAS interface must also use COCOS 11
    factory = IDSFactory("3.41.0+dipole")
    eq = factory.equilibrium()
    eq.ids_properties.homogeneous_time = ids_defs.IDS_TIME_MODE_HOMOGENEOUS
    wall = factory.wall()
    wall.ids_properties.homogeneous_time = ids_defs.IDS_TIME_MODE_HOMOGENEOUS
    return ImasDS(eq), ImasDS(wall)


def imas_input_params(equilibrium: IDSBase) -> dict[str, Any] | None:
    """Create a DipolEQ input data from an IMAS equilibrium IDS"""
    return mas_input_params(ImasDS(equilibrium))


def to_imas(
    m: Machine, db: DBEntry, time_index: int | None = None, time: float = 0.0
) -> DBEntry:
    """Add the equilibrium data to an IMAS database"""
    eq, wall = prepare_imas_ds()

    fill_ds(m, eq, wall, time_index, time)

    eq = convert_ids(eq.inner, None, factory=db.factory)
    wall = convert_ids(wall.inner, None, factory=db.factory)

    eq.validate()
    wall.validate()

    db.put(eq)
    db.put(wall)

    return db
