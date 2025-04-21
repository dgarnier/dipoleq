"""
"""

from pathlib import Path

from imas.dd_zip import Traversable, lru_cache, Dict, Tuple, Union, _open_zipfile, ZIPFILE_LOCATIONS, re
@lru_cache
def _read_dd_versions() -> Dict[str, Tuple[Union[Path, Traversable], str]]:
    """Traverse all possible DD zip files and return a map of known versions.

    Returns:
        version_map: version -> (zipfile path, filename)
    """
    versions = {}
    xml_re = re.compile(r"^data-dictionary/(\S+)\.xml$")
    for path in ZIPFILE_LOCATIONS:
        if not path.is_file():
            continue
        with _open_zipfile(path) as zipfile:
            for fname in zipfile.namelist():
                match = xml_re.match(fname)
                if match:
                    version = match.group(1)
                    if version not in versions:
                        versions[version] = (path, fname)
    if not versions:
        raise RuntimeError(
            "Could not find any data dictionary definitions. "
            f"Looked in: {', '.join(map(repr, ZIPFILE_LOCATIONS))}."
        )
    return versions
import imas
imas.dd_zip._read_dd_versions = _read_dd_versions

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
