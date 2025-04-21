"""
"""

import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import numpy as np

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
from imas.ids_toplevel import IDSToplevel
from imas.ids_structure import IDSStructure
from json2xml.json2xml import Json2xml

from ._version import __version__, __version_tuple__
from .input import MachineIn
from .post_process import Machine
from .util import is_polygon_closed
from .mas import add_boundary, add_mas_code_info, add_boundary_separatrix, add_limiters, get_or_append, set_or_append

# export (or re-export) these functions
__all__ = [
    "to_imas"
]
# TODO: Make sure we are inputting with COCOS 11

def prepare_imas_things(db: DBEntry):
    eq = db.factory.equilibrium()
    eq.ids_properties.homogeneous_time = ids_defs.IDS_TIME_MODE_HOMOGENEOUS
    wall = db.factory.wall()
    wall.ids_properties.homogeneous_time = ids_defs.IDS_TIME_MODE_HOMOGENEOUS
    return eq, wall


def to_imas(
    m: Machine, db: DBEntry, time_index: int | None = None, time: float = 0.0
) -> DBEntry:
    """Add the equilibrium data to an OMAS data structure."""
    pl = m.Plasma
    pg = m.PsiGrid
    eq, wall = prepare_imas_things(db)

    add_limiters(m, wall, True)

    # add the structure and the wall from the machine
    input_data = getattr(m, "input_data", None)
    add_mas_code_info(eq, input_data=input_data)

    eq["vacuum_toroidal_field"]["r0"] = pl.R0

    if time_index is None:
        # Append to the end
        time_index = len(eq["time_slice"])
    eqt = get_or_append(eq["time_slice"], time_index)

    # the values here are taken from the IMAS/OMAS schema
    # https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/equilibrium.html
    # https://gafusion.github.io/omas/schema/schema_equilibrium.html

    psi = np.array(pl.Psi_pr)  # flux values, 1d

    # Set the time array
    eqt["time"] = time
    set_or_append(eq["time"], time_index, time)

    # 0D quantities
    glob = eqt["global_quantities"]
    glob["psi_axis"] = pg.PsiMagAxis
    glob["psi_boundary"] = pg.PsiLim
    glob["psi_inner_boundary"] = pg.PsiAxis
    glob["magnetic_axis"]["r"] = pg.RMagAxis
    glob["magnetic_axis"]["z"] = pg.ZMagAxis
    glob["ip"] = pl.Ip

    # B0, R0 is weird
    # set_or_append(eq["vacuum_toroidal_field"]["b0"], time_index, pl.B0)

    # 1D quantities
    eq1d = eqt["profiles_1d"]
    eq1d["psi"] = psi
    eq1d["f"] = np.asarray(pl.G_pr) * pl.B0R0
    eq1d["f_df_dpsi"] = np.asarray(pl.G2p_pr) * (pl.B0R0) ** 2
    eq1d["pressure"] = np.array(pl.P_pr)
    eq1d["dpressure_dpsi"] = np.array(pl.Pp_pr)
    eq1d["q"] = np.array(pl.q_pr)

    # 2D quantities
    MU0 = 4.0e-7 * 3.14159265358979323846
    eqt.profiles_2d.resize(1)
    eq2d = eqt["profiles_2d[0]"]
    eq2d["type"]["index"] = (
        0  # total fields.. could also be broken down into components
    )
    eq2d["grid_type"]["index"] = 1  # regular R,Z grid
    eq2d["grid_type"]["name"] = "RZ"
    eq2d["grid"]["dim1"] = R = np.array(m.PsiGrid.R)
    eq2d["grid"]["dim2"] = np.array(m.PsiGrid.Z)
    eq2d["psi"] = np.array(m.PsiGrid.Psi)
    eq2d["j_tor"] = np.asarray(m.PsiGrid.Current) / MU0
    eq2d["b_field_r"] = np.asarray(pl.GradPsiZ) / (2 * np.pi * R)
    eq2d["b_field_z"] = -np.asarray(pl.GradPsiR) / (2 * np.pi * R)
    eq2d["b_field_tor"] = np.array(pl.Bt)
    # others to add
    # eq2d['grid.volume_element']
    # eq2d['phi']   # the toroidal flux

    # boundaries
    add_boundary(m, eqt)
    add_boundary_separatrix(m, eqt, True)
    # add_inner_boundary_separatrix(m, eqt)  # TODO: Add this back in

    eq.validate()
    wall.validate()

    db.put(eq)
    db.put(wall)

    return db
