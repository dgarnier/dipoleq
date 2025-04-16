"""
"""

import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import numpy as np
from imas import DBEntry, ids_defs
from imas.ids_toplevel import IDSToplevel
from imas.ids_structure import IDSStructure
from json2xml.json2xml import Json2xml

from ._version import __version__, __version_tuple__
from .input import MachineIn
from .post_process import Machine
from .util import is_polygon_closed
from .mas import add_boundary, add_mas_code_info, add_boundary_separatrix, add_limiters

# export (or re-export) these functions
__all__ = [
    "to_imas"
]
# TODO: Make sure we are inputting with COCOS 11


def to_imas(
    m: Machine, db: DBEntry, time_index: int | None = None, time: float = 0.0
) -> None:
    """Add the equilibrium data to an OMAS data structure."""
    pl = m.Plasma
    pg = m.PsiGrid
    equilibrium_ids = db.factory.equilibrium()
    equilibrium_ids.ids_properties.homogeneous_time = ids_defs.IDS_TIME_MODE_HOMOGENEOUS
    wall_ids = db.factory.wall()
    wall_ids.ids_properties.homogeneous_time = ids_defs.IDS_TIME_MODE_HOMOGENEOUS

    add_limiters(m, wall_ids, True)

    input_data = getattr(m, "input_data", None)
    add_mas_code_info(equilibrium_ids, input_data)

    equilibrium_ids['vacuum_toroidal_field/r0'] = pl.R0

    if time_index is None:
        time_index = len(equilibrium_ids['time_slice'])
        equilibrium_ids['time_slice'].resize(time_index+1)
    eqt = equilibrium_ids['time_slice'][time_index]

    equilibrium_ids.time_slice[time_index].time = time
    equilibrium_ids.time = [0.0]

    # the values here are taken from the IMAS schema
    # TODO: Link

    psi = np.array(pl.Psi_pr)  # flux values, 1d
    coordsio = {f"equilibrium.time_slice.{time_index}.profiles_1d.psi": psi}

    # with omas_environment(ods, cocosio=11, coordsio=coordsio):
    if True:
        # Set the time array
        # eqt["time"] = time
        # ods.set_time_array("equilibrium.time", time_index, eqt["time"])

        # 0D quantities
        glob = eqt["global_quantities"]
        glob["psi_axis"] = pg.PsiMagAxis
        glob["psi_boundary"] = pg.PsiLim
        # glob["psi_inner_boundary"] = pg.PsiAxis  # TODO: Put this back in
        glob["magnetic_axis"]["r"] = pg.RMagAxis
        glob["magnetic_axis"]["z"] = pg.ZMagAxis
        glob["ip"] = pl.Ip

        # B0, R0 is weird
        # ods["equilibrium.vacuum_toroidal_field/r0"] = pl.R0
        # ods.set_time_array("equilibrium.vacuum_toroidal_field/b0", time_index, pl.B0)
        # eqt['vacuum_toroidal_field/r0'] = pl.R0

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
    
    equilibrium_ids.validate()
    wall_ids.validate()

    db.put(equilibrium_ids)
    db.put(wall_ids)
