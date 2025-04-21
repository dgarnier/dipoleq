"""
Map to IMAS data structures. Since IMAS is not yet available
under public license, use OMAS.
"""

import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import numpy as np
from json2xml.json2xml import Json2xml  # type: ignore[import-untyped]
from omas import ODS, omas_environment  # type: ignore[import-untyped]

# from .core import Machine
from ._version import __version__, __version_tuple__
from .omas_dipole_extras import add_inner_boundary_to_omas
from .input import MachineIn
from .post_process import Machine
from .util import is_polygon_closed
from .mas import add_boundary, add_mas_code_info, add_boundary_separatrix, add_limiters, add_inner_boundary_separatrix, mas_input_params

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


def prepare_ods(ods: ODS | None):
    if ods is None:
        # need this for consistency check to pass
        # extend the data dictionary for the inner boundary
        add_inner_boundary_to_omas()
        # by default, consistency check is on, and cocos=11
        ods = ODS(cocos=11)
    eq = ods["equilibrium"]
    wall = ods["wall"]
    return ods, eq, wall


def to_omas(
    m: Machine, ods: ODS | None = None, time_index: int | None = None, time: float = 0.0
) -> ODS:
    """Add the equilibrium data to an OMAS data structure."""
    pl = m.Plasma
    pg = m.PsiGrid
    ods, eq, wall = prepare_ods(ods)

    add_limiters(m, wall)

        # add the structure and the wall from the machine
    input_data = getattr(m, "input_data", None)
    add_imas_code_info(eq, input_data=input_data)

    if time_index is None:
        time_index = len(eq["time_slice"])
    eqt = eq["time_slice"][time_index]

    # the values here are taken from the IMAS/OMAS schema
    # https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/equilibrium.html
    # https://gafusion.github.io/omas/schema/schema_equilibrium.html

    psi = np.array(pl.Psi_pr)  # flux values, 1d

    # Set the time array
    eqt["time"] = time
    ods.set_time_array("equilibrium.time", time_index, time)

    # 0D quantities
    glob = eqt["global_quantities"]
    glob["psi_axis"] = pg.PsiMagAxis
    glob["psi_boundary"] = pg.PsiLim
    glob["psi_inner_boundary"] = pg.PsiAxis
    glob["magnetic_axis.r"] = pg.RMagAxis
    glob["magnetic_axis.z"] = pg.ZMagAxis
    glob["ip"] = pl.Ip

    # B0, R0 is weird
    eq["vacuum_toroidal_field.r0"] = pl.R0
    ods.set_time_array("equilibrium.vacuum_toroidal_field.b0", time_index, pl.B0)

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
    eq2d = eqt["profiles_2d.0"]
    eq2d["type.index"] = (
        0  # total fields.. could also be broken down into components
    )
    eq2d["grid_type.index"] = 1  # regular R,Z grid
    eq2d["grid_type.name"] = "RZ"
    eq2d["grid.dim1"] = R = np.array(m.PsiGrid.R)
    eq2d["grid.dim2"] = np.array(m.PsiGrid.Z)
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
    add_boundary_separatrix(m, eqt)
    add_inner_boundary_separatrix(m, eqt)

    return ods
