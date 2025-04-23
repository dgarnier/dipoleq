import xml.etree.ElementTree as ET
from abc import ABC, abstractmethod
from typing import Any

import numpy as np
from json2xml.json2xml import Json2xml  # type: ignore[import-untyped]

from ._version import __version__, __version_tuple__
from .input import MachineIn
from .post_process import Machine
from .util import is_polygon_closed


class DS(ABC):
    def __getitem__(self, key: str | int) -> Any:
        return self._getitem(self._separate_key(key))

    def __setitem__(self, key: str | int, value: Any) -> None:
        self._setitem(self._separate_key(key), value)

    def __getattr__(self, name: str) -> Any:
        if name.startswith("_") and hasattr(super(), "__getattr__"):
            return super().__getattr__(name)  # type: ignore [misc]
        return self[name]

    def __setattr__(self, name: str, value: Any) -> None:
        if name.startswith("_"):
            super().__setattr__(name, value)
        else:
            self[name] = value

    def _separate_key(self, key: str | int) -> list[str | int]:
        if isinstance(key, int):
            return [key]
        split_key = key.replace("/", ".").replace("]", "").replace("[", ".").split(".")
        return [int(k) if k.isdecimal() else k for k in split_key]

    @abstractmethod
    def _wrap_object(self, o: Any, force: bool) -> Any: ...

    @abstractmethod
    def _getitem(self, key: list[str | int]) -> Any: ...

    @abstractmethod
    def _setitem(self, key: list[str | int], value: Any) -> None: ...

    @abstractmethod
    def __len__(self) -> int: ...

    @property
    @abstractmethod
    def inner(self) -> Any: ...


def add_imas_code_info(ds: DS, input_data: MachineIn | None = None) -> None:
    """Add the code info to an OMAS data structure"""
    code = ds["code"]
    code["name"] = "DipolEQ"
    code["description"] = "Dipole equilibrium solver"
    code["repository"] = "https://github.com/dgarnier/dipoleq"
    # return tag of release or commit hash depending on if its a release or not
    code["version"] = (
        f"v{__version__}"
        if len(__version_tuple__) < 5
        else str(__version_tuple__[4]).split(".", maxsplit=1)[0]
    )
    # could add the input data as XML, but... ugh
    # needs pydantic_xml, xmltodict doesn't preserve types going roundtrip.
    # but pydantic works.. so lets just save as json.
    if input_data:
        input_dict = input_data.model_dump(mode="json", exclude_unset=True)
        xml = Json2xml(input_dict, wrapper="Machine", pretty=True).to_xml()
        code["parameters"] = xml

    # input_data.model_dump(mode=)
    # could add which libraries and commits to add, but, really
    # this is enough to reproduce the results
    # one would encode the input data as XML into a string
    # code['output_flag'][eq_time_index] = 0 if run is successful,
    # negative if not to be used


def mas_input_params(equilibrium: DS) -> dict[str, Any] | None:
    """Create a DipolEQ input data from an equilibrium data structure"""
    # need to create the MachineIn object from the input data
    # which is stored in the code section of the data structure.
    # this is a bit of a pain, but it is what it is.
    # the input data is stored as XML in the code section
    # so we need to extract it and then load it into the MachineIn object.
    if equilibrium["code.name"] != "DipolEQ":
        return None

    def xml_to_dict(node: ET.Element) -> Any:
        """Compatible conversion from XML based on Json2xml output"""
        match node.attrib.get("type", "dict"):
            case "dict":
                return {child.tag: xml_to_dict(child) for child in node}
            case "list":
                return [xml_to_dict(child) for child in node]
            case "str":
                return node.text
            case "int":
                return int(node.text) if node.text is not None else None
            case "float":
                return float(node.text) if node.text is not None else None
            case "bool":
                return node.text and node.text.lower() in ("true", "yes", "t", "1")

    data = xml_to_dict(ET.fromstring(equilibrium["code.parameters"]))
    return data if isinstance(data, dict) else None


def add_limiters(m: Machine, wall: DS) -> None:
    """Add the limiter as a IMAS/OMAS wall structure
    For limiters, there are different types.
    Each type can have a name, index, and description.
    index = 0 = single contour
            1 = PFC structure block
    the limiter is also made of units.. which can
    also have component types.
    """
    # add data providence to wall structure
    add_imas_code_info(wall)

    # IMAS units should be contiguous and clockwise
    # of course, the inner limiter will be done in reverse
    # since it limits on the outside of its surface.
    wall0 = wall["description_2d[0]"]
    wall0["type.name"] = "Limiter(s)"
    wall0["type.index"] = 1  # multiple limiter units, no vessel
    wall0["type.description"] = (
        "As a dipole, this will contain outer and inner limiters"
    )
    wall0["limiter.type.index"] = 0  # official single contour limiter
    wall0["limiter.type.name"] = "Limiter Contour(s)"

    # outer limiter first
    outline = m.Limiters.olim_outline()
    unit = wall0["limiter.unit[0]"]
    unit["name"] = "Outer limiter"
    unit["identifier"] = "outer_limiter"
    if is_polygon_closed(outline):  # determine closed by last
        # unit['outline.closed'] = 1  .. not in IMAS schema
        outline = outline[:-1]  # remove last point for closed
    # else:
    #    unit['outline.closed'] = 0
    unit["outline.r"] = outline[:, 0]
    unit["outline.z"] = outline[:, 1]
    unit["component_type.name"] = "outer_limiter"
    unit["component_type.index"] = 5  # 5 = limiter

    # inner limiter
    outline = m.Limiters.ilim_outline()
    unit = wall0["limiter.unit[1]"]
    unit["name"] = "inner_limiter"
    unit["identifier"] = "inner_limiter"
    if is_polygon_closed(outline):  # determine closed by last
        #    unit['outline.closed'] = 1
        outline = outline[:-1]
    # else:
    #    unit['outline.closed'] = 0
    unit["outline.r"] = outline[:, 0]
    unit["outline.z"] = outline[:, 1]
    unit["component_type.name"] = "inner_limiter"
    unit["component_type.index"] = -5  # 5 = limiter, but "private/custom" type


def add_boundary(m: Machine, ts: DS) -> None:
    """Add the boundary data to an IMAS/OMAS equilibrium time slice"""
    # the values here are taken from the IMAS/OMAS schema
    # https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/equilibrium.html
    # https://gafusion.github.io/omas/schema/schema_equilibrium.html
    # make boundary not quite at the separatrix
    # this is the 99.5% flux surface
    pg = m.PsiGrid
    psi_norm = m.Plasma.PsiXmax  # set by input
    bound = ts["boundary"]
    bnd_r, bnd_z = pg.get_contour(psi_norm)
    bound["outline.r"] = bnd_r
    bound["outline.z"] = bnd_z
    bound["psi_norm"] = psi_norm
    bound["psi"] = (pg.PsiLim - pg.PsiAxis) * psi_norm + pg.PsiAxis


def add_boundary_separatrix(m: Machine, ts: DS) -> None:
    """Add the boundary separatrix data to an IMAS/OMAS equilibrium time slice"""
    # the values here are taken from the IMAS/OMAS schema
    # https://imas-data-dictionary.readthedocs.io/en/latest/generated/ids/equilibrium.html
    # https://gafusion.github.io/omas/schema/schema_equilibrium.html
    sep_r, sep_z = m.PsiGrid.get_contour(1.0)
    bsep = ts["boundary_separatrix"]
    bsep["outline.r"] = sep_r
    bsep["outline.z"] = sep_z
    bsep["psi"] = m.PsiGrid.PsiLim
    if m.is_diverted():  # type: ignore[attr-defined]
        bsep["type"] = 1
    else:
        bsep["type"] = 0
        active_point = m.get_outer_limiter_contact_point()  # type: ignore[attr-defined]
        if active_point is not None:
            bsep["active_limiter_point.r"] = active_point[0]
            bsep["active_limiter_point.z"] = active_point[1]
    # not sure how useful this is if not diverted, but
    # decided to add them anyway.. would be more helpful
    # with the flux values there.
    for i, xp in enumerate(m.get_x_points()):  # type: ignore[attr-defined]
        bsep[f"x_point[{i}].r"] = xp.Rs
        bsep[f"x_point[{i}].z"] = xp.Zs
        # this doesn't exist in IMAS schema
        # bsep[f'x_point[{i}].psi'] = xp.Psi

    # FIXME:  other things to calculate..
    # closest_wall_point
    # dr_dz_zero_point, geometric_axis
    # elongation, elongation_upper, elongation_lower
    # triangularity, triangularity_upper, triangularity_lower
    # triangularity_inner, triangularity_outer
    # squareness (in alpha)
    # strikepoint(s), gap(s)


def add_inner_boundary_separatrix(m: Machine, ts: DS) -> None:
    """Add the inner boundary separatrix data to an equilibrium time slice"""
    if m.PsiGrid.PsiAxis == m.PsiGrid.PsiMagAxis:
        # no fcfs when the inner boundary is at the axis.
        return
    sep_r, sep_z = m.PsiGrid.get_contour(0.0)
    bsep = ts["inner_boundary_separatrix"]
    bsep["outline.r"] = sep_r
    bsep["outline.z"] = sep_z
    bsep["psi"] = m.PsiGrid.PsiAxis
    bsep["type"] = 0
    active_point = m.get_inner_limiter_contact_point()  # type: ignore[attr-defined]
    if active_point is not None:
        bsep["active_limiter_point.r"] = active_point[0]
        bsep["active_limiter_point.z"] = active_point[1]


def fill_ds(m: Machine, eq: DS, wall: DS, time_index: int | None, time: float) -> None:
    """Add all the equilibrim and wall information from a DipolEQ Machine
    into the equilibrium and wall data structures.
    """
    pl = m.Plasma
    pg = m.PsiGrid

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
    eq[f"time.{time_index}"] = time

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
    eq[f"vacuum_toroidal_field.b0.{time_index}"] = pl.B0

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
    eq2d["type.index"] = 0  # total fields.. could also be broken down into components
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
