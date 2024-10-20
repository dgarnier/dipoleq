"""
Extra analysis routines for after a solve.

These are pure python functions that extend the machine object functionality.
"""

import numpy as np
import numpy.typing as npt
from scipy.integrate import cumulative_trapezoid

from .core import Limiters, Machine, Separatrix
from .util import (
    ArrayN2,
    add_method,
    area_of_polygon,
    is_polygon_clockwise,
    segments_to_polygon,
)

# will re-export extensions to these core classes
# this is only necessary if including from this file
# but its not necessary to do that since __init__.py will
__all__ = ["Machine", "Limiters", "Separatrix"]


@add_method(Machine)
def is_outer_limited(m: Machine) -> bool:
    """Check if the equilibrium is diverted"""
    return any(lim.Enabled > 0 and lim.PsiLim <= m.PsiGrid.PsiLim for lim in m.Limiters)


@add_method(Machine)
def is_diverted(m: Machine) -> bool:
    """Check if the equilibrium is diverted."""
    return not is_outer_limited(m)


@add_method(Machine)
def get_outer_limiter_contact_point(m: Machine) -> tuple[float, float] | None:
    """Get the outer limiter contact point"""
    for lim in m.Limiters:
        if lim.Enabled > 0 and np.isclose(lim.PsiLim, m.PsiGrid.PsiLim):
            return lim.RLim, lim.ZLim
    return None


@add_method(Machine)
def get_inner_limiter_contact_point(m: Machine) -> tuple[float, float] | None:
    """Get the inner limiter contact point"""
    for lim in m.Limiters:
        if lim.Enabled < 0 and np.isclose(lim.PsiLim, m.PsiGrid.PsiAxis):
            return lim.RLim, lim.ZLim
    return None


@add_method(Machine)
def get_x_points(m: Machine) -> list[Separatrix]:
    """Get the X-points (separatrices)"""
    valid_seps = [
        sep
        for sep in m.Seps
        if sep.Enabled and sep.IsSeparatrix and sep.Psi > m.PsiGrid.PsiAxis
    ]
    valid_seps.sort(key=lambda sep: sep.Psi)
    return valid_seps


@add_method(Separatrix)
def __repr__(sep: Separatrix) -> str:
    """Pretty print a separatrix"""
    if sep.IsSeparatrix and sep.Enabled:
        return (
            f"<X-point: {sep.Name} Psi: {sep.Psi:.3f}, "
            f"Rs: {sep.Rs:.3f}, Zs: {sep.Zs:.3f}>"
        )
    return (
        f"<X-point: {sep.Name} Enabled: {sep.Enabled}, " f"Found: {sep.IsSeparatrix}>"
    )


@add_method(Limiters)
def olim_outline(lims: Limiters) -> ArrayN2:
    """Get the outer limiter outline, if it exists
    return as a 2D array of R, Z points
    and arranged in a clockwise fashion.
    """
    olim = np.array(
        [[[lim.R1, lim.Z1], [lim.R2, lim.Z2]] for lim in lims if lim.Enabled > 0]  # type: ignore[attr-defined]
    )
    outline = segments_to_polygon(olim)
    return outline if is_polygon_clockwise(outline) else np.flipud(outline)


@add_method(Limiters)
def ilim_outline(lims: Limiters) -> ArrayN2:
    """Get the inner limiter outline, if it exists
    return as a 2D array of R, Z points
    and arranged in a counter-clockwise fashion.
    """
    ilim = np.array(
        [[[lim.R1, lim.Z1], [lim.R2, lim.Z2]] for lim in lims if lim.Enabled < 0]  # type: ignore[attr-defined]
    )
    outline = segments_to_polygon(ilim)
    return outline if not is_polygon_clockwise(outline) else np.flipud(outline)


# because we need both plasma and psigrid.. (even though it should be just in psigrid)
@add_method(Machine)
def integral_flux_function(
    m: Machine,
    data: npt.ArrayLike,
    average: bool = False,
    psi_norm: npt.ArrayLike | None = None,
) -> np.ndarray:
    """Calculate the flux surface integral of data on the psi grid as a flux function.
    Parameters
    ----------
    psi : PsiGrid
        The psi grid.

    data : npt.ArrayLike
        The function to integrate, defined on the psi grid.

    Returns
    -------
    np.ndarray
        The flux surface average of data on the 1d normalized psi grid.
    """
    integrand_mv = m.PsiGrid.get_new_integrand()
    integrand = np.asarray(integrand_mv)
    if not hasattr(data, "shape"):
        data = np.asarray(data)
    if integrand.shape != data.shape:
        raise ValueError("Data shape does not match psigrid shape")
    one_over_b_mv = m.PsiGrid.get_new_integrand()
    one_over_b = np.asarray(one_over_b_mv)
    np.copyto(one_over_b, 1.0 / np.sqrt(np.asarray(m.Plasma.B2)))
    np.copyto(integrand, data * one_over_b)
    if psi_norm is None:
        psi_norm = np.asarray(m.Plasma.PsiX_pr)
    return np.array(
        [
            m.PsiGrid.contour_integral(integrand_mv, psi_n, False)
            / m.PsiGrid.contour_integral(one_over_b_mv, psi_n, False)
            if average
            else m.PsiGrid.contour_integral(integrand_mv, psi_n, False)
            for psi_n in psi_norm
        ]
    )


def area_of_first_FCFS(m: Machine) -> float:
    """Calculate the area of the first flux surface.
    Parameters
    ----------
    psi : PsiGrid
        The psi grid.

    Returns
    -------
    float
        The area of the first flux surface.
    """
    if np.isclose(m.PsiGrid.PsiFCFS, m.PsiGrid.PsiMagAxis):
        return 0.0

    fcfsr, fcfsx = m.PsiGrid.get_contour(0)
    return area_of_polygon(np.column_stack((fcfsr, fcfsx)))


def rho_toroidal(m: Machine) -> np.ndarray:
    """Calculate the toroidal flux coordinate rho_toroidal from the normalized psi grid.
    Parameters
    ----------
    psi : PsiGrid
        The psi grid.

    Returns
    -------
    np.ndarray
        The toroidal flux coordinate rho_toroidal on the 1d normalized psi grid.
    """

    G = np.asarray(m.Plasma.G)
    R, _ = np.meshgrid(m.PsiGrid.R, m.PsiGrid.Z)
    dphi_dpsi_over_B0R0 = m.integral_flux_function(G / R**2) / (2 * np.pi)
    # I think this isn't quite right, but it is close for the first point
    phi_over_B0R0 = cumulative_trapezoid(
        dphi_dpsi_over_B0R0, m.Plasma.Psi_pr, initial=area_of_first_FCFS(m)
    )
    return np.sqrt(phi_over_B0R0 * m.Plasma.R0 / np.pi)
