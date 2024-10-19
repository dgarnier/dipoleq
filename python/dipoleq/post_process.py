"""
Extra analysis routines for after a solve.

These are pure python functions that extend the machine object functionality.
"""

import numpy as np

from .core import Limiters, Machine, Separatrix
from .util import ArrayN2, add_method, is_polygon_clockwise, segments_to_polygon

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
def flux_surface_average(m: Machine, data: np.ndarray) -> np.ndarray:
    """Calculate the flux surface average of f.
    Parameters
    ----------
    psi : PsiGrid
        The psi grid.

    f : np.ndarray
        The function to average, defined on the psi grid.

    Returns
    -------
    np.ndarray
        The flux surface average of f on the 1d normalized psi grid.
    """
    integrand_mv = m.PsiGrid.get_new_integrand()
    integrand = np.asarray(integrand_mv)
    if integrand.shape != data.shape:
        raise ValueError("Value shape does not match psigrid shape")
    one_over_b_mv = m.PsiGrid.get_new_integrand()
    one_over_b = np.asarray(one_over_b_mv)
    one_over_b[:] = (
        2
        * np.pi
        * np.array(m.PsiGrid.R)[:, None]
        / np.sqrt(np.array(m.Plasma.GradPsi2))
    )
    return np.array(
        [
            m.PsiGrid.contour_integral(integrand_mv, psi_n, False)
            / m.PsiGrid.contour_integral(one_over_b_mv, psi_n, False)
            for psi_n in m.Plasma.PsiX_pr
        ]
    )
