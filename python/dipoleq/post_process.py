"""
Extra analysis routines for after a solve.

These are pure python functions that extend the machine object functionality.
"""

from collections.abc import Callable
from typing import Any, TypeVar

from .core import Machine, Separatrix

# lets pirate onto the c-class with extra functions
# just have to import this file... nice.
_T = TypeVar("_T")


def _add_method(cls: type[_T]) -> Callable[..., Any]:
    def decorator(func: Callable[[_T, Any], Any]) -> Callable[[_T, Any], Any]:
        setattr(cls, func.__name__, func)
        return func

    return decorator


@_add_method(Machine)
def is_outer_limited(m: Machine) -> bool:
    """Check if the equilibrium is diverted"""
    return any(lim.Enabled > 0 and lim.PsiLim <= m.PsiGrid.PsiLim for lim in m.Limiters)


@_add_method(Machine)
def is_diverted(m: Machine) -> bool:
    """Check if the equilibrium is diverted."""
    return not is_outer_limited(m)


@_add_method(Machine)
def get_outer_limiter_contact_point(m: Machine) -> tuple[float, float] | None:
    """Get the outer limiter contact point"""
    for lim in m.Limiters:
        if lim.Enabled > 0 and lim.PsiMin <= m.PsiGrid.PsiLim:
            return lim.Rmin, lim.Zmin
    return None


@_add_method(Machine)
def get_inner_limiter_contact_point(m: Machine) -> tuple[float, float] | None:
    """Get the inner limiter contact point"""
    for lim in m.Limiters:
        if lim.Enabled < 0 and lim.PsiMin >= m.PsiGrid.PsiLim:
            return lim.Rmin, lim.Zmin
    return None


@_add_method(Machine)
def get_xpoints(m: Machine) -> list[Separatrix]:
    """Get the X-points"""
    valid_seps = [
        sep
        for sep in m.Seps
        if sep.Enabled and sep.IsSeparatrix and sep.Psi > m.PsiGrid.PsiAxis
    ]
    valid_seps.sort(key=lambda sep: sep.Psi)
    return valid_seps
