"""
Plotting functions for dipoleq
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from matplotlib.axes import Axes

    from . import Machine


def plot_eq(m: Machine, ax: Axes | None = None, show_peak: bool = True) -> Axes | None:
    """Plot the equilibrium"""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib is required for plotting.")
        return None

    if ax is None:
        ax = plt.gca()

    pg = m.PsiGrid
    ax.contour(pg.R, pg.Z, pg.Psi, 100)
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.set_aspect("equal")
    LCFS = np.array(pg.get_contour(1.0)).T
    FCFS = np.array(pg.get_contour(0.0)).T if pg.PsiAxis != pg.PsiMagAxis else None
    olim = np.array(
        [[[lim.R1, lim.Z1], [lim.R2, lim.Z2]] for lim in m.Limiters if lim.Enabled > 0]
    )
    ilim = np.array(
        [[[lim.R1, lim.Z1], [lim.R2, lim.Z2]] for lim in m.Limiters if lim.Enabled < 0]
    )

    ax.plot(LCFS[:, 0], LCFS[:, 1], "b--")
    if FCFS is not None:
        ax.plot(FCFS[:, 0], FCFS[:, 1], "b--")

    if olim is not None:
        ax.plot(olim[:, :, 0].flatten(), olim[:, :, 1].flatten(), "k-.")

    if ilim is not None:
        ax.plot(ilim[:, :, 0].flatten(), ilim[:, :, 1].flatten(), "k-.")

    if show_peak:
        pl = m.Plasma
        ipsi_peak = np.argmax(pl.P_pr)
        psix = np.array(pl.PsiX_pr)
        psix_peak = psix[ipsi_peak]
        Rpeak, Zpeak = pg.get_contour(psix_peak)
        ax.plot(Rpeak, Zpeak, "m--")

    return ax
