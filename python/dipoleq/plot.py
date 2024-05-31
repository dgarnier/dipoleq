"""
Plotting functions for dipoleq
"""

import numpy as np

from . import Machine


def plot_eq(m: Machine, ax=None, **kwargs):
    """Plot the equilibrium"""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib is required for plotting.")
        return

    if ax is None:
        fig, ax = plt.subplots()

    pg = m.PsiGrid
    ax.contour(pg.R, pg.Z, pg.Psi, 100)
    ax.set_xlabel(f"R [m]")
    ax.set_ylabel(f"Z [m]")
    ax.set_aspect("equal")
    LCFS = pg.get_contour(1.0)
    FCFS = pg.get_contour(0.0) if pg.PsiAxis != pg.PsiMagAxis else None
    olim = np.array([[[l.R1, l.Z1], [l.R2, l.Z2]] for l in m.Limiters if l.Enabled > 0])
    ilim = np.array([[[l.R1, l.Z1], [l.R2, l.Z2]] for l in m.Limiters if l.Enabled < 0])

    ax.plot(LCFS[:, 0], LCFS[:, 1], "b--")
    if FCFS is not None:
        ax.plot(FCFS[:, 0], FCFS[:, 1], "b--")

    if olim is not None:
        ax.plot(olim[:, 0], olim[:, 1], "k-")

    if ilim is not None:
        ax.plot(ilim[:, 0], ilim[:, 1], "k-")

    return ax
