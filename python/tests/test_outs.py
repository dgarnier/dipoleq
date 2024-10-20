# Tests weird outputs...
# GS2 for now, but there are others
# for now just to include coverage
# should add more tests

import os
from pathlib import Path

import numpy as np
from dipoleq import Machine
from dipoleq.util import area_of_polygon

data_dir = Path(os.path.realpath(__file__)).parent / "data"


def test_gs2_out() -> None:  # pylint: disable=missing-function-docstring
    m = Machine.from_file(data_dir / "beta1.yaml")
    m.solve()
    m.write_GS2_geo()


def test_flux_surface_average() -> None:
    m = Machine.from_file(data_dir / "beta1.yaml")
    m.solve()

    B2_avg = m.integral_flux_function(m.Plasma.B2, average=True)
    assert np.asarray(B2_avg.shape) == np.asarray(m.Plasma.B2_pr).shape
    assert np.allclose(B2_avg, m.Plasma.B2_pr)


def test_polygon_area() -> None:
    """Test the area of a polygon, Copilot got it wrong."""
    poly = np.array([[0, 0], [2, 0], [2, 1], [0, 1]])
    assert np.isclose(area_of_polygon(poly), 2)

    poly2 = np.array([[0, 0], [2, 0], [2, 1], [0, 1], [0, 0]])
    assert np.isclose(area_of_polygon(poly2), 2)

    poly3 = np.array([[-3, -2], [-1, 4], [6, 1], [3, 10], [-4, 9]])
    assert np.isclose(area_of_polygon(poly3), 60)
