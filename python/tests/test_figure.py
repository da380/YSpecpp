import numpy as np

import yspec
from yspec import _core


def test_default_and_override_figures_match_independent_prem(model):
    figure = model.figure()
    assert figure.enabled
    assert figure.rotation_rate == _core.EARTH_ROTATION_RATE
    assert np.isclose(figure.surface_ellipticity, 0.003334826529126116,
                      rtol=0, atol=2e-11)
    slower = model.figure(figure.rotation_rate / 2)
    assert np.isclose(slower.surface_ellipticity /
                      figure.surface_ellipticity, 0.25, rtol=1e-12)


def test_disabled_and_zero_spin_figures_are_identity(model):
    earth = model.figure()
    sphere = model.figure(latitude_correction=False)
    zero_spin = model.figure(rotation_rate_rad_s=0)
    assert not sphere.enabled
    assert zero_spin.enabled
    assert sphere.surface_ellipticity == 0
    assert zero_spin.surface_ellipticity == 0
    for figure in (sphere, zero_spin):
        assert figure.geocentric_latitude(45) == 45
        assert figure.geographic_latitude(-90) == -90
    assert earth.geocentric_latitude(45) != 45
    assert np.isclose(earth.geographic_latitude(
        earth.geocentric_latitude(45)), 45, atol=2e-14)


def test_calculation_exposes_coordinate_settings():
    calculation = yspec.Calculation(0, (1, 2), 1)
    assert calculation.latitude_correction
    assert calculation.rotation_rate_rad_s == _core.EARTH_ROTATION_RATE
    assert calculation._core().latitude_correction
    disabled = yspec.Calculation(0, (1, 2), 1, latitude_correction=False,
                                 rotation_rate_rad_s=0)
    core = disabled._core()
    assert not core.latitude_correction
    assert core.rotation_rate == 0


def test_geographic_great_circle_requires_explicit_figure(model):
    source = (10.0, 20.0)
    end_geo = yspec.great_circle(source, 90, [10], figure=model.figure())[0]
    end_geo_coords = (end_geo.latitude, end_geo.longitude)
    end_sphere = yspec.great_circle(source, 90, [10])[0]
    assert end_geo_coords != (end_sphere.latitude, end_sphere.longitude)
    recovered = model.figure().geocentric_latitude(end_geo.latitude)
    assert -90 <= recovered <= 90
