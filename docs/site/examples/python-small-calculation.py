"""A small source-built Python calculation used by the YSpecpp manual."""

import numpy as np
import yspec


yspec.set_threads(2)
model = yspec.Model.prem()
calculation = yspec.Calculation(
    degrees=60,
    frequencies_mhz=(0.5, 6.0),
    length_min=20,
    time_step_s=20,
    filter_mhz=(0.5, 1.0, 4.0, 6.0),
    quantity="displacement",
)
source = yspec.Source(
    latitude=0,
    longitude=0,
    depth_km=29.5,
    moment_tensor={
        "rr": 1.0e18,
        "rt": 0.0,
        "rp": 0.0,
        "tt": 0.0,
        "tp": 0.0,
        "pp": 0.0,
    },
)
receiver = yspec.Receiver(latitude=0, longitude=30, depth_km=3)
result = yspec.seismograms(model, calculation, source, receiver, progress=True)

assert result.data.shape == (1, 1, 3, 61)
assert result.components == ("Z", "N", "E")
assert result.time_s.shape == (61,)
assert np.isfinite(result.data).all()
assert np.max(np.abs(result.trace("Z"))) > 0
print("data shape [source, receiver, component, sample]:", result.data.shape)
print("components:", result.components)
print("time_s shape:", result.time_s.shape)
print("vertical peak displacement (m):", np.max(np.abs(result.trace("Z"))))
