#!/usr/bin/env python3
"""Writes gallery.json, the input of the gallery of kernels: the traveltimes
of P, PP, S, ScS and SS at one receiver, each measured in a window about its
arrival in PREM, from TauP, one dominant period before it and one and a half
after, with tapers of half a period.

    poetry -P ../../plotting run python make_input.py

The gallery.json beside it was written so, and the example runs without it;
this is how to change the phases, the distance or the band.
"""

import json
import pathlib

from obspy.taup import TauPyModel

DISTANCE = 60.0  # degrees
AZIMUTH = 70.0   # degrees, clockwise from north at the source
DEPTH = 29.5     # km
PERIOD = 25.0    # s, the dominant period of the band below
PHASES = [("P", "Z"), ("PP", "Z"), ("S", "T"), ("ScS", "T"), ("SS", "T")]
# The arrivals that may share a window, which are printed.
OTHERS = ["P", "pP", "sP", "PcP", "pPcP", "PP", "pPP", "sPP", "PPP", "PS", "SP", "S",
          "sS",
          "ScS", "sScS", "SS", "sSS", "SSS", "SKS", "SKKS", "ScP", "PcS"]


def main():
    model = TauPyModel("prem")
    measurements = []
    for phase, component in PHASES:
        arrivals = model.get_travel_times(source_depth_in_km=DEPTH,
                                          distance_in_degree=DISTANCE,
                                          phase_list=[phase])
        time = min(a.time for a in arrivals)
        start = round(time - PERIOD)
        end = round(time + 1.5 * PERIOD)
        others = sorted({(round(a.time), a.name) for a in model.get_travel_times(
            source_depth_in_km=DEPTH, distance_in_degree=DISTANCE,
            phase_list=OTHERS) if start <= a.time <= end and a.name != phase})
        print(f"{phase:4s} on {component}: arrives at {time:7.1f} s, "
              f"window {start} to {end} s; with "
              + (", ".join(f"{name} {t} s" for t, name in others) or "nothing else"))
        measurements.append({
            "name": phase,
            "component": component,
            "traveltime": {"window_s": [start, end], "taper_s": PERIOD / 2},
        })
    gallery = {
        "model": "../prem.200",
        "gravity": "full",
        "output": {"prefix": "gallery", "quantity": "velocity"},
        # Enough degrees for the fundamental surface waves up to the top of
        # the band, l ~ 2 pi f a / c ~ 410 at 40 mHz: cut off below that,
        # they would be cut off in frequency where their branches pass the
        # last degree, and their sinc tails run through the whole record,
        # windows and all. Nothing below the band is computed.
        "degrees": {"max": 450},
        "frequencies_mHz": {"min": 7, "max": 40},
        # A time step at which the samples are just within a power of two,
        # for the fewest frequencies the damping allows.
        "time_series": {"length_min": 30, "time_step_s": 1.2},
        "filter_mHz": [8, 12, 35, 40],
        "sources": [
            {"depth_km": DEPTH, "latitude": 36.87, "longitude": 35.58,
             "moment_tensor_Nm": {"rr": 4.22e17, "rt": 3.08e17,
                                  "rp": 6.99e17, "tt": -2.976e18,
                                  "tp": 6.54e17, "pp": 2.55e18}}
        ],
        "receiver": {"depth_km": 3, "distance_deg": DISTANCE,
                     "azimuth_deg": AZIMUTH},
        "measurements": measurements,
        "kernels": {
            "parameters": "velocities",
            "scaling": "relative",
            "mesh": {"order": 4, "max_element_km": 100},
            "sections": {"margin_deg": 15, "half_width_deg": 15,
                         "spacing_deg": 0.25, "radial_spacing_km": 30},
        },
    }
    path = pathlib.Path(__file__).with_name("gallery.json")
    path.write_text(json.dumps(gallery, indent=2) + "\n")
    print(f"wrote {path.name}")


if __name__ == "__main__":
    main()
