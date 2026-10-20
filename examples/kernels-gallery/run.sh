#!/bin/sh
# A gallery of three-dimensional traveltime kernels at 60 degrees, between 25
# and 125 s: P and PP on the vertical component, S, ScS and SS on the
# transverse, each measured by cross-correlation in a window about its arrival
# in PREM, from TauP (make_input.py wrote gallery.json). One run computes
# them all, from one forward field and an adjoint field for each component;
# each measurement has its own kernels, drawn on sections along the path and
# across it at its middle, below its seismogram and window. About two minutes
# on four threads.
set -e
cd "$(dirname "$0")"
. ../common.sh
run yspec-kernels gallery.json
plot gallery gallery.json -o gallery.png
plot banana gallery.json --measurement S -o banana-S.png
