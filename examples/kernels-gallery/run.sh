#!/bin/sh
# A gallery of three-dimensional traveltime kernels, each phase measured by
# cross-correlation in a window about its arrival in PREM, from TauP
# (make_input.py wrote the inputs): at 60 degrees, between 25 and 125 s, P and
# PP on the vertical component and S, ScS and SS on the transverse; at 110
# degrees, between 40 and 200 s, the diffracted Pdiff and Sdiff. Each run
# computes its kernels from one forward field and an adjoint field for each
# component; each measurement has its own kernels, drawn on sections along
# the path and across it at its middle, below its seismogram and window.
# About three minutes on four threads, two of them for the first run.
set -e
cd "$(dirname "$0")"
. ../common.sh
run yspec-kernels gallery.json
plot gallery gallery.json -o gallery.png
plot banana gallery.json --measurement S -o banana-S.png
run yspec-kernels diffracted.json
plot gallery diffracted.json -o diffracted.png
