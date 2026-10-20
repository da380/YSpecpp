#!/bin/sh
# The radial kernels of two traveltimes at 60 degrees, measured by
# cross-correlation, relative and in density, the isotropic speeds and the
# quality factors: the S wave on the transverse component, sensitive to the
# shear speed down through the mantle, and the Rayleigh wave on the vertical,
# sensitive to the top few hundred kilometres.
set -e
cd "$(dirname "$0")"
. ../common.sh
run yspec-kernels s-wave.json
run yspec-kernels rayleigh.json
plot kernels s-wave.json -o s-wave.png
plot kernels rayleigh.json --deepest 700 -o rayleigh.png
