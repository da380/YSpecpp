#!/bin/sh
# The three-dimensional kernel of an S-wave traveltime at 60 degrees, measured
# by cross-correlation on the transverse component, relative, in the
# isotropic speeds and quality factors at fixed density: on a section along
# the great circle through the source and the receiver, where it is the
# banana of finite-frequency theory, and on maps at four depths, where its
# Fresnel zones lie along the path. About a minute on four threads; the
# kernels are written as a planetmodel model file, s-wave.nc.
set -e
cd "$(dirname "$0")"
. ../common.sh
run yspec-kernels s-wave.json
plot kernels s-wave.json -o s-wave.png
plot volume s-wave.json --parameter vs -o s-wave-volume.png
