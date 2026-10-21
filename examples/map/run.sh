#!/bin/sh
# The vertical velocity at the sea floor around the source, in
# azimuthal equidistant projection with the source at the centre and its
# antipode on the rim: the surface waves spread out, with the radiation
# pattern of the source, and converge again on the antipode.
set -e
cd "$(dirname "$0")"
. ../common.sh
run yspec-snapshot snapshot.json
plot snapshot snapshot.json --component up -o map.png
