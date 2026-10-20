#!/usr/bin/env python3
"""Converts a parameter file of the original (Fortran) YSpec to the JSON input
of this one.

    tools/legacy_to_json.py input.txt > input.json

The original file is read as the Fortran read it, by counting lines: two
header lines, then each value preceded by one blank line and one comment line
(whatever they hold), and after the number of receivers two more lines and
then a latitude and longitude per receiver. As in a Fortran list-directed
read, blank lines before a value are skipped. Paths may be quoted; unlike the
Fortran, an unquoted path is not cut short at its first '/'.
"""
import json
import sys


class Reader:
    """Reads values by counting lines, as the Fortran does."""

    def __init__(self, path):
        self.lines = open(path).read().splitlines()
        self.position = 0

    def tokens(self, skip):
        """Skips lines, then returns the tokens of the next non-blank one."""
        self.position += skip
        while not self.lines[self.position].strip():
            self.position += 1
        line = self.lines[self.position]
        self.position += 1
        return line.split()

    def word(self, skip=2):
        token = self.tokens(skip)[0]
        return token.strip('"\'')


def main(path):
    r = Reader(path)
    prefix = r.word(4)
    model = r.word()
    attenuation = int(r.word()) == 1
    gravity = ['none', 'cowling', 'full'][int(r.word())]
    quantity = ['displacement', 'velocity', 'acceleration'][int(r.word())]
    corrections = int(r.word()) == 1
    lmin, lmax = int(r.word()), int(r.word())
    fmin, fmax = float(r.word()), float(r.word())
    length, step = float(r.word()), float(r.word())
    corners = [float(r.word()) for _ in range(4)]
    depth, lat, lon = float(r.word()), float(r.word()), float(r.word())
    moment = dict(zip(['rr', 'rt', 'rp', 'tt', 'tp', 'pp'],
                      [float(r.word()) for _ in range(6)]))
    receiver_depth = float(r.word())
    count = int(r.word())
    locations = []
    skip = 2
    for _ in range(count):
        tokens = r.tokens(skip)
        skip = 0
        locations.append([float(tokens[0]), float(tokens[1])])

    document = {
        'model': model,
        'attenuation': attenuation,
        'gravity': gravity,
        'output': {
            'prefix': prefix,
            'quantity': quantity,
            'potential_and_tilt_corrections': corrections,
        },
        'degrees': {'min': lmin, 'max': lmax},
        'frequencies_mHz': {'min': fmin, 'max': fmax},
        'time_series': {'length_min': length, 'time_step_s': step},
        'filter_mHz': corners,
        'source': {
            'depth_km': depth,
            'latitude': lat,
            'longitude': lon,
            'moment_tensor_Nm': moment,
        },
        'receivers': {'depth_km': receiver_depth, 'locations': locations},
    }
    text = json.dumps(document, indent=2)
    # Keep each receiver on one line.
    import re
    text = re.sub(r'\[\s+(-?[\d.e+-]+),\s+(-?[\d.e+-]+)\s+\]', r'[\1, \2]', text)
    print(text)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
