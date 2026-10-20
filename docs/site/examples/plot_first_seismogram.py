#!/usr/bin/env python3
"""Plot the vertical trace from a yspec four-column output file."""
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

source, destination = map(Path, sys.argv[1:3])
rows = []
for line in source.read_text().splitlines():
    fields = line.split()
    if fields:
        rows.append(tuple(map(float, fields)))
if not rows or any(len(row) != 4 for row in rows):
    raise SystemExit('expected four numeric columns per output row')
if any(not math.isfinite(value) for row in rows for value in row):
    raise SystemExit('output contains a non-finite value')
times, vertical = [r[0] for r in rows], [r[1] for r in rows]
if max(map(abs, vertical)) == 0:
    raise SystemExit('vertical trace is identically zero')
matplotlib.rcParams['svg.hashsalt'] = 'yspecpp-first-seismogram'
fig, ax = plt.subplots(figsize=(8.2, 3.6), constrained_layout=True)
ax.plot(times, vertical, color='#176b70', linewidth=1.1)
ax.set(title='Vertical displacement · 30° epicentral distance',
       xlabel='Time after origin (s)', ylabel='Up displacement (m)')
ax.grid(True, color='#e3e9e7', linewidth=.7)
ax.margins(x=0)
destination.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(destination, format=destination.suffix.lstrip('.') or 'svg',
            metadata={'Creator': 'YSpecpp tutorial plot script', 'Date': None})
svg = destination.read_text()
destination.write_text(''.join(
    line.rstrip('\r\n').rstrip(' \t') + line[len(line.rstrip('\r\n')):]
    for line in svg.splitlines(keepends=True)
))
print(f'{len(rows)} samples; vertical peak absolute amplitude {max(map(abs, vertical)):.9g} m; wrote {destination}')
