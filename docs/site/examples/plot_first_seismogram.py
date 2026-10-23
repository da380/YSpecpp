#!/usr/bin/env python3
"""Plot the vertical trace from a yspec four-column output file.

Usage: plot_first_seismogram.py TRACE OUTPUT [INPUT_JSON]
The JSON is optional for compatibility; when supplied, it labels the quantity
and derives the plotted receiver distance from the run's receiver table.
"""
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv) not in (3, 4):
    raise SystemExit(__doc__.splitlines()[2])
source, destination = map(Path, sys.argv[1:3])
input_path = Path(sys.argv[3]) if len(sys.argv) == 4 else None
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

quantity = None
unit = None
if input_path is not None:
    settings = json.loads(input_path.read_text())
    quantity = settings.get('output', {}).get('quantity', 'displacement')
    unit = {'displacement': 'm', 'velocity': 'm/s',
            'acceleration': 'm/s²'}.get(quantity)
    if unit is None:
        raise SystemExit(f'unsupported output quantity: {quantity}')
    receiver_groups = settings.get('receivers', [])
    if isinstance(receiver_groups, dict):
        receiver_groups = [receiver_groups]
    if any(group.get('instrument') for group in receiver_groups):
        instrument_path = input_path.parent / next(
            group['instrument'] for group in receiver_groups
            if group.get('instrument'))
        for line in instrument_path.read_text().splitlines():
            if line.upper().startswith('* OUTPUT UNIT :'):
                unit = line.split(':', 1)[1].strip()
                break
        else:
            unit = 'instrument output (unit unspecified)'

distance = None
receiver_table = source.with_name(f'{source.name.rsplit(".", 1)[0]}.receivers')
if receiver_table.exists():
    lines = [line.split() for line in receiver_table.read_text().splitlines()
             if line.strip()]
    if lines and lines[0][0] == '#':
        header = lines[0][1:]
        if 'file' in header and 'distance_deg' in header:
            file_column, distance_column = (header.index('file'),
                                            header.index('distance_deg'))
            match = next((row for row in lines[1:]
                          if len(row) > max(file_column, distance_column) and
                          Path(row[file_column]).name == source.name), None)
            if match is not None:
                distance = float(match[distance_column])

quantity_label = quantity or 'vertical component'
title = quantity_label.capitalize()
if distance is not None:
    title += f' · {distance:g}° epicentral distance'
ylabel = 'Up component' + (f' ({unit})' if unit else '')
matplotlib.rcParams['svg.hashsalt'] = 'yspecpp-first-seismogram'
fig, ax = plt.subplots(figsize=(8.2, 3.6), constrained_layout=True)
ax.plot(times, vertical, color='#176b70', linewidth=1.1)
ax.set(title=title, xlabel='Time after origin (s)', ylabel=ylabel)
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
peak_unit = f' {unit}' if unit else ' (output unit unspecified)'
print(f'{len(rows)} samples; vertical peak absolute amplitude '
      f'{max(map(abs, vertical)):.9g}{peak_unit}; wrote {destination}')
