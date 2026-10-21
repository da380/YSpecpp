#!/usr/bin/env python3
"""Run the small downloadable manual fixtures and compile its library example."""
import argparse
import math
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def run(command, env=None):
    print('+', ' '.join(map(str, command)), flush=True)
    subprocess.run(command, check=True, cwd=ROOT, env=env)


def finite_rows(path, columns, skip=0):
    lines = path.read_text().splitlines()[skip:]
    if not lines:
        raise RuntimeError(f'{path}: no data rows')
    for number, line in enumerate(lines, skip + 1):
        values = line.split()
        if len(values) != columns:
            raise RuntimeError(f'{path}:{number}: expected {columns} columns, got {len(values)}')
        if not all(map(lambda value: re.fullmatch(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?', value), values)):
            raise RuntimeError(f'{path}:{number}: expected finite numeric values')
        if not all(math.isfinite(float(value)) for value in values):
            raise RuntimeError(f'{path}:{number}: non-finite value')
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True, help='configured Release/OpenMP CMake build directory')
    args = parser.parse_args()
    build = args.build.resolve()
    for binary in ('bin/yspec', 'bin/yspec-snapshot', 'libYSpec.a'):
        if not (build / binary).is_file():
            parser.error(f'missing {build / binary}; build yspec, yspec-snapshot and YSpec first')
    cache = (build / 'CMakeCache.txt').read_text()
    compiler = re.search(r'^CMAKE_CXX_COMPILER:(?:FILEPATH|STRING)=(.+)$', cache, re.M)
    openmp = re.search(r'^YSPEC_WITH_OPENMP:BOOL=(ON|OFF)$', cache, re.M)
    if not compiler or not openmp or openmp.group(1) != 'ON':
        parser.error('this focused example check needs the Release/OpenMP build used by docs CI')
    env = os.environ | {'OMP_NUM_THREADS': '2'}
    with tempfile.TemporaryDirectory(prefix='yspec-doc-examples-', dir=build) as temporary:
        work = Path(temporary)
        fixtures = {
            'snapshot': 'snapshot-small',
            'workflows': 'workflows',
            'first-seismogram': 'first-seismogram',
        }
        for name, fixture in fixtures.items():
            (work / name).mkdir()
            shutil.copy2(ROOT / f'docs/site/examples/{fixture}.json', work / name / 'input.json')
            shutil.copy2(ROOT / 'examples/prem.200', work / name / 'prem.200')
        shutil.copy2(ROOT / 'docs/site/examples/station.pz', work / 'workflows/station.pz')
        run([build / 'bin/yspec-snapshot', work / 'snapshot/input.json'], env)
        snapshot = work / 'snapshot/snapshot.1'
        rows = finite_rows(snapshot, 8, skip=2)
        if len(rows) != 4:
            raise RuntimeError(f'{snapshot}: expected four grid rows, got {len(rows)}')
        first_run = work / 'first-seismogram'
        run([build / 'bin/yspec', first_run / 'input.json'], env)
        first_trace = first_run / 'first-seismogram.out.1'
        # The traces are the files with a numbered suffix; beside them yspec
        # writes first-seismogram.out.receivers, a table of each trace's
        # receiver, checked below.
        trace_files = sorted(path for path in first_run.glob('first-seismogram.out.*')
                             if path.suffix[1:].isdigit())
        if trace_files != [first_trace] or not first_trace.is_file():
            raise RuntimeError('first-seismogram run did not create exactly first-seismogram.out.1')
        receivers = first_run / 'first-seismogram.out.receivers'
        if not receivers.is_file():
            raise RuntimeError(f'{receivers}: not written')
        rows = receivers.read_text().splitlines()
        if len(rows) != 2 or not rows[0].startswith('#') or rows[1].split()[0] != first_trace.name:
            raise RuntimeError(f'{receivers}: expected a header and one row naming {first_trace.name}')
        rows = finite_rows(first_trace, 4)
        if len(rows) != 91:
            raise RuntimeError(f'{first_trace}: expected 91 time samples, got {len(rows)}')
        times = [float(row.split()[0]) for row in rows]
        if any(not math.isclose(time, 20.0 * index, rel_tol=0.0, abs_tol=1e-9)
               for index, time in enumerate(times)):
            raise RuntimeError(f'{first_trace}: expected times from 0 to 1800 s at 20 s intervals')
        if not any(float(row.split()[1]) != 0.0 for row in rows):
            raise RuntimeError(f'{first_trace}: expected a nonzero vertical displacement trace')
        run([build / 'bin/yspec', work / 'workflows/input.json'], env)
        traces = sorted((work / 'workflows').glob('multi.event-*.*'))
        expected = [work / 'workflows' / f'multi.{source}.{receiver}'
                    for source in ('event-a', 'event-b') for receiver in range(1, 5)]
        if sorted(traces) != sorted(expected):
            raise RuntimeError('multi-source run did not create the eight documented filenames')
        for trace in expected:
            rows = finite_rows(trace, 4)
            if len(rows) != 31:
                raise RuntimeError(f'{trace}: expected 31 time samples, got {len(rows)}')

        executable = work / 'library-workflows'
        command = [compiler.group(1), '-std=c++23', '-O3', '-I', str(ROOT),
                   str(ROOT / 'docs/site/examples/library-workflows.cpp'),
                   str(build / 'libYSpec.a'), '-lfftw3', '-lfftw3f', '-lfftw3l',
                   '-fopenmp', '-o', str(executable)]
        run(command)
        run([executable, ROOT / 'examples/prem.200'], env)
    print('Checked first-seismogram and snapshot fixtures, multi-source/profile/instrument traces and strain/kernel library consumer.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
