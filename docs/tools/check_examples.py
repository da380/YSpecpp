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


def numeric_rows(path, columns=None, tagged=False):
    rows = []
    for number, line in enumerate(path.read_text().splitlines(), 1):
        if not line or line.startswith('#'):
            continue
        fields = line.split()
        if tagged:
            if fields[0] not in ('path', 'across'):
                raise RuntimeError(f'{path}:{number}: expected path/across section tag')
            fields = fields[1:]
        if columns is not None and len(fields) != columns:
            raise RuntimeError(f'{path}:{number}: expected {columns} numeric columns, got {len(fields)}')
        try:
            values = [float(value) for value in fields]
        except ValueError as error:
            raise RuntimeError(f'{path}:{number}: expected numeric data') from error
        if not all(map(math.isfinite, values)):
            raise RuntimeError(f'{path}:{number}: non-finite value')
        rows.append(values)
    if not rows:
        raise RuntimeError(f'{path}: no data rows')
    return rows


def copy_kernel_fixture(directory, fixture):
    directory.mkdir()
    shutil.copy2(ROOT / f'docs/site/examples/{fixture}', directory / 'input.json')
    shutil.copy2(ROOT / 'examples/prem.200', directory / 'prem.200')


def check_kernel_outputs(directory, prefix):
    measurement = directory / f'{prefix}.measurement'
    data = numeric_rows(measurement, 4)
    if len(data) != 61 or not any(row[1] != 0 for row in data):
        raise RuntimeError(f'{measurement}: expected 61 samples and a nonzero seismogram')
    if not any(row[2] != 0 for row in data) or not any(row[3] != 0 for row in data):
        raise RuntimeError(f'{measurement}: expected a nonzero window and adjoint source')

    radial = directory / f'{prefix}.kernels'
    headers = [line for line in radial.read_text().splitlines() if line.startswith('#')]
    at_source = next((line.split() for line in headers
                      if line.startswith('# at_source ')), None)
    if not at_source:
        raise RuntimeError(f'{radial}: missing source-point P_X terms')
    columns = next((line.split() for line in reversed(headers)
                    if line.startswith('# radius_km depth_km weight_m3')), None)
    if not columns or len(columns) < 4:
        raise RuntimeError(f'{radial}: missing radial kernel parameter columns')
    parameter_count = len(columns) - 4
    try:
        point_terms = [float(value) for value in at_source[2:]]
    except ValueError as error:
        raise RuntimeError(f'{radial}: invalid source-point values') from error
    if len(point_terms) != parameter_count or not all(map(math.isfinite, point_terms)):
        raise RuntimeError(f'{radial}: source-point P_X count/values do not match the parameter columns')
    values = numeric_rows(radial, len(columns) - 1)
    if not any(value != 0 for row in values for value in row[3:]):
        raise RuntimeError(f'{radial}: expected nonzero radial kernel values')

    section_file = directory / f'{prefix}.sections'
    section_lines = section_file.read_text().splitlines()
    section_header = next((line.split() for line in section_lines
                           if line.startswith('# section layer depth_km')), None)
    if not section_header:
        raise RuntimeError(f'{section_file}: missing path/across column header')
    sampled_parameters = len(section_header) - 8  # '# section' + six coordinates
    if sampled_parameters < 1:
        raise RuntimeError(f'{section_file}: missing sampled-kernel columns')
    sections = numeric_rows(section_file, columns=6 + sampled_parameters,
                            tagged=True)
    section_tags = {line.split()[0] for line in section_lines
                    if line and not line.startswith('#')}
    if section_tags != {'path', 'across'}:
        raise RuntimeError(f'{prefix}.sections: expected both path and across grids')
    if not any(value != 0 for row in sections for value in row[6:]):
        raise RuntimeError(f'{prefix}.sections: expected nonzero sampled kernels')
    map_file = directory / f'{prefix}.maps'
    map_header = next((line.split() for line in map_file.read_text().splitlines()
                       if line.startswith('# layer depth_km')), None)
    if not map_header or len(map_header) - 7 != sampled_parameters:
        raise RuntimeError(f'{map_file}: sampled columns differ from sections')
    maps = numeric_rows(map_file, columns=6 + sampled_parameters)
    if not any(row[6] != 0 for row in maps):
        raise RuntimeError(f'{prefix}.maps: expected nonzero sampled kernels')


def check_png(path):
    if not path.is_file() or path.stat().st_size < 1000 or path.read_bytes()[:8] != b'\x89PNG\r\n\x1a\n':
        raise RuntimeError(f'{path}: missing or invalid PNG plot')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True, help='configured Release/OpenMP CMake build directory')
    parser.add_argument('--python', type=Path, required=True, help='Python executable with this checkout installed using the plot extra')
    args = parser.parse_args()
    build = args.build.resolve()
    python = args.python.absolute()  # Preserve the venv executable symlink.
    for binary in ('bin/yspec', 'bin/yspec-snapshot', 'bin/yspec-kernels', 'libYSpec.a'):
        if not (build / binary).is_file():
            parser.error(f'missing {build / binary}; build yspec, yspec-snapshot and YSpec first')
    cache = (build / 'CMakeCache.txt').read_text()
    compiler = re.search(r'^CMAKE_CXX_COMPILER:(?:FILEPATH|STRING|UNINITIALIZED)=(.+)$', cache, re.M)
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

        # Check the documented Python API and package entry points from the
        # exact checkout installed by the caller.
        python_env = env | {'MPLBACKEND': 'Agg'}
        python_script = ROOT / 'docs/site/examples/python-small-calculation.py'
        run([python, python_script], python_env)
        python_seismo = work / 'python-seismogram'
        python_seismo.mkdir()
        shutil.copy2(ROOT / 'docs/site/examples/first-seismogram.json', python_seismo / 'input.json')
        shutil.copy2(ROOT / 'examples/prem.200', python_seismo / 'prem.200')
        run([python.parent / 'yspec', python_seismo / 'input.json'], python_env)
        python_trace = python_seismo / 'first-seismogram.out.1'
        rows = finite_rows(python_trace, 4)
        if len(rows) != 91 or not any(float(row.split()[1]) != 0 for row in rows):
            raise RuntimeError(f'{python_trace}: expected 91 samples and a nonzero vertical trace')

        python_snapshot = work / 'python-snapshot'
        python_snapshot.mkdir()
        # The snapshot package wrapper is exercised on the separate snapshot fixture.
        shutil.copy2(ROOT / 'docs/site/examples/snapshot-small.json', python_snapshot / 'input.json')
        shutil.copy2(ROOT / 'examples/prem.200', python_snapshot / 'prem.200')
        run([python.parent / 'yspec-snapshot', python_snapshot / 'input.json'], python_env)
        if len(finite_rows(python_snapshot / 'snapshot.1', 8, skip=2)) != 4:
            raise RuntimeError('Python snapshot wrapper did not produce the four-row fixture')

        cpp_kernels = work / 'cpp-kernels'
        copy_kernel_fixture(cpp_kernels, 'kernels-small.json')
        run([build / 'bin/yspec-kernels', cpp_kernels / 'input.json'], env)
        check_kernel_outputs(cpp_kernels, 'kernel-small')

        python_kernels = work / 'python-kernels'
        copy_kernel_fixture(python_kernels, 'kernels-small.json')
        run([python.parent / 'yspec-kernels', python_kernels / 'input.json'], python_env)
        check_kernel_outputs(python_kernels, 'kernel-small')

        plotter = python.parent / 'yspec-plot'
        run([plotter, 'traces', python_seismo / 'input.json', '--receivers', '1',
             '--components', 'ZRT', '-o', python_seismo / 'traces.png'], python_env)
        run([plotter, 'kernels', python_kernels / 'input.json', '--measurement', 'S',
             '--deepest', '1500', '-o', python_kernels / 'radial.png'], python_env)
        run([plotter, 'banana', python_kernels / 'input.json', '--measurement', 'S',
             '--parameter', 'vs', '-o', python_kernels / 'section.png'], python_env)
        for image in (python_seismo / 'traces.png', python_kernels / 'radial.png',
                      python_kernels / 'section.png'):
            check_png(image)

        # The Python source build deliberately disables NetCDF. Confirm that
        # volume writing is refused there, while the C++ build writes and
        # exposes the expected model variables when its option is enabled.
        volume_off = work / 'python-volume-off'
        copy_kernel_fixture(volume_off, 'kernels-volume-small.json')
        refused = subprocess.run([python.parent / 'yspec-kernels',
                                  volume_off / 'input.json'], cwd=ROOT,
                                 env=python_env, capture_output=True, text=True)
        if refused.returncode == 0 or 'built without netCDF' not in refused.stderr:
            raise RuntimeError('Python package volume route did not report its NetCDF-disabled boundary')

        netcdf = re.search(r'^YSPEC_WITH_NETCDF:BOOL=(ON|OFF)$', cache, re.M)
        if not netcdf:
            raise RuntimeError('cannot determine C++ build NetCDF setting')
        if netcdf.group(1) == 'ON':
            volume_on = work / 'cpp-volume-on'
            copy_kernel_fixture(volume_on, 'kernels-volume-small.json')
            run([build / 'bin/yspec-kernels', volume_on / 'input.json'], env)
            ncfile = volume_on / 'kernel-volume-small.nc'
            inspect_netcdf = ("import netCDF4,numpy as np,sys; f=netCDF4.Dataset(sys.argv[1]); "
                "g=next(iter(f.groups.values())); "
                "assert {k:len(f.dimensions[k]) for k in ('layer','interface','bound')} == {'layer':13,'interface':13,'bound':2}; "
                "assert {k:len(g.dimensions[k]) for k in ('node','element','latitude','longitude','bound','pole')} == {'node':26,'element':13,'latitude':13,'longitude':25,'bound':2,'pole':2}; "
                "names=('K_vs','K_vp','K_qmu','K_qkappa'); "
                "assert all(n in g.variables and n+'_poles' in g.variables for n in names); "
                "assert all(g.variables[n].dimensions == ('node','latitude','longitude') for n in names); "
                "assert all(g.variables[n+'_poles'].dimensions == ('node','pole') for n in names); "
                "assert all(g.variables[n].units == 'm-3 s' for n in names); "
                "assert all(np.isfinite(g.variables[n][:]).all() and np.any(g.variables[n][:] != 0) for n in names); f.close()")
            run([python, '-c', inspect_netcdf, str(ncfile)], python_env)

        executable = work / 'library-workflows'
        command = [compiler.group(1), '-std=c++23', '-O3', '-I', str(ROOT),
                   str(ROOT / 'docs/site/examples/library-workflows.cpp'),
                   str(build / 'libYSpec.a'), '-lfftw3', '-lfftw3f', '-lfftw3l',
                   '-fopenmp', '-o', str(executable)]
        run(command)
        run([executable, ROOT / 'examples/prem.200'], env)
    print('Checked C++ and Python first-seismogram/snapshot/kernel routes, bounded Python results, radial/section/map/measurement outputs, three headless plots, NetCDF build boundaries, and the strain/kernel library consumer.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
