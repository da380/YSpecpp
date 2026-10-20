# YSpec

YSpec computes synthetic seismograms for a spherically symmetric earth model
by direct radial integration of the equations of motion in the frequency
domain (Friederich & Dalkolmo 1995), including self-gravitation (Al-Attar &
Woodhouse 2008). This is a C++23 re-implementation of the original Fortran
YSpec: the same physics, inputs and outputs, with better numerics and a
structure that is easier to extend.

The earth model may be transversely isotropic and anelastic, with any
sequence of solid and fluid regions (a fluid or solid inner core, an ocean or
none). Attenuation is the constant-Q logarithmic dispersion of PREM,
implemented exactly at complex frequency rather than to first order in 1/Q.
Gravity can be left out, included in the Cowling approximation, or included
in full.

## Building

YSpec needs a C++23 compiler (GCC 13 or Clang 18 and later), CMake 3.24 or
later, FFTW3 (in double, single and long double precision, as FFTWpp requires),
and optionally OpenMP. Everything else is fetched by CMake:
[NumericConcepts](https://github.com/da380/NumericConcepts),
[Interpolation](https://github.com/da380/Interpolation),
[FFTWpp](https://github.com/da380/FFTWpp),
[GaussQuad](https://github.com/da380/GaussQuad),
[GSHTrans](https://github.com/da380/GSHTrans),
[planetmodel](https://github.com/da380/planetmodel) and, for the tests,
GoogleTest.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build
```

To build against local checkouts of the dependencies instead of fetching
them, pass `-DFETCHCONTENT_SOURCE_DIR_<NAME>=/path/to/checkout` for each
(`GSHTRANS`, `INTERPOLATION`, ...). `-DYSPEC_WITH_OPENMP=OFF` builds a serial
version, which gives the same results.

## Running

There is a program for each kind of output: `yspec` computes seismograms at
receivers, and `yspec-snapshot` the wavefield at chosen times on a grid
around each source. They read JSON files that share most of their members.

```sh
cd build
bin/yspec examples/input.json
bin/yspec-snapshot examples/snapshot.json
```

The build copies the example inputs and models into `build/examples`, so the
outputs are written there, beside their inputs, and not in the sources.

The number of threads is set with `OMP_NUM_THREADS`. For `yspec`, one file is
written per receiver, `<prefix>.1`, `<prefix>.2`, ..., in the order the
receivers are listed; with several sources, one per source and receiver,
`<prefix>.<source>.1`, .... Each line holds the time after the origin time (s)
and the vertical (up), north and east components, in SI units of the quantity
asked for, or in the units of the instrument's output where there is one. A
table `<prefix>.receivers` lists the files, one line for each source and
receiver: the receiver's latitude, longitude and depth (km), and its distance,
its azimuth at the source and the source's azimuth at it (degrees).

### Examples and plotting

Each directory of `examples/` holds an example, with its input and a script
that runs it and plots the result:

| example | program | what it shows |
|---------|---------|---------------|
| `record-section` | `yspec` | velocity every 5 degrees along a great circle, vertical and transverse |
| `instrument` | `yspec` | one place recorded twice: the ground displacement, and a seismometer's output in counts |
| `depths` | `yspec` | one place at depths from the sea floor to the inner core |
| `gravity` | `yspec` | a long-period seismogram without gravity, in the Cowling approximation and in full |
| `force` | `yspec` | the waves from a point force, of the kind a landslide exerts |
| `cross-section` | `yspec-snapshot` | the wavefield on a great circle through the whole earth at six times |
| `map` | `yspec-snapshot` | the wavefield at the sea floor around the source, over the whole globe |

```sh
cd build
examples/record-section/run.sh     # or examples/run-all.sh for all of them
```

Each takes seconds on four threads (the scripts use four unless
`OMP_NUM_THREADS` says otherwise), and writes its output and figures in its
own directory. The programs are found in the build, on the `PATH`, or in
`YSPEC_BIN`.

The figures are drawn by `yspec-plot`, a small Python package in `plotting/`
built on [planetmodel](https://github.com/da380/planetmodel), whose sections
and colour scales it shares, with maps through cartopy. Its environment is
made once with poetry:

```sh
poetry -P plotting install
poetry -P plotting run yspec-plot section build/examples/record-section/input.json
```

`yspec-plot` takes the input file of a run and finds the output from it:
`section` draws a record section, `traces` the components at a few receivers,
`compare` one receiver in several runs, and `snapshot` a section through the
earth or a map around the source. The components of seismograms are Z, N and
E, or R and T, radial and transverse to the source. Without `-o` the figure is
shown.

### Input

The calculation is read from a JSON file:

```json
{
  "model": "prem.200",
  "attenuation": true,
  "gravity": "full",
  "output": {
    "prefix": "yspec.out",
    "quantity": "velocity",
    "potential_and_tilt_corrections": false
  },
  "degrees": { "min": 0, "max": 600 },
  "frequencies_mHz": { "min": 0.2, "max": 40 },
  "time_series": { "length_min": 60, "time_step_s": 1.0 },
  "filter_mHz": [1.0, 5.0, 38, 40],
  "sources": [
    { "name": "event", "depth_km": 29.5, "latitude": 36.87, "longitude": 35.58,
      "moment_tensor_Nm": { "rr": 4.22e17, "rt": 3.08e17, "rp": 6.99e17,
                            "tt": -2.976e18, "tp": 6.54e17, "pp": 2.55e18 } }
  ],
  "receivers": { "depth_km": 3, "locations": [[0, 4], [0, 6], [0, 8]] },
  "solver": { "tolerance": 1e-5 }
}
```

- `model`: an earth model in the deck format of Minos (six or nine columns).
  Relative paths, here and in `output.prefix`, are relative to the input
  file.
- `attenuation`: whether to include attenuation (default `true`).
- `gravity`: `"none"`, `"cowling"` or `"full"` (default `"full"`).
- `output.quantity`: `"displacement"`, `"velocity"` or `"acceleration"`
  (default displacement). `potential_and_tilt_corrections` applies to
  acceleration: it adds the change in gravity seen by an instrument that moves
  with the ground, and the tilt of its axis (default `false`).
- `degrees`: the spherical harmonic degrees summed (`min` defaults to 0).
- `frequencies_mHz`: the band computed; `time_series` the length of the
  seismograms and their sampling. The Nyquist frequency must lie above the
  highest frequency; without `time_step_s` the step is a quarter of the
  shortest period.
- `filter_mHz`: an optional raised-cosine taper, rising from zero to one
  between its first two corners and falling back between its last two.
- `sources`: one or more, each with its hypocentre (geographic latitude) and
  its moment tensor `moment_tensor_Nm`, a point force `force_N` (members
  `r`, `t`, `p`, each defaulting to zero), or both, in spherical coordinates
  at the source (r up, theta south, phi east), as a step at the origin time.
  Sources at one depth share the radial solves, so extra ones cost little. The optional `name`
  labels the output files when there are several sources (by default the
  sources are numbered from 1).
- `receivers`: a group of receivers, or an array of groups (see below). They
  may lie anywhere but at the centre: in the solid region of the sources, in
  an ocean or in the core. The sources must lie in one solid region of the
  model: in `prem.200`, below its 3 km ocean and above the core.
- `solver.tolerance`: the local error allowed per integration step, relative
  to the size of the solution (default `1e-5`; see below), and
  `solver.start_decay`, the decay in e-folds below a turning point at which an
  integration starts (default 11).

A member not listed here is an error, so that a misspelt one is not ignored.

Each group of receivers has its depth, or several, and its places: either a
list of `locations`, or a `profile` along a great circle, from one point `to`
another or from a point at an azimuth (clockwise from north) between two
distances, with `count` points equally spaced, both ends included:

```json
"receivers": [
  { "depth_km": 3, "locations": [[34.95, -106.46]], "instrument": "ANMO.pz" },
  { "depth_km": [3, 100, 200],
    "profile": { "from": [36.87, 35.58], "azimuth_deg": 60,
                 "distances_deg": [10, 90], "count": 41 } }
]
```

Receivers are numbered group by group and, within a group, depth by depth,
with every place at each depth. Profiles are great circles of the spherical
model, so a profile that starts at a source has points at exactly the
epicentral distances asked for. The optional `instrument` is a SAC pole-zero
file, applied to every component at the group's receivers; its output (counts,
for the files of the IRIS DMC) is then recorded. A comment `* INPUT UNIT : M`
(or `M/S`, or `M/S**2`) in the file says what motion the instrument responds
to, which is otherwise displacement, as in SAC; the conversion from the
quantity computed is exact. The response is applied at the same complex
frequencies as the rest of the calculation, which gives the exact causal
response.

### Snapshots

`yspec-snapshot` reads the same members as `yspec`, except that its `output`
has only `prefix` and `quantity`, and that `snapshots` replaces `receivers`:

```json
"snapshots": {
  "times_s": [300, 600, 900, 1200],
  "depths_km": { "min": 25, "max": 2875, "count": 58 },
  "distances_deg": { "min": 0, "max": 180, "count": 181 },
  "azimuths_deg": [30, 210]
}
```

The grid is in the frame of each source: depths, distances from it and
azimuths at it (clockwise from north). Each of the four is one number, an
array of them, or a range of `count` values from `min` to `max`, both
included. Two azimuths 180 degrees apart give a cross-section through the
source along a great circle; one depth and many azimuths, a map around it.
The points may lie anywhere but at the centre, through the oceans and the
core as well as the solid region of the sources.

One file is written for each time, `<prefix>.1`, `<prefix>.2`, ... (or
`<prefix>.<source>.1`, ... with several sources). After two comment lines,
each line holds a point: its depth (km), distance and azimuth (degrees), its
latitude and longitude, and the up, away and transverse components, away
being along the great circle from the source and transverse 90 degrees
anticlockwise from it, seen from above.

The cost grows with the numbers of depths and distances, while further
azimuths cost almost nothing: around the source only orders |m| <= 2 are
excited, so each depth and distance needs five sums over degree, which any
number of azimuths then share. The values at each time are summed directly
from the spectra, and agree with the samples of seismograms at the same
points to rounding error. The example above, at degrees to 300 and up to 15
mHz, takes about 3 s on 4 threads, and the same section through the whole
earth, from the sea surface to near the centre, about 4 s.

With gravity, the free surface of an ocean carries surface-gravity waves
(tsunamis), which are slow and outlast most time series: on and near the sea
surface they can wrap around in time, and appear early. A longer time series,
or points below the sea surface, avoid that.

## Accuracy and speed

The integration controls its step size, so its accuracy is set by
`solver.tolerance`. On the example input (degrees to 600, up to 40 mHz, 60
receivers), against seismograms from the Fortran with sixteen times its usual
step density, the largest error relative to the peak of each component is:

| tolerance       | error   | time, 16 threads |
|-----------------|---------|------------------|
| 1e-4            | 1.1e-3  | 4.9 s            |
| 1e-5 (default)  | 2.2e-4  | 8.0 s            |
| 1e-6            | 2.9e-5  | 12.3 s           |
| Fortran YSpec   | 2.3e-3  | 5.9 s            |

(on an 8-core laptop).

## How it works

For each frequency (complex, ω − iε, so that the response is damped and its
time series does not wrap around) and each degree, the radial
equations are solved for unit jumps in displacement and traction at the source
radius. The source is rotated to the north pole, where only orders |m| ≤ 2
are excited; the moment tensor and the receivers enter only
afterwards, through the weights of the jumps and the angular functions.

- **Model** (`YSpec/Model.hpp`, `src/Model.cpp`). The deck is read with
  planetmodel and made non-dimensional. Each property is a natural cubic
  spline within each layer, and gravity is the exact integral of the density
  spline.
- **Equations** (`src/Equations.hpp`). Toroidal and radial systems are
  2x2. Spheroidal systems are 4x4 without self-gravitation and 6x6 with it,
  integrated as minor vectors (second-order minors of the 4x4 systems, the
  fourteen independent third-order minors of the 6x6), with the source carried
  to the receiver by the b-vector system. Solid and fluid regions have their
  own forms, joined at every boundary between them.
- **Solver** (`YSpec/RadialSolver.hpp`, `src/RadialSolver.cpp`). One
  two-point solver serves the toroidal, radial and spheroidal systems, all
  three treatments of gravity and any sequence of regions. The integrations
  from below and from above are shared by every source and receiver depth:
  each is kept as it passes them, and each source then needs only its own
  determinant and the b-vectors carried out to the receivers, which change
  form with the minors at each boundary between solid and fluid, so that a
  receiver may be in an ocean or the core. Integration starts where a WKB
  estimate says the solution has decayed by e^-11 below its turning point,
  beneath the deepest source, from starting values for a homogeneous sphere;
  below that start the response is negligible and is taken to be zero.
- **Integration** (`src/RungeKutta.hpp`). A Dormand-Prince 5(4) method
  with step-size control, restarted at every discontinuity. Solutions are
  carried with a running logarithmic scale, since they grow and decay by far
  more than a double can hold.
- **Receivers** (`src/Geometry.hpp`, `src/AngularFunctions.cpp`). The angular
  functions come from the generalised Legendre functions of GSHTrans, which
  are regular at both poles.
- **Time series** (`src/TimeSeries.cpp`, `src/Synthetics.cpp`). Spectra
  are rotated to north and east, passed through any instrument, filtered,
  transformed with FFTW and undamped.
- **Snapshots** (`src/Snapshots.cpp`). The coefficients of each degree and
  order are gathered at each frequency, summed over degree at each depth and
  distance a block of distances at a time, and summed over frequency directly
  at each time.
- **Kernels** (`YSpec/Kernels.hpp`, `src/Kernels.cpp`). The forward and
  adjoint fields come from one solver, with the source and the receiver as its
  two source radii and the nodes of the radial mesh as its receivers. The
  adjoint field is expanded in the source's frame, where only orders |m| ≤ 2
  meet the forward field, so the integrals over each sphere are sums over
  degree and order of products of radial functions. The density's effect on
  gravity is a radial integral on the same nodes.

The calculation is threaded over frequency, with a loop over degree inside
each thread. The library's public headers, in `YSpec/`, use only the standard
library; the numerics and the libraries they build on are in `src/`.

## Using it as a library

Link the CMake target `YSpec::YSpec`. The public headers in `YSpec/` use only
the standard library, and the library works in SI units.

```cpp
#include <YSpec/YSpec.hpp>

const auto model = YSpec::Model::ReadDeck("prem.200");

auto calculation = YSpec::Calculation{};
calculation.degrees = {.min = 0, .max = 200};
calculation.band = {.min = 0.2e-3, .max = 15e-3};           // Hz
calculation.sampling = {.duration = 3600, .timeStep = 20};  // s

const auto sources = std::vector<YSpec::Source>{
    {.location = {0, 0}, .depth = 29.5e3, .moment = {.rr = 1e18}}};
const auto receivers = std::vector<YSpec::Receiver>{
    {.location = {0, 30}, .depth = 3e3}, {.location = {0, 60}, .depth = 3e3}};

const auto seismograms =
    YSpec::ComputeSeismograms(model, calculation, sources, receivers);
// The vertical component from the first source at the first receiver.
const auto z = seismograms.Trace(0, 0, YSpec::Vertical);
```

Sources may be at any depths in one solid region of the model (below an
ocean, for instance, and above the core), and receivers anywhere but at the
centre. The radial integrations are shared by all of them: several sources,
or receivers at several depths, cost little more than one. A source may be a
moment tensor, a point force (`YSpec::Force`) or both, and
`YSpec::ComputeStrains` gives seismograms of the strain tensor instead of the
displacement, as strainmeters and fibre-optic cables record; strain and
displacement satisfy reciprocity with each other to 1e-8, which the tests
check. A receiver may carry an `Instrument`, read from a
SAC pole-zero file with `YSpec::ReadSacPoleZero`, and `YSpec::GreatCircle`
gives lines of receivers along great circles.

`YSpec::ComputeRadialKernels` gives the sensitivity of a measurement on a
seismogram to spherically symmetric perturbations of the model, by the adjoint
method: of the density (with the gravity that follows from it), of the five
Love moduli and of their quality factors. A measurement is given by its
adjoint source, and `YSpec::Traveltime` forms that of a cross-correlation
traveltime in a window. The adjoint field is the forward response, at the same
complex frequencies, to a force at the receiver, so the kernels cost about
what two sets of seismograms at every node of a radial mesh would; they
predict finite differences of the measurement to about 1e-6, which the tests
check. The derivation is in `docs/theory/adjoint-kernels.tex`.

`YSpec::ComputeSnapshots` gives the wavefield on a `YSpec::SnapshotGrid` at
chosen times. `YSpec::RunSeismograms("input.json")` and
`YSpec::RunSnapshots("snapshot.json")` do what the programs do: read an input
file, compute what it describes and write it. `YSpec::RadialSolver` gives the
responses to unit jumps at a set of source radii, at a set of receiver radii,
for one degree and frequency.

## Differences from the Fortran YSpec

- Input is JSON instead of a parameter file read by counting lines.
- The integration controls its error, and is about ten times more accurate at
  the default tolerance.
- Several faults of the original are not reproduced: with gravity switched
  off the degree-zero equations no longer keep a gravity term; an ocean is no
  longer ignored at degree zero when the receiver is below the source; the
  horizontal components at the epicentre no longer have the wrong sign; and
  gravity near the centre is exact.
- Any sequence of solid and fluid regions is handled, including a solid
  region that reaches the centre, and a source or receiver on a solid-fluid
  boundary is placed on its solid side.
- The free-air and tilt corrections use gravity at the receiver rather than at
  the model knot below it, and with full self-gravitation the potential
  correction includes the gravity-flux term that the original omits; it
  vanishes at a free surface but not at a buried receiver (about 1e-3 of the
  peak acceleration on the sea floor).
- Receivers, and the points of snapshots, may be anywhere: in an ocean, in the
  outer core or in the inner core, as well as in the solid region of the
  sources. Their responses agree with direct shooting to 1e-8.
- Without self-gravitation and in the Cowling approximation, the starting
  solution of the S type lacked a factor of zeta on the Bessel-function term
  of its tangential displacement. Its error decays upwards from the start, so
  it only mattered for receivers near the start, which the original could not
  have.
- Degrees above 1827 are allowed: the angular functions are then tabulated in
  long double.

## Tests

`ctest --test-dir build` runs unit tests; comparisons of the solvers with the
Fortran's for single degrees and frequencies over a range of models and
options; checks independent of the Fortran (the normal modes of PREM, the
toroidal modes of a homogeneous sphere against their closed form, responses
in the oceans and the core against direct shooting of the solutions,
continuity at the epicentre, snapshots against seismograms, strain against
displacement by differences and reciprocity, kernels against finite
differences of their measurement in every treatment of gravity); and four
end-to-end cases compared with the Fortran's seismograms.
The scripts in `tests/fortran` regenerate the reference data from a checkout
of the Fortran YSpec.

## References

- Friederich, W. & Dalkolmo, J., 1995. Complete synthetic seismograms for a
  spherically symmetric earth by a numerical computation of the Green's
  function in the frequency domain. Geophys. J. Int., 122, 537-550.
- Al-Attar, D. & Woodhouse, J. H., 2008. Calculation of seismic displacement
  fields in self-gravitating earth models -- applications of minors vectors
  and symplectic structure. Geophys. J. Int., 175, 1176-1208.

## Licence

BSD 3-Clause; see `LICENSE`.
