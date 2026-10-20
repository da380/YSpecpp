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

```sh
build/bin/yspec examples/input.json
```

The number of threads is set with `OMP_NUM_THREADS`. One file is written per
receiver, `<prefix>.1`, `<prefix>.2`, ..., in the order the receivers are
listed; with several sources, one per source and receiver,
`<prefix>.<source>.1`, .... Each line holds the time after the origin time (s)
and the vertical (up), north and east components, in SI units of the quantity
asked for.

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
  highest frequency.
- `filter_mHz`: an optional raised-cosine taper, rising from zero to one
  between its first two corners and falling back between its last two.
- `sources`: one or more, each with its hypocentre (geographic latitude) and
  its moment tensor in spherical coordinates at the source (r up, theta south,
  phi east), for a step in moment at the origin time. Sources at one depth
  share the radial solves, so extra ones cost little. The optional `name`
  labels the output files when there are several sources (by default the
  sources are numbered from 1).
- `receivers`: their depth and locations as `[latitude, longitude]` pairs.
  The sources and the receivers must lie in the same solid region of the
  model: in `prem.200`, below its 3 km ocean.
- `solver.tolerance`: the local error allowed per integration step, relative
  to the size of the solution (default `1e-5`; see below), and
  `solver.start_decay`, the decay in e-folds below a turning point at which an
  integration starts (default 11).

A member not listed here is an error, so that a misspelt one is not ignored.

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
  determinant and the b-vectors carried out to the receivers. Integration
  starts where a WKB estimate says the solution has decayed by e^-11 below its
  turning point, from starting values for a homogeneous sphere.
- **Integration** (`src/RungeKutta.hpp`). A Dormand-Prince 5(4) method
  with step-size control, restarted at every discontinuity. Solutions are
  carried with a running logarithmic scale, since they grow and decay by far
  more than a double can hold.
- **Receivers** (`src/Geometry.hpp`, `src/AngularFunctions.cpp`). The angular
  functions come from the generalised Legendre functions of GSHTrans, which
  are regular at both poles.
- **Time series** (`src/TimeSeries.cpp`, `src/Synthetics.cpp`). Spectra
  are filtered, transformed with FFTW, undamped, and rotated to north and east.

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

Sources and receivers may be at any depths in one solid region of the model
(below an ocean, for instance, and above the core). The radial integrations
are shared by all of them: several sources, or receivers at several depths
beneath the same places, as for a snapshot of the wavefield, cost little more
than one. `YSpec::Run("input.json")` does what the `yspec` program does: it
reads an input file, computes the seismograms and writes them.
`YSpec::RadialSolver` gives the responses to unit jumps at a set of source
radii, at a set of receiver radii, for one degree and frequency.

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
- Degrees above 1827 are allowed: the angular functions are then tabulated in
  long double.

## Tests

`ctest --test-dir build` runs unit tests; comparisons of the solvers with the
Fortran's for single degrees and frequencies over a range of models and
options; checks independent of the Fortran (the normal modes of PREM, the
toroidal modes of a homogeneous sphere against their closed form, continuity
at the epicentre); and four end-to-end cases compared with the Fortran's
seismograms.
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
