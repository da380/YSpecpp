# Handoff: notes from the review of the Fortran YSpec

Written 2026-10-19 by the Claude session that reviewed and is cleaning up the
original Fortran code in `../YSpec`. It records what that session learned so
that the C++ version does not have to rediscover it. Treat it as background:
David's instructions take precedence.

## Goal (as agreed with David)

- A standalone, modern C++ re-implementation of YSpec in this repository.
- Functionally equivalent in *intent*, not in execution: same physics, inputs
  and outputs, but free to use better numerics.
- Keep the YSpec name visible: namespace `YSpec`, CMake target `YSpec::YSpec`,
  an executable called `yspec`, ideally able to read the legacy parameter file.
- No Numerical Recipes code (its licence does not allow redistribution of the
  source). No Boost. Write a small local ODE integrator instead.
- OpenMP from the start: every (l, ω) solve is independent.
- David makes all git commits; do not commit unless asked. Refer to David's
  student Alex (author of `../DSpecM1D`) as they/them.

## The method in one page

Synthetic seismograms for a spherically symmetric, transversely isotropic,
anelastic, (optionally) self-gravitating earth by direct radial integration
(Friederich & Dalkolmo 1995; Al-Attar & Woodhouse 2008):

- Work at complex frequency ω − iε with ε = mex/T (mex = 5, T the length of the
  time series); undo the damping with exp(εt) after the inverse FFT. See
  `fcal` in `src/module_fourier.f90`.
- Rotate so the source is at the pole: only orders |m| ≤ 2 are excited, and the
  receiver sum uses the epicentral distance and azimuth.
- For each degree l and frequency ω, solve the radial ODEs for the response to
  *unit* source jumps (the identity matrices `s2tmp`, `s4tmp`, `s6tmp` in
  `src/yspec.f90`). The solutions do not depend on the moment tensor, which is
  contracted afterwards. This makes many sources at one depth, or all six
  elementary moment tensors, almost free.
- Toroidal: 2×2 system. Radial (l = 0): 2×2 (no gravity) or with the potential.
  Spheroidal: 4×4 without gravity or with Cowling, 6×6 with full
  self-gravitation, integrated as **minor vectors** (5 components for 4×4,
  14 reduced components for 6×6), with solid↔fluid maps at the ICB, CMB and
  sea floor, and the inhomogeneous ("b-vector") system carried from the source
  to the receiver. See `src/module_int.f90` and `src/module_solver.f90`.
- Attenuation: exact constant-Q log dispersion at complex ω,
  M(ω) = M0 [1 + q (2/π) ln(iω/ω0)], with Minos-style effective Q for C, F.
- Starting radius: a WKB estimate of where the solution has decayed by e^−11
  below its turning point (`src/module_start.f90`); start vectors from
  homogeneous-sphere spherical Bessel solutions.
- Source time function: a Heaviside step in moment (the 1/(iω) in the source
  vectors). Output: Z, N, E in SI units after cosine band-pass filtering.

## Using the Fortran as the reference

- Use the **cleaned Fortran** in `../YSpec` (branch `develop`, being worked on
  and not yet committed as of 2026-10-19), never the original release: the
  original has the bugs listed below.
- Build: `cmake -S ../YSpec -B build-yspec -DCMAKE_BUILD_TYPE=Release && cmake
  --build build-yspec && ctest --test-dir build-yspec`. The cleaned version
  will need FFTW3 (it replaces the Numerical Recipes FFT).
- `../YSpec/tests/` has four small regression cases (input files plus reference
  outputs) that between them cover all three gravity options, solid and fluid
  inner core, displacement/velocity/acceleration, potential and tilt
  corrections, receivers at both poles, a receiver below the source, and an
  l = 0-only case. They run in about a second and make good first targets.
- Parameter-file quirks to support in a legacy reader: values are found by
  counting lines (two header lines, then each value preceded by exactly one
  blank line and one comment line); unquoted paths are truncated at `/`.
- Output files `prefix.N`: columns `t  Z  N  E`, format `5e15.6`.
- Latitudes are geographic, converted to geocentric with flattening about
  1/297 (`delaz` in `src/module_function.f90`); points exactly at a pole are
  moved 0.02° off it.

## Bugs found in the original Fortran (do not reproduce)

All fixed in the cleaned Fortran; sizes measured on PREM test cases.

1. Toroidal u_θ used X_lm instead of X_lm/sinθ (`x` for `xc`): radial-horizontal
   component wrong by 0% at 90° up to ~100% at 5° and ~40% at 175°.
2. l = 1 toroidal modes skipped (0.1–1% effect).
3. Spheroidal source moduli lost their imaginary (1/Q) part (declared real):
   ~0.03% for crustal sources, ~1% in the low-Q zone.
4. No free-air correction for l = 0 with potential/tilt corrections on. For
   l = 0 the correction is exactly −2gU/r: the 4πGρU terms from ∂rφ1 and from
   the gravity gradient cancel, so the (gauge-dependent, trapezoid-rule)
   l = 0 potential computed in `rad_solver_wg` is not needed.
5. Receiver at exactly 90°N placed at the south pole (truncated deg→rad).
6. Double-precision π etc. were written with single-precision literals
   (π accurate to 3e-8); several `sqrt(real(...))` in single precision.
7. Uninitialised arrays (harmonics for l = 1, |m| = 2), uninitialised function
   result in `find_radius_index`, out-of-bounds write for `lmax = 0`.
8. Spherical Bessel starting values: the `sbessj2` "ratio" is only the
   small-argument limit and is off by one in its order for its callers; `bfs`
   uses sqrt(abs(x²)) (wrong phase for complex arguments) and the convergence
   test |abs(Δ)−1| < ε, which is wrong for complex Δ. These turn out to be
   harmless in practice (see below), but a clean version needs only one
   function: q_l(z) = x j_{l+1}(x)/j_l(x) with z = x², which satisfies
   q_l(z) = z/(2l+3 − q_{l+1}(z)). It is analytic in z, so it needs no square
   roots or branch choices and works for complex k²r².
9. The Minos frequency thresholds in `sph_start_level` (aw, bw, dw) are
   PREM-specific shortcuts; the WKB decay criterion alone is model-independent
   (the cleaned Fortran tries the mantle, then the outer core, then the inner
   core).
10. The original FFT (Numerical Recipes `fourrow`) built its twiddle factors
   from the single-precision π, an error of ~1e-6 of peak in the seismograms.
11. The number of output samples depended on rounding when the series length
   is a multiple of the time step; the cleaned code writes 0 ≤ t ≤ tout.

## Measured facts worth knowing

- Accuracy of the Fortran integrator: fixed-step Cash–Karp RK5 with step
  α/(4ω) and the error estimate discarded. Doubling the step density changes
  seismograms by ~1e-3 of peak, so that is the reference's accuracy. A C++
  version with error control should beat it; compare with tolerances to match.
- Start level (lmax 600, 40 mHz example): compared with always starting at the
  bottom (3× slower), an e^−11 cut-off differs by ~1e-6 of peak, e^−15 by
  ~5e-9 (+4% run time) and e^−20 by ~1e-9 (+10%). The cleaned Fortran uses
  e^−15. All of these are far below the integration error.
- OpenMP over frequency (per l) gives output bitwise identical to serial;
  the shipped example (lmax 600, 40 mHz, 60 receivers) went 22.8 s → 6.5 s on
  16 cores, with the remaining time in the serial FFT and output.
- `../DSpecM1D` (Alex's spectral-element code, GPL-3) is an independent check:
  its published comparisons with (Alex's fixed) YSpec agree well below 1%.
  It evaluates receivers with GSHTrans canonical components. Comparing raw
  time series between the codes is confounded by different filter/taper
  chains; compare spectra, or pass both through the same filter.

## Suggested building blocks (David's libraries in `../`)

All header-only C++20/23, BSD-3, CMake targets `Name::Name`, consumed with
`FetchContent_Declare(... FIND_PACKAGE_ARGS)`; look at them for house style.

- **NumericConcepts**: concepts (`Real`, `Complex`, `RealOrComplex`, ranges).
- **Interpolation**: `CubicSpline` (natural/clamped/not-a-knot, real or complex
  ordinates, any derivative), `CubicSplineSystem` (factor a grid once, solve
  many profiles), `Piecewise` + `SplitAtRepeats` for the deck's repeated radii,
  `Primitive(...).Integral(a, b)` for the gravity integral.
- **planetmodel** (C++ part in `cpp/`): `ReadDeck(path, MINEOS())` reads the
  YSpec/Minos deck format (6 or 9 columns) with layers and named boundaries.
  Its C++ side has no dispersion (the Python `rheology.py` does).
- **GSHTrans**: `Wigner<...>` evaluates sqrt((2l+1)/4π) d^l_{Nm}(θ) for all l at
  one θ. Using N ∈ {−1, 0, +1} (and ±2 for strain) gives the receiver terms in
  canonical components with no 1/sinθ, no pole special cases, and no room for
  bug 1. Check the convention table in `docs/gshtrans-reference.pdf`.
- **FFTWpp**: batched c2r transforms over receivers × components.
- **GaussQuad**: Gauss–Legendre rules if needed.
- Not covered by any library (write locally): the complex ODE integrator, the
  Bessel ratio q_l(z), the minor-vector systems and interface maps.

## Integrator suggestion

A small header templated on a Butcher tableau, so that Dormand–Prince 5(4) and
DOP853 can be compared; step-size control on a scaled norm; restart at every
model discontinuity; track the solution's scale as a running log-norm instead
of the Fortran's 1e±20 rescaling counters. The coefficients are smooth within
layers and the systems are small (2–14 complex components), so fixed-size
`std::array` state is natural.

## Suggested order of work

1. Model: deck → per-layer splines of ρ, A, C, F, L, N, Qκ, Qμ; gravity;
   complex moduli at complex ω.
2. Toroidal end to end for one (l, ω), checked against the Fortran.
3. Radial, then spheroidal (no gravity → Cowling → full; fluid core, ocean).
4. Source jumps, receiver terms (GSHTrans), frequency grid, FFT, filters,
   output; legacy input reader; OpenMP.
5. Regression tests against the cleaned Fortran's `tests/` cases.

Fortran routines to read first: `src/yspec.f90` (driver), `tor_solver` and
`source_vector_tor` in `src/module_solver.f90`, `yint` and the `*_amat`/
`*_derivs` routines in `src/module_int.f90`, `start_level` in
`src/module_start.f90`, `read_model` in `src/module_model.f90`.
