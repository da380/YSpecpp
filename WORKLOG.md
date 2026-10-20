# YSpecpp work log

Running notes from the overnight session that started 2026-10-19. Nothing has
been committed.

## Summary for the morning

- `yspec input.json` works end to end. It reproduces the Fortran: per solve
  to 1e-6 or better, and seismograms to ~2e-4 of the peak against an accurate
  (16x step density) Fortran, while the shipped Fortran is off by ~2e-3. Speed
  at equal accuracy is a little better than the Fortran's.
- 54 tests, also clean under ASan/UBSan, -Werror, the serial build and a
  Doxygen build with warnings as errors.
- Two checks need no Fortran: the toroidal modes of a homogeneous sphere
  (closed form, 1e-5) and the normal modes of PREM.
- An independent review found two real bugs that no Fortran comparison could
  catch, both fixed and tested: a toroidal start at the centre of a solid
  region, and the sign of N/E at the epicentre (inherited from the Fortran).
- **Please look at:**
  1. The confirmed Fortran bugs listed under "Open questions"; two of them
     change l = 0 results.
  2. The change to the potential correction for buried receivers (derivation
     under "Independent review").
  3. The licence: there is no LICENSE file yet.
- `examples/input.json` is the Fortran example converted;
  `tools/legacy_to_json.py` converts others.

## Suggested next steps (need your decisions)

1. Several moment tensors per run. The solves do not depend on the moment
   tensor, so e.g. the six elementary tensors cost hardly more than one. It
   needs an input form (an array under `moment_tensor_Nm`?) and output naming
   (`prefix.<source>.<receiver>`?).
2. Receivers at several depths (one radial solve per depth, or several
   outputs from one integration).
3. Output: a header line, SAC, or one file per run.
4. Releases: pin planetmodel to a tag once one contains `cpp/`; with an
   install/export in planetmodel, YSpec could export a CMake package too.
5. GSHTrans: relax `MaxSafeDegree` when mMax is small (see above).
6. Licence (BSD-3 like the other libraries?).

## Status

The program works end to end and reproduces the Fortran.

- [x] Scaffold: CMake (dependencies pinned, fetched), GoogleTest, .clang-format
- [x] Model: deck -> layers, regions, natural splines, exact gravity
- [x] Bessel ratio q_l(z); Runge-Kutta integrator (Dormand-Prince 5(4), Cash-Karp)
      with step control and a running log-scale
- [x] Toroidal, radial and spheroidal solvers; one generic spheroidal solver
      over (None, Cowling, Full) x any sequence of solid/fluid regions
- [x] Source vectors; receiver geometry; angular functions from GSHTrans
      (no 1/sin(theta) anywhere); frequency grid; FFT (FFTWpp); filters
- [x] JSON input; `yspec` program; `tools/legacy_to_json.py` for old inputs
- [x] Tests: 54 (unit, per-(l, omega) against the Fortran solvers, normal
      modes of PREM, toroidal modes of a homogeneous sphere, epicentre
      continuity, the Fortran's four end-to-end cases, the README example)
- [x] README, CMake presets (incl. `local` against ~/dev checkouts), CI
      workflow (GCC 13/14, Clang 18, Debug/Release, serial, clang-format)

`ctest --test-dir build` runs everything in about two seconds.

## How the C++ compares with the Fortran

References are the Fortran at a3ae638 with its step density raised 16x
(nw = 64), which makes it accurate to ~1e-6; differences are relative to the
peak of each component.

- Per (l, omega), responses to unit jumps (tests/TestSolver.cpp, 11
  configurations: all gravity options, prem.200/.fic/.no, receiver above and
  below the source, attenuation on/off, l = 0..600, 0.4-40 mHz): agreement to
  1e-9 (toroidal), 3e-8..2e-6 (spheroidal), 1e-7 (radial).
- End to end (tests/regression): the four Fortran test cases agree to
  0.9-1.8e-4 at the default tolerance (with the two Fortran l = 0 bugs below
  patched in the reference). The shipped Fortran itself is off by 1.8e-3 on
  the same cases.
- Shipped example (lmax 600, 40 mHz, 60 receivers), 16 threads on this
  8-core laptop:

  | run                     | error vs accurate Fortran | time  |
  |-------------------------|---------------------------|-------|
  | Fortran as shipped      | 2.3e-3                    | 5.9 s |
  | C++, tolerance 1e-4     | 1.1e-3                    | 4.9 s |
  | C++, tolerance 3e-5     | 4.4e-4                    | 6.5 s |
  | C++, tolerance 1e-5 (default) | 2.2e-4              | 8.0 s |
  | C++, tolerance 1e-6     | 2.9e-5                    | 12.3 s|

  Serial: C++ 21.5 s at 1e-4 (Fortran 22.8 s per the handoff). Scaling is
  ~4.4x on 8 cores + SMT for both codes.
- Performance notes. All the time is in the solves (spheroidal ~90%); the
  angular sums, FFT and output are negligible. One right-hand side of the
  14-minor system costs ~70 ns (material 13, coefficient block 30, minors
  40). Thread scaling: 4.1x on 8 threads, no gain from SMT (the solver is
  floating-point bound; a latency-bound toy loop gets 5.4x/9x here), so this
  laptop caps it. Tried without effect: handing out high frequencies first,
  per-frequency local accumulation (kept: avoids false sharing in principle),
  degree-major angular table (kept). Real-valued entries in the coefficient
  blocks gave ~5%.
- The four b-vectors are integrated from source to receiver together
  (`RungeKutta::IntegrateTogether`: one step sequence, the coefficients formed
  once per stage, each vector with its own error control and scale). A 600 km
  source now costs 2.5x a 29.5 km one, against 3.6x before.
- Future speed-ups, not done: (1) solve several degrees together on one
  radial grid to share the material evaluation and vectorise; (2) a
  per-frequency cut-off in degree (risky: shallow sources have slowly
  decaying near-field terms).
- A residual ~8e-5 does not shrink with tolerance or start depth; it comes
  from the Fortran's gravity near the centre (below).
- Independent of the Fortran (tests/TestNormalModes.cpp): the responses for
  prem.200 peak at 0T2 = 0.37916 mHz (PREM 0.37917), 0S0 = 0.81403
  (0.8143), 0S2 = 0.30891 (0.30928), 0S3 = 0.46812 (0.46856), from memory of
  the Minos values. Spheroidal modes come out ~0.1% low; the Fortran gives the
  same (it agrees per solve to 1e-8), so this is YSpec's physics, perhaps
  exact versus first-order dispersion, or my reference values. Without
  gravity 0S2 is at 0.242 mHz, and the Cowling approximation barely moves it.
- The whole suite passes under AddressSanitizer and UBSan (preset `asan`) and
  with -Wall -Wextra -Werror.

## Independent review (a subagent, read-only, against the Fortran)

Found and fixed, each with a test that fails without the fix:
- A toroidal integration starting at the centre of a solid region walked from
  r = 0 (where the system is singular) and threw: any source in a solid inner
  core, or a model whose solid region reaches the centre, at low l. The
  Fortran never starts below the CMB, so no Fortran test could see it. New
  analytic test: the toroidal modes of a homogeneous solid sphere agree with
  the roots of (l-1) j_l(x) = x j_{l+1}(x) to 1e-5.
- North and east had the wrong sign at the epicentre (inherited from the
  Fortran's delaz, which returns zero azimuths there). Geometry now uses
  atan2 (azimuths exact near 0 and pi, where acos lost 1e-8) and the limits of
  a receiver approaching the source. Test: seismograms at the epicentre equal
  those 0.001 degrees away in every direction.
Also: solver options and input values are validated (tolerance in (0, 1),
positive start decay, latitudes, non-negative depths, depths inside the
model, a non-empty filter pass band); a source or receiver exactly on a
solid-fluid boundary is put on the solid side (the Fortran refuses a source
at the CMB depth), and a zero solution (some jumps there have none) no longer
stalls the integrator; the integrator no longer restarts each layer with the
tiny step that ended the previous one; `legacy_to_json.py` counts lines as
the Fortran does.

Changed, **please confirm**: the potential correction for acceleration at a
buried receiver. The Fortran adds -((l+1) phi + 2 g U)/r (in its r-scaled
variables), which is r (phi' + U g') only where the flux variable psi
vanishes, i.e. at a free surface. The variables are y1 = r U (the Fortran
divides by r for the displacement), y3 = r phi and
y6 = psi = r (phi' + (l+1) phi/r + 4 pi G rho U)/(4 pi G): with y3 = r phi and
psi = 0, the 6x6 system's phi equation is exactly the condition for matching
the exterior potential, and the Fortran's surface start imposes
P = S = psi = 0. Then r (phi' + U g') = 4 pi G psi - ((l+1) phi + 2 g U)/r
anywhere, so the C++ adds 4 pi G psi. A test checks psi = 0 at a free
surface. Effect: up to 1.3e-3 of the peak of the vertical acceleration on the
sea floor (3 km) and 2.8e-3 at 50 km depth; nothing at a free surface, none
of the regression cases changes (none has full gravity, corrections and l > 0
together). Also noted, not changed: in the Cowling approximation the
correction is -2 g U / r, which omits 4 pi G rho U of g' at a buried receiver;
in the full case that term cancels against phi'.

## Reference

- The reference is the cleaned Fortran at ../YSpec commit a3ae638 (develop,
  "modernisation and clean up of the code"; David is pushing it to main). It
  is exported with `git archive` into the session's scratch directory.
- `tests/fortran/` holds the driver and the scripts that regenerate
  `tests/data/solver` and `tests/regression/*/reference` from a YSpec checkout.

## Decisions

- House style of GSHTrans; header-only library in `YSpec/`, program in `app/`.
- JSON input via planetmodel's own reader (`planetmodel/json.hpp`), so no new
  dependency. It is strict RFC 8259: no comments in input files. Unknown
  members are errors. Relative paths are relative to the input file.
- G = 6.6723e-11 (Minos and the Fortran) and the Fortran's normalisation
  (pi G = 1), so that intermediate values compare directly.
- Natural cubic splines per layer, as in the Fortran. Gravity is integrated
  exactly from the density splines.
- planetmodel pinned to a commit on main (fa6e230), netCDF off. GSHTrans
  pinned to 8dd116d, BLAS off.
- The loop over frequency is the threaded one, with degree inside it: each
  thread adds into its own frequency of the receiver spectra, so there is no
  shared mutable state and no barrier per degree (the Fortran threads over
  frequency inside a loop over degree).
- Start level: WKB decay criterion only (no PREM-specific thresholds), 11
  e-folds as in the Fortran; reset at each solid/fluid boundary. Toroidal
  starts are confined to the source's solid region.
- Free-air/tilt corrections use g at the receiver radius (the Fortran uses g
  at the knot below it: ~1e-4 effect on the acceleration cases).
- Tried and dropped: controlling only the error perpendicular to the solution
  for the long integrations (whose overall scale cancels). Fewer steps, but
  only ~20% faster at equal accuracy and less predictable, because the norm of
  a minor vector oscillates too.

## Open questions / follow-ups for David

- GSHTrans `Wigner` refuses lmax > 1827 in double. As you suspected it is a
  large-mMax issue: `MaxSafeDegree` guards the underflow of the seed of the
  boundary orders m = +-l, which a table with mMax = 2 never forms. So the
  check could depend on min(lMax, mMax). Meanwhile YSpec builds its table in
  long double above 1827 (limit 30747); a test checks it to l = 2500.
- **Fortran bug (confirmed):** `rad_derivs_ng` calls `rad_amat`, which
  includes -4 rho g / r, so with gravity switched off l = 0 still has gravity
  (while `rad_ystart_ng` starts from a gravity-free solution). Changes the
  l = 0 response by up to ~20% at 0.5 mHz. The C++ omits the term.
- **Fortran bug (confirmed):** `rad_solver_ng`, receiver-below-source branch,
  starts the downward integration at the sea floor (`r(nsl)`) rather than the
  surface, so the ocean is ignored for l = 0 there (~10% at 6 mHz); the other
  branch, and `rad_solver_wg`, use the surface. Affects the `fic_nograv_vel`
  regression case (~1e-3 of the peak of Z).
- **Fortran inaccuracy:** gravity near the centre. The Fortran splines
  rho r^2 with natural end conditions, forcing zero curvature at r = 0 where
  the true curvature is 2 rho: g is 15% off at the first knot above the
  centre, 5e-6 in mid inner core, 1e-7 in the mantle.
- Fortran fragility (harmless in yspec.f90): `rad_solver_ng`/`_wg` call
  `set_omega` before `set_ats`, so the dispersion uses whatever `ats` the
  thread last set. yspec.f90 always calls `source_vector_rad` first, which
  sets it.
- Licence: I have not added a LICENSE file; your other libraries are BSD-3.
