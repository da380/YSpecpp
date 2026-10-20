# First-seismogram illustration provenance

The illustration `../assets/first-seismogram.svg` was generated from this exact
fixture and output. It is an explanatory trace, not a scientific reference
result.

- Code revision: `9552d2bde73b11551027fb829b17d7b4813253d9` (`docs/phase1-user-guide`; documentation changes were uncommitted during the run).
- Platform/compiler: Ubuntu 24.04; GCC 13.3.0 (`g++-13`).
- Build: CMake 3.28.3, Release, `YSPEC_BUILD_TESTS=OFF`, OpenMP enabled; two build jobs. Five dependencies were selected from the existing read-only cache via explicit source-directory overrides; planetmodel was fetched at its pinned SHA. This was not a cold-network configure.
- FFTW: 3.3.10, with single-, double- and long-double components found.
- Dependency revisions: NumericConcepts `e03a2889d43a03660d0879f46dda1406d90c9517` (v1.0.0); Interpolation `d5a1711ca6a89fc8ad8148fe13212398fe569bc0` (v1.1.2); FFTWpp `cb152c0e85f3822493c2d9e85a56307ac7b4f1c0` (v1.0.1); GaussQuad `7d4786d192f68077cd99a34dcdf1379ed1eeac9a`; GSHTrans `8dd116df8df4222fe9b10a40c7213d34f2a274fb`; planetmodel `fa6e230fde25152d40f27d9839de44c4f3182a24`.
- Runtime command: `OMP_NUM_THREADS=2 ../guide/bin/yspec input.json`, from `build/first-seismogram/`.
- Measured solver runtime: 0.50 s. Progress reported 26 frequencies, one source and one receiver.
- Input: `first-seismogram.json`; SHA-256 `72a30885258026dc1662f7713ac3674a37985bf045165965d41dd308c1f219a4`.
- Output: `first-seismogram.out.1`; 91 rows, four columns, times 0–1800 s in 20 s steps; SHA-256 `606e7c0eece9336b83e3247c97587d5e560c871fef84c90de9d43fe9dc4b5ece`.
- Column meanings: time (s), up displacement (m), north displacement (m), east displacement (m). The vertical peak absolute amplitude is `4.344679e-6 m`; all values were finite and the trace was nonzero.
- Plot library: Matplotlib 3.6.3. The plot script sets a fixed SVG hash salt for reproducible element IDs and removes only trailing horizontal whitespace while preserving line endings. Plot SHA-256: `925f2609ad95bab8d498edcdca73a3ba1add470780bb8c267ff7f7e5878abe61`.
- Plot command from repository root: `python3 docs/site/examples/plot_first_seismogram.py build/first-seismogram/first-seismogram.out.1 build/first-seismogram/first-seismogram.svg`.
