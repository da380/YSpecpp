# yspec-seismo

Synthetic seismograms, snapshots of the wavefield and adjoint sensitivity
kernels for spherically symmetric earth models, from Python. It is the
Python package of [YSpec](https://github.com/da380/YSpecpp), which integrates
the equations of motion of a self-gravitating, anelastic, transversely
isotropic earth directly in radius, in the frequency domain, for every degree
and frequency, and sums the responses: complete seismograms, with every body
wave, surface wave and normal mode the band holds.

```sh
pip install yspec-seismo            # the library and the programs
pip install "yspec-seismo[plot]"    # and yspec-plot, to draw what they write
```

It is imported as `yspec`:

```python
import numpy as np
import yspec

model = yspec.Model.prem()
source = yspec.Source(latitude=36.87, longitude=35.58, depth_km=29.5,
                      moment_tensor=dict(rr=4.22e17, rt=3.08e17, rp=6.99e17,
                                         tt=-2.976e18, tp=6.54e17, pp=2.55e18))
receivers = yspec.great_circle(source, azimuth_deg=60,
                               distances_deg=np.arange(5, 180, 5))
calculation = yspec.Calculation(degrees=300, frequencies_mhz=(0.2, 15),
                                length_min=60, filter_mhz=(1, 2, 12, 15),
                                quantity="velocity")

seismograms = yspec.seismograms(model, calculation, source, receivers,
                                progress=True)
seismograms.time_s           # the time of each sample (s)
seismograms.data             # [source, receiver, component, sample]: Z, N, E
seismograms.rotated()        # the same with Z, R and T
seismograms.to_obspy()       # an ObsPy Stream (with ObsPy installed)
```

Lengths are in km, angles in degrees, frequencies in mHz and times in s, as
in the programs' input files; the results are numpy arrays in SI units.

- `yspec.seismograms` and `yspec.strains`: at receivers anywhere but the
  centre, in the ocean, the crust or the core, from moment tensors and forces,
  through seismometers given by their poles and zeros.
- `yspec.snapshots`: the wavefield at chosen times on a grid of depths,
  distances and azimuths around each source: cross-sections through the
  earth, or maps.
- `yspec.traveltimes`, `yspec.radial_kernels` and `yspec.sampled_kernels`:
  cross-correlation traveltimes and their sensitivity kernels by the adjoint
  method, to spherically symmetric perturbations of density, the elastic
  parameters and the quality factors, and to laterally varying ones of the
  elastic parameters and quality factors, through the earth, on sections
  about the path or on maps.
- The programs `yspec`, `yspec-snapshot` and `yspec-kernels`, which read
  JSON input files and write their results as text, and `yspec-plot`, which
  draws them (with the `plot` extra).

The calculations run on as many threads as OpenMP gives them
(`OMP_NUM_THREADS`, or `yspec.set_threads`), except on macOS, where the
wheels have one. Ctrl-C stops a calculation.

## References

- Al-Attar, D. & Woodhouse, J. H., 2008. Calculation of seismic displacement
  fields in self-gravitating earth models -- applications of minors vectors
  and symplectic structure. Geophys. J. Int., 175, 1176-1208.
- Friederich, W. & Dalkolmo, J., 1995. Complete synthetic seismograms for a
  spherically symmetric earth by a numerical computation of the Green's
  function in the frequency domain. Geophys. J. Int., 122, 537-550.

BSD 3-Clause licence.
