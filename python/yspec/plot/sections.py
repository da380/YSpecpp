"""Pictures of three-dimensional kernels on sections and maps about the
path from the source to the receiver, with the seismogram measured: the
banana of finite-frequency theory along the path, and its doughnut across
it, drawn as planetmodel draws sections.

Near the source and the receiver the kernels are far larger than elsewhere;
the colours are scaled to the kernel away from them, at more than `clear`
degrees, and saturate there.
"""

import math

import numpy as np
from planetmodel.plotting import section
from planetmodel.slices import SectionSlice, Shown, SliceLayer

from .runs import KernelRun
from .seismograms import COLOURS
from .snapshots import GUIDE, STAR, _phase_boundaries
from .volume import TRIANGLE


def _parameter(run, parameter):
    """The kernel drawn: as asked, or vs for a phase with a leg in S, else
    vp."""
    if parameter is None:
        parameter = "vs" if "S" in run.name else "vp"
    if parameter not in run.parameters:
        raise ValueError(f"there is no kernel of {parameter}: there are "
                         f"{', '.join(run.parameters)}")
    return parameter


def _slice(sampled, coordinate, parameter, name, surface):
    """A section as planetmodel's slice, drawn in a plane of its own: the
    distance along it from its start (radians), and each layer's radii (km)
    by places."""
    angle = np.radians(coordinate - coordinate[0])
    path = np.stack([np.cos(angle), np.sin(angle), np.zeros_like(angle)],
                    axis=-1)
    layers = []
    for layer in sampled.layers:
        radius = np.broadcast_to((surface - layer.depth)[:, None],
                                 layer.values[parameter].shape)
        points = radius[..., None] * path[None, :, :]
        layers.append(SliceLayer(index=layer.index, radius=np.array(radius),
                                 values=layer.values[parameter],
                                 reference=points, physical=points))
    plane = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    return SectionSlice(distance=angle, path=path, plane=plane,
                        layers=tuple(layers), shown=Shown(name))


def _limit(run, sampled, parameter, keep, clip, shallowest=150.0):
    """The colour limit: the clip percentile of the kernel's size at the
    places kept, below the shallowest depth (km), above which the crust's
    large sensitivity would set it."""
    values = np.concatenate([
        layer.values[parameter][layer.depth > shallowest][:, keep].ravel()
        for layer in sampled.layers])
    values = np.abs(values[np.isfinite(values)])
    return float(np.percentile(values, clip)) if values.size else 1.0


def _deepest(run):
    """How deep a section is drawn: through the core for a phase with a leg
    in it, otherwise to just below the core-mantle boundary."""
    return 0.0 if "K" in run.name else 3000.0


def _draw(ax, run, sampled, coordinate, parameter, limit, boundaries,
          marks):
    """A section in its plane, with the boundaries between solid and fluid
    and the marks (angle along it, radius, style) of the source and
    receiver, cropped to the depths the phase samples."""
    label = f"K, {'ln ' if run.relative else ''}{parameter} (s m$^{{-3}}$)"
    cut = _slice(sampled, coordinate, parameter, label, run.surface)
    section(ax, cut, physical=True, diverging=True, vmin=-limit,
            vmax=limit, interfaces=False)
    half = 0.5 * cut.length
    deepest = _deepest(run)
    arc = np.linspace(-half, half, 721)
    for r in [run.surface] + [b for b in boundaries if b > deepest]:
        ax.plot(r * np.sin(arc), r * np.cos(arc), **GUIDE)
    for side in (-half, half):
        ax.plot([deepest * math.sin(side), run.surface * math.sin(side)],
                [deepest * math.cos(side), run.surface * math.cos(side)],
                **GUIDE)
    for angle, radius, style in marks:
        turn = math.radians(angle - coordinate[0]) - half
        ax.plot(radius * math.sin(turn), radius * math.cos(turn), **style)
    width = run.surface * math.sin(min(half, math.pi / 2))
    bottom = deepest * math.cos(half) if half < math.pi / 2 else -run.surface
    if deepest == 0:
        bottom = min(0.0, run.surface * math.cos(half))
    ax.set_xlim(-1.02 * width, 1.02 * width)
    ax.set_ylim(bottom - 0.02 * run.surface, 1.02 * run.surface)
    ax.set_aspect("equal")
    ax.set_axis_off()


def _trace(ax, run):
    """The seismogram measured, its window shaded."""
    t = run.time / 60
    peak = np.abs(run.seismogram).max() or 1.0
    inside = run.taper > 0
    ax.plot(t, run.seismogram / peak, lw=0.5, color="0.6")
    ax.plot(t[inside], run.seismogram[inside] / peak, lw=0.9,
            color=COLOURS.get(run.component, "k"))
    ax.fill_between(t, -run.taper, run.taper, color="0.9", zorder=0)
    ax.set_xlim(t[0], t[-1])
    ax.set_ylim(-1.05, 1.05)
    ax.set_yticks([])
    ax.set_xlabel("time (min)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.text(0.01, 0.92, f"{run.name} on {run.component}",
            transform=ax.transAxes, va="top", fontsize=9)


def _marks(run):
    surface = run.surface
    return ([(0.0, run.source_radius, STAR),
             (run.distance, surface - run.receiver_depth, TRIANGLE)],
            [])


def banana(run, parameter=None, *, clear=5.0, clip=99.0):
    """One measurement: its seismogram and window, and its kernel on the
    sections along and across the path."""
    import matplotlib.pyplot as plt

    parameter = _parameter(run, parameter)
    path, across = run.sections()
    boundaries, _ = _phase_boundaries(run.model)
    keep = (np.abs(path.along) > clear) & (np.abs(path.along - run.distance) >
                                           clear)
    limit = _limit(run, path, parameter, keep, clip)
    fig = plt.figure(figsize=(14, 9))
    grid = fig.add_gridspec(2, 2, height_ratios=[1, 4],
                            width_ratios=[2.2, 1])
    _trace(fig.add_subplot(grid[0, :]), run)
    on_path, _ = _marks(run)
    ax = fig.add_subplot(grid[1, 0])
    _draw(ax, run, path, path.along, parameter, limit, boundaries, on_path)
    ax.set_title("along the path", fontsize=9)
    ax = fig.add_subplot(grid[1, 1])
    _draw(ax, run, across, across.across, parameter,
          _limit(run, across, parameter, np.ones(len(across.across), bool),
                 clip), boundaries, [])
    ax.set_title("across it, at its middle", fontsize=9)
    measured = run.description.removeprefix("radial kernels of ")
    fig.suptitle(f"The kernel of {'ln ' if run.relative else ''}{parameter} "
                 f"for {measured}; colours saturate near the source (star) "
                 f"and the receiver (triangle)", fontsize=10)
    fig.tight_layout()
    return fig


def gallery(path, parameter=None, *, clear=5.0, clip=99.0):
    """Every measurement of a run, a row for each: its seismogram and
    window, and its kernel along the path and across it."""
    import matplotlib.pyplot as plt

    names = KernelRun(path).names
    fig = plt.figure(figsize=(16, 3.0 * len(names)))
    grid = fig.add_gridspec(len(names), 3, width_ratios=[1.1, 2.0, 0.8])
    for row, name in enumerate(names):
        run = KernelRun(path, name)
        chosen = _parameter(run, parameter)
        along, across = run.sections()
        boundaries, _ = _phase_boundaries(run.model)
        keep = (np.abs(along.along) > clear) & (
            np.abs(along.along - run.distance) > clear)
        _trace(fig.add_subplot(grid[row, 0]), run)
        on_path, _ = _marks(run)
        ax = fig.add_subplot(grid[row, 1])
        _draw(ax, run, along, along.along, chosen,
              _limit(run, along, chosen, keep, clip), boundaries, on_path)
        ax.set_title(f"{run.name}: {'ln ' if run.relative else ''}{chosen}, "
                     f"along the path", fontsize=9)
        ax = fig.add_subplot(grid[row, 2])
        _draw(ax, run, across, across.across, chosen,
              _limit(run, across, chosen, np.ones(len(across.across), bool),
                     clip), boundaries, [])
        ax.set_title("across it, at its middle", fontsize=9)
    fig.suptitle("Traveltime kernels by cross-correlation in the windows "
                 "shaded; colours saturate near the source (star) and "
                 "receiver (triangle)", fontsize=11)
    fig.tight_layout()
    return fig
