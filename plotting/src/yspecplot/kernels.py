"""Pictures of the radial kernels of a measurement."""

import numpy as np
from planetmodel.deck import MINEOS, read_deck
from planetmodel.plotting import BOUNDARY_STYLE, profile
from planetmodel.slices import ProfileSlice, Shown, SliceLayer

from .seismograms import COLOURS
from .snapshots import GUIDE, STAR, _phase_boundaries


def _profile(radius, values, boundaries, name):
    """The kernel at its radii as a profile slice, a layer for each span
    between the deck's boundaries."""
    layers = []
    which = np.searchsorted(boundaries, radius)
    for index in np.unique(which):
        inside = which == index
        r = radius[inside]
        points = np.stack([np.zeros_like(r), np.zeros_like(r), r], axis=-1)
        layers.append(SliceLayer(index=int(index), radius=r,
                                 values=values[inside], reference=points,
                                 physical=points))
    return ProfileSlice(theta=None, phi=None, layers=tuple(layers),
                        shown=Shown(name))


def radial_kernels(run, names=None, *, deepest=None):
    """The measured seismogram and its window above, and the kernel of each
    parameter (by default all of them) against radius below, each with its
    point part at the source marked, down to `deepest` km if given."""
    import matplotlib.pyplot as plt

    names = names or run.names
    deck = read_deck(run.model, MINEOS)
    boundaries = deck.boundaries / 1e3
    phase, surface = _phase_boundaries(run.model)
    fig = plt.figure(figsize=(2.6 * len(names) + 1.5, 9))
    grid = fig.add_gridspec(2, len(names), height_ratios=[1, 3.2])
    top = fig.add_subplot(grid[0, :])
    t = run.time / 60
    peak = np.abs(run.seismogram).max() or 1.0
    top.plot(t, run.seismogram, lw=0.7,
             color=COLOURS.get(run.component, "k"))
    top.fill_between(t, -peak * run.taper, peak * run.taper, color="0.85",
                     zorder=0, label="window")
    top.set_xlim(t[0], t[-1])
    top.set_xlabel("time after the origin (min)")
    top.set_ylabel(f"{run.component}, {run.quantity}")
    top.legend(fontsize=8, loc="upper right")
    top.set_title(run.description, fontsize=9)
    unit = "s m$^{-3}$" if run.relative else "s m$^{-3}$ per unit"
    lowest = 0.0 if deepest is None else surface - deepest
    axes = []
    for j, name in enumerate(names):
        ax = fig.add_subplot(grid[1, j], sharey=axes[0] if axes else None)
        axes.append(ax)
        values = run.kernels[name]
        keep = run.radius >= lowest
        profile(ax, _profile(run.radius[keep], values[keep], boundaries,
                             name), joiners=False, boundaries=False,
                color="k", lw=0.9)
        for b in boundaries:
            if b >= lowest:
                ax.axhline(b, **BOUNDARY_STYLE)
        for b in phase:
            if b >= lowest:
                ax.axhline(b, **{k: v for k, v in GUIDE.items()
                                 if k != "zorder"})
        ax.axvline(0, color="0.6", lw=0.6)
        point = run.at_source[name]
        ax.plot(0, run.source_radius, **STAR)
        if point != 0:
            ax.annotate(f"P = {point:.2g}", (0, run.source_radius),
                        xytext=(4, -12), textcoords="offset points",
                        fontsize=7)
        largest = np.abs(values[keep]).max()
        ax.set_xlim(-1.1 * largest or -1, 1.1 * largest or 1)
        ax.ticklabel_format(axis="x", style="sci", scilimits=(-2, 3))
        ax.tick_params(labelsize=8)
        ax.xaxis.get_offset_text().set_fontsize(8)
        ax.set_title(("ln " if run.relative else "") + name)
        ax.set_xlabel(unit, fontsize=8)
        if j > 0:
            ax.tick_params(labelleft=False)
    axes[0].set_ylabel("radius (km)")
    axes[0].set_ylim(lowest, surface)
    fig.tight_layout()
    return fig
