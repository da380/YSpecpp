"""Pictures of snapshots of the wavefield: sections through the earth and
maps around the source.

A grid of depths along a great circle (one azimuth, or two opposite) is a
section, drawn as planetmodel draws one, through a `SectionSlice` made from
the snapshot; a grid at one depth and many azimuths is a map, drawn in
azimuthal equidistant projection about the source with cartopy, so that its
distances and azimuths are those of the grid. Anything else is drawn as
lines against distance.
"""

import math

import numpy as np
from planetmodel.deck import MINEOS, read_deck
from planetmodel.plotting import colour_scale, section
from planetmodel.slices import SectionSlice, Shown, SliceLayer

from .runs import DIMENSIONS

#: The guides drawn at the boundaries between solid and fluid, over the
#: picture.
GUIDE = {"color": "0.25", "lw": 0.6, "alpha": 0.7, "zorder": 3}

#: The mark of the source.
STAR = {"marker": "*", "color": "gold", "markeredgecolor": "k",
        "markersize": 11, "zorder": 4}

#: The radius of the sphere on which the maps are drawn (m).
RADIUS = 6371e3


def _unit(distance, azimuth):
    """Unit vectors at angular distances and azimuths from a source at the
    pole, as YSpec places them: longitude pi minus the azimuth."""
    phi = np.pi - azimuth
    return np.stack([np.sin(distance) * np.cos(phi),
                     np.sin(distance) * np.sin(phi),
                     np.cos(distance) * np.ones_like(phi)], axis=-1)


def _section_slice(run, field, surface, shown):
    """The snapshot on a great circle as a slice: along the second azimuth
    from its far end to the source, then along the first, radii upward."""
    delta = np.radians(run.distances)
    pieces = [(delta, np.full_like(delta, math.radians(run.azimuths[0])),
               field[:, 0, :])]
    if len(run.azimuths) == 2:
        pieces.insert(0, (delta[::-1],
                          np.full_like(delta, math.radians(run.azimuths[1])),
                          field[:, 1, ::-1]))
    path = np.concatenate([_unit(d, a) for d, a, _ in pieces])
    along = (np.concatenate([np.pi - pieces[0][0], np.pi + delta])
             if len(pieces) == 2 else delta)
    along = along - along[0]
    values = np.concatenate([v for _, _, v in pieces], axis=1)[::-1]
    radius = (surface - run.depths)[::-1]
    grid = np.broadcast_to(radius[:, None], values.shape)
    points = grid[..., None] * path[None, :, :]
    travel = path[1] - path[0]
    travel = travel - np.dot(travel, path[0]) * path[0]
    plane = np.stack([path[0], travel / np.linalg.norm(travel)])
    layer = SliceLayer(index=0, radius=np.array(grid), values=values,
                       reference=points, physical=points)
    return SectionSlice(distance=along, path=path, plane=plane,
                        layers=(layer,), shown=shown)


def _phase_boundaries(model):
    """The radii (km) of the surface and of every boundary between solid and
    fluid in a deck."""
    deck = read_deck(model, MINEOS)
    shear = deck["vsv"] if "vsv" in deck else deck["vs"]
    r = deck.radius
    radii = [r[-1]]
    for i in range(1, len(r)):
        if r[i] == r[i - 1] and (shear[i] == 0) != (shear[i - 1] == 0):
            radii.append(r[i])
    return np.array(radii) / 1e3, r[-1] / 1e3


def _draw_section(ax, run, field, surface, boundaries, shown, limits):
    slice_ = _section_slice(run, field, surface, shown)
    section(ax, slice_, diverging=True, vmin=-limits, vmax=limits,
            interfaces=False)
    half = 0.5 * slice_.length
    arc = np.linspace(-half, half, 721)
    for r in boundaries:
        ax.plot(arc, np.full_like(arc, r), **GUIDE)
    at = 0.0 if len(run.azimuths) == 2 else -half
    ax.plot(at, surface - run.hypocentre["depth_km"], **STAR)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)


def _draw_map(ax, run, field, projection, limits, cmap, coastlines):
    azimuth = np.radians(run.azimuths)
    image = field[0]
    if not math.isclose(run.azimuths[-1] - run.azimuths[0], 360):
        azimuth = np.append(azimuth, azimuth[0] + 2 * np.pi)
        image = np.concatenate([image, image[:1]], axis=0)
    a, d = np.meshgrid(azimuth, np.radians(run.distances) * RADIUS,
                       indexing="ij")
    mesh = ax.pcolormesh(d * np.sin(a), d * np.cos(a), image,
                         transform=projection, shading="gouraud", cmap=cmap,
                         vmin=-limits, vmax=limits)
    if coastlines:
        ax.coastlines(lw=0.4, color="0.3")
    ax.plot(0, 0, transform=projection, **STAR)
    ax.set_global()
    return mesh


def snapshots(run, component="up", *, common=False, clip=99.5,
              coastlines=True):
    """The snapshots of a run at each of its times, a panel for each, each
    scaled to the `clip` percentile of its own values (or all alike, with
    `common`)."""
    import matplotlib.pyplot as plt

    shown = Shown(f"{component} {run.quantity}",
                  dimensions=DIMENSIONS[run.quantity])
    fields = [run.field(k, component) for k in range(len(run.times))]
    cmap, _, common_limit = colour_scale(fields, shown, diverging=True,
                                         clip=clip)
    columns = min(len(run.times), 3)
    rows = math.ceil(len(run.times) / columns)
    is_section = len(run.depths) > 1 and len(run.azimuths) <= 2
    is_map = len(run.depths) == 1 and len(run.azimuths) > 2
    fig = plt.figure(figsize=(4.2 * columns, 4.4 * rows + 0.8))
    if is_section:
        boundaries, surface = _phase_boundaries(run.model)
        kw = {"projection": "polar"}
    elif is_map:
        import cartopy.crs as ccrs
        globe = ccrs.Globe(ellipse=None, semimajor_axis=RADIUS,
                           semiminor_axis=RADIUS)
        projection = ccrs.AzimuthalEquidistant(
            central_longitude=run.hypocentre["longitude"],
            central_latitude=run.hypocentre["latitude"], globe=globe)
        kw = {"projection": projection}
    else:
        kw = {}
    axes = []
    for k, (time, field) in enumerate(zip(run.times, fields)):
        ax = fig.add_subplot(rows, columns, k + 1, **kw)
        axes.append(ax)
        limits = (common_limit if common else
                  colour_scale([field], shown, diverging=True, clip=clip)[2])
        if is_section:
            _draw_section(ax, run, field, surface, boundaries, shown, limits)
            ax.set_title(f"{time:g} s")
        elif is_map:
            mesh = _draw_map(ax, run, field, projection, limits, cmap,
                             coastlines)
            ax.set_title(f"{time:g} s, {run.depths[0]:g} km deep")
        else:
            for d, depth in enumerate(run.depths):
                for a, azimuth in enumerate(run.azimuths):
                    ax.plot(run.distances, field[d, a],
                            label=f"{depth:g} km, {azimuth:g}°")
            ax.set_xlabel("distance (degrees)")
            ax.set_ylabel(shown.label)
            ax.legend(fontsize=7)
            ax.set_title(f"{time:g} s")
    scaling = ("common scale" if common
               else "each panel scaled to its own range")
    fig.suptitle(f"{shown.label}, source {run.source} ({scaling})")
    if common and (is_section or is_map):
        mappable = plt.cm.ScalarMappable(
            norm=plt.Normalize(-common_limit, common_limit), cmap=cmap)
        fig.colorbar(mappable, ax=axes, shrink=0.6, label=shown.label)
    return fig
