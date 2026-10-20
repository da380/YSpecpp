"""Pictures of seismograms: record sections, traces and comparisons."""

import numpy as np

#: A colour for each component.
COLOURS = {"Z": "#1f4e79", "N": "#7f2704", "E": "#2d6a4f",
           "R": "#7f2704", "T": "#2d6a4f"}


def record_section(run, component="Z", *, common=False, gain=1.0):
    """One component at every receiver against epicentral distance, each
    trace scaled to its peak (or all alike, with `common`), amplified by
    `gain` and clipped where it would reach the next."""
    import matplotlib.pyplot as plt

    receivers = sorted(run.receivers, key=lambda r: r.distance)
    traces = [r.trace(component) for r in receivers]
    distances = np.array([r.distance for r in receivers])
    spacing = np.median(np.diff(distances)) if len(receivers) > 1 else 1.0
    peak = max(np.abs(u).max() for _, u in traces) or 1.0
    fig, ax = plt.subplots(figsize=(8, 9))
    for distance, (t, u) in zip(distances, traces):
        scale = peak if common else (np.abs(u).max() or 1.0)
        u = np.clip(gain * u / scale, -1, 1)
        ax.plot(t / 60, distance + 0.9 * spacing * u, lw=0.6, color="k")
    t = traces[0][0]
    ax.set_xlim(t[0] / 60, t[-1] / 60)
    ax.set_ylim(distances[0] - spacing, distances[-1] + spacing)
    ax.set_xlabel("time after the origin (min)")
    ax.set_ylabel("epicentral distance (degrees)")
    scaling = "common scale" if common else "each trace scaled to its peak"
    if gain != 1:
        scaling += f", times {gain:g} and clipped"
    ax.set_title(f"{component}, {run.quantity}, source {run.source} "
                 f"({scaling})")
    return fig


def traces(run, numbers=None, components="ZRT"):
    """The components at some receivers (by default the first eight), a row
    for each."""
    import matplotlib.pyplot as plt

    numbers = numbers or list(range(1, min(len(run.receivers), 8) + 1))
    receivers = [run.receiver(k) for k in numbers]
    components = list(components)
    fig, axes = plt.subplots(len(receivers), len(components), sharex=True,
                             squeeze=False,
                             figsize=(3.4 * len(components) + 1.2,
                                      1.5 * len(receivers) + 1.2))
    for i, receiver in enumerate(receivers):
        for j, component in enumerate(components):
            t, u = receiver.trace(component)
            ax = axes[i][j]
            ax.plot(t / 60, u, lw=0.7, color=COLOURS[component])
            ax.set_xlim(t[0] / 60, t[-1] / 60)
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 3))
            ax.tick_params(labelsize=8)
            ax.yaxis.get_offset_text().set_fontsize(8)
            if i == 0:
                ax.set_title(component)
            if j == 0:
                ax.set_ylabel(f"{receiver.number}: {receiver.depth:g} km\n"
                              f"{receiver.distance:.1f}°, {receiver.units}",
                              fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel("time after the origin (min)")
    fig.suptitle(f"{run.quantity}, source {run.source}")
    fig.tight_layout()
    return fig


def compare(runs, labels=None, *, receiver=1, component="Z"):
    """One receiver in several runs, overlaid, and their differences from
    the last."""
    import matplotlib.pyplot as plt

    labels = labels or [run.path.stem for run in runs]
    if len(labels) != len(runs):
        raise ValueError("there must be a label for each run")
    receivers = [run.receiver(receiver) for run in runs]
    traces = [r.trace(component) for r in receivers]
    fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(9, 6))
    for (t, u), label in zip(traces, labels):
        top.plot(t / 60, u, lw=0.8, label=label)
    t, reference = traces[-1]
    for (t, u), label in zip(traces[:-1], labels[:-1]):
        bottom.plot(t / 60, u - reference, lw=0.8,
                    label=f"{label} $-$ {labels[-1]}")
    for ax in (top, bottom):
        ax.legend(fontsize=8, loc="upper right")
        ax.set_ylabel(receivers[-1].units)
    bottom.set_xlim(t[0] / 60, t[-1] / 60)
    bottom.set_xlabel("time after the origin (min)")
    last = receivers[-1]
    top.set_title(f"{component}, {runs[-1].quantity} at "
                  f"{last.distance:.1f}°, {last.depth:g} km deep")
    fig.tight_layout()
    return fig
