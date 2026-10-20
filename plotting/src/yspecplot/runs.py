"""The output of a run, found from its input file."""

import json
import math
import pathlib
from dataclasses import dataclass

import numpy as np
from planetmodel.units import Dimensions

#: The dimensions of each quantity a run may record.
DIMENSIONS = {
    "displacement": Dimensions(length=1),
    "velocity": Dimensions(length=1, time=-1),
    "acceleration": Dimensions(length=1, time=-2),
}

#: The SI unit of each, for labels.
UNITS = {"displacement": "m", "velocity": "m/s", "acceleration": "m/s$^2$"}


def values(spec):
    """The values a number, an array or a range {min, max, count} gives, as
    the programs read them."""
    if isinstance(spec, (int, float)):
        return [float(spec)]
    if isinstance(spec, list):
        return [float(v) for v in spec]
    count = spec["count"]
    if count == 1:
        return [float(spec["min"])]
    return list(np.linspace(spec["min"], spec["max"], count))


def source_names(sources):
    """The names of the sources, numbered from 1 where not given."""
    return [s.get("name", str(i + 1)) for i, s in enumerate(sources)]


def _choose(names, source, path):
    name = names[0] if source is None else source
    if name not in names:
        raise ValueError(f"there is no source {name} in {path}")
    return name


@dataclass(frozen=True)
class Receiver:
    """A receiver as the table beside the seismograms gives it, as seen
    from one source."""

    number: int
    file: pathlib.Path
    depth: float          # km
    distance: float       # degrees
    azimuth: float        # degrees, at the source
    back_azimuth: float   # degrees, of the source at the receiver
    units: str

    def trace(self, component):
        """The times (s) and one component: Z (up), N, E, or R and T,
        radial (away from the source) and transverse (90 degrees
        anticlockwise from radial, seen from above, as in the snapshots).
        """
        t, z, n, e = np.loadtxt(self.file).T
        back = math.radians(self.back_azimuth)
        components = {
            "Z": lambda: z,
            "N": lambda: n,
            "E": lambda: e,
            "R": lambda: -n * math.cos(back) - e * math.sin(back),
            "T": lambda: -n * math.sin(back) + e * math.cos(back),
        }
        if component not in components:
            raise ValueError(f"unknown component {component}: use Z, N, E, R "
                             "or T")
        return t, components[component]()


class SeismogramRun:
    """A run of yspec: its input, and the receivers of one source."""

    def __init__(self, path, source=None):
        self.path = pathlib.Path(path)
        self.input = json.loads(self.path.read_text())
        output = self.input["output"]
        prefix = self.path.parent / output["prefix"]
        self.quantity = output.get("quantity", "displacement")
        self.source = _choose(source_names(self.input["sources"]), source,
                              path)
        table = prefix.with_name(prefix.name + ".receivers")
        rows = [line.split() for line in table.read_text().splitlines()
                if not line.startswith("#")]
        rows = [row for row in rows if row[1] == self.source]
        self.receivers = []
        for number, (row, instrumented) in enumerate(
                zip(rows, self._instruments()), start=1):
            file, _, _, _, depth, distance, azimuth, back = row
            self.receivers.append(Receiver(
                number=number, file=prefix.parent / file, depth=float(depth),
                distance=float(distance), azimuth=float(azimuth),
                back_azimuth=float(back),
                units=("instrument output" if instrumented
                       else UNITS[self.quantity])))

    def receiver(self, number):
        """A receiver, numbered from 1 as the files are."""
        if not 1 <= number <= len(self.receivers):
            raise ValueError(f"there is no receiver {number} in {self.path}")
        return self.receivers[number - 1]

    def _instruments(self):
        """Whether each receiver has an instrument, group by group, depth by
        depth, place by place, as the program numbers them."""
        groups = self.input["receivers"]
        if isinstance(groups, dict):
            groups = [groups]
        for group in groups:
            depths = group["depth_km"]
            depths = len(depths) if isinstance(depths, list) else 1
            places = (len(group["locations"]) if "locations" in group
                      else group["profile"]["count"])
            for _ in range(depths * places):
                yield "instrument" in group


class SnapshotRun:
    """A run of yspec-snapshot: its grid and the fields at each time, for
    one source, as arrays (depth, azimuth, distance)."""

    COLUMNS = {"up": 5, "away": 6, "transverse": 7}

    def __init__(self, path, source=None):
        self.path = pathlib.Path(path)
        self.input = json.loads(self.path.read_text())
        prefix = self.path.parent / self.input["output"]["prefix"]
        self.quantity = self.input["output"].get("quantity", "displacement")
        names = source_names(self.input["sources"])
        self.source = _choose(names, source, path)
        self.hypocentre = self.input["sources"][names.index(self.source)]
        self.model = self.path.parent / self.input["model"]
        self.times = values(self.input["snapshots"]["times_s"])
        self._data = []
        for k in range(1, len(self.times) + 1):
            suffix = f".{k}" if len(names) == 1 else f".{self.source}.{k}"
            self._data.append(np.loadtxt(prefix.with_name(prefix.name +
                                                          suffix)))
        first = self._data[0]
        self.depths = np.unique(first[:, 0])
        self.distances = np.unique(first[:, 1])
        self.azimuths = np.array(list(dict.fromkeys(first[:, 2])))

    def field(self, time_index, component):
        """One component at one time, (depth, azimuth, distance)."""
        shape = (len(self.depths), len(self.azimuths), len(self.distances))
        return self._data[time_index][:, self.COLUMNS[component]].reshape(shape)


class KernelRun:
    """A run of yspec-kernels: the measurement and its radial kernels."""

    def __init__(self, path):
        self.path = pathlib.Path(path)
        self.input = json.loads(self.path.read_text())
        output = self.input["output"]
        prefix = self.path.parent / output["prefix"]
        self.quantity = output.get("quantity", "displacement")
        self.model = self.path.parent / self.input["model"]
        self.hypocentre = self.input["sources"][0]
        measurement = self.input["measurement"]
        self.component = measurement["component"]
        self.window = measurement["traveltime"]["window_s"]
        kernels = self.input.get("kernels", {})
        self.relative = kernels.get("scaling", "absolute") == "relative"
        self.time, self.seismogram, self.taper, self.adjoint = np.loadtxt(
            prefix.with_name(prefix.name + ".measurement")).T
        file = prefix.with_name(prefix.name + ".kernels")
        header = [line for line in file.read_text().splitlines()
                  if line.startswith("#")]
        self.description = header[0].removeprefix("# yspec-kernels: ")
        self.at_source = [float(v) for v in next(
            line for line in header
            if line.startswith("# at_source")).split()[2:]]
        self.names = header[-1].split()[4:]
        table = np.loadtxt(file)
        self.radius, self.depth, self.weight = table[:, :3].T
        self.kernels = {name: table[:, 3 + j]
                        for j, name in enumerate(self.names)}
        self.at_source = dict(zip(self.names, self.at_source))
        self.source_radius = self.radius[-1] + self.depth[-1] - \
            self.hypocentre["depth_km"]
