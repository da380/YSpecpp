#include <YSpec/Model.hpp>
#include <YSpec/Run.hpp>
#include <format>
#include <fstream>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#include "Geometry.hpp"

namespace YSpec {

Seismograms ComputeSeismograms(const SeismogramInput& input,
                               const Progress& progress) {
  const auto model = Model::ReadDeck(input.model);
  return ComputeSeismograms(model, input.calculation, input.sources,
                            input.receivers, progress);
}

namespace {

// The file of the seismograms of a source at a receiver.
std::filesystem::path SeismogramPath(const SeismogramInput& input,
                                     std::size_t sources, std::size_t source,
                                     std::size_t receiver) {
  auto path = input.prefix;
  path += sources == 1
              ? std::format(".{}", receiver + 1)
              : std::format(".{}.{}", input.names.at(source), receiver + 1);
  return path;
}

// A table of the seismogram files, one line for each source and receiver:
// where the receiver is, and where it is seen from the source.
void WriteReceivers(const SeismogramInput& input, std::size_t sources) {
  auto path = input.prefix;
  path += ".receivers";
  auto file = std::ofstream(path);
  if (!file) throw std::runtime_error("cannot write " + path.string());
  file << "# file source latitude longitude depth_km distance_deg "
          "azimuth_deg back_azimuth_deg\n";
  for (auto source = std::size_t{0}; source < sources; ++source) {
    for (auto receiver = std::size_t{0}; receiver < input.receivers.size();
         ++receiver) {
      const auto& at = input.receivers[receiver];
      const auto view = Geometry(at.location, input.sources[source].location);
      const auto degree = std::numbers::pi / 180;
      file << std::format(
          "{} {} {:10.4f}{:10.4f}{:10.3f}{:10.4f}{:10.4f}{:10.4f}\n",
          SeismogramPath(input, sources, source, receiver).filename().string(),
          input.names.at(source), at.location.latitude, at.location.longitude,
          at.depth / 1e3, view.distance / degree,
          (std::numbers::pi - view.phi) / degree,
          (2 * std::numbers::pi - view.rotation) / degree);
    }
  }
  if (!file) throw std::runtime_error("cannot write " + path.string());
}

}  // namespace

void WriteSeismograms(const Seismograms& seismograms,
                      const SeismogramInput& input) {
  WriteReceivers(input, seismograms.Sources());
  for (auto source = std::size_t{0}; source < seismograms.Sources(); ++source) {
    for (auto receiver = std::size_t{0}; receiver < seismograms.Receivers();
         ++receiver) {
      const auto path =
          SeismogramPath(input, seismograms.Sources(), source, receiver);
      auto file = std::ofstream(path);
      if (!file) throw std::runtime_error("cannot write " + path.string());
      const auto vertical = seismograms.Trace(source, receiver, Vertical);
      const auto north = seismograms.Trace(source, receiver, North);
      const auto east = seismograms.Trace(source, receiver, East);
      for (auto i = std::size_t{0}; i < seismograms.Samples(); ++i) {
        file << std::format("{:15.6E}{:15.6E}{:15.6E}{:15.6E}\n",
                            static_cast<double>(i) * seismograms.TimeStep(),
                            vertical[i], north[i], east[i]);
      }
      if (!file) throw std::runtime_error("cannot write " + path.string());
    }
  }
}

Seismograms RunSeismograms(const std::filesystem::path& path,
                           const Progress& progress) {
  const auto input = ReadSeismogramInput(path);
  auto seismograms = ComputeSeismograms(input, progress);
  WriteSeismograms(seismograms, input);
  return seismograms;
}

Snapshots ComputeSnapshots(const SnapshotInput& input,
                           const Progress& progress) {
  const auto model = Model::ReadDeck(input.model);
  return ComputeSnapshots(model, input.calculation, input.sources, input.grid,
                          input.times, progress);
}

void WriteSnapshots(const Snapshots& snapshots, const SnapshotInput& input) {
  const auto& grid = snapshots.Grid();
  const auto quantity = [&] {
    switch (input.calculation.quantity) {
      case Quantity::Displacement:
        return "displacement (m)";
      case Quantity::Velocity:
        return "velocity (m/s)";
      case Quantity::Acceleration:
        return "acceleration (m/s^2)";
    }
    return "";
  }();
  for (auto source = std::size_t{0}; source < snapshots.Sources(); ++source) {
    // Where the points are, around this source.
    auto places = std::vector<Location>{};
    for (auto azimuth : grid.azimuths) {
      for (auto distance : grid.distances) {
        places.push_back(GreatCircle(input.sources[source].location, azimuth,
                                     distance, distance, 1)
                             .front());
      }
    }
    for (auto time = std::size_t{0}; time < snapshots.Times().size(); ++time) {
      auto path = input.prefix;
      path += snapshots.Sources() == 1
                  ? std::format(".{}", time + 1)
                  : std::format(".{}.{}", input.names.at(source), time + 1);
      auto file = std::ofstream(path);
      if (!file) throw std::runtime_error("cannot write " + path.string());
      file << std::format("# source {}, {} s after the origin, {}\n",
                          input.names.at(source), snapshots.Times()[time],
                          quantity);
      file << "# depth_km distance_deg azimuth_deg latitude longitude up away "
              "transverse\n";
      for (auto depth = std::size_t{0}; depth < grid.depths.size(); ++depth) {
        for (auto a = std::size_t{0}; a < grid.azimuths.size(); ++a) {
          for (auto i = std::size_t{0}; i < grid.distances.size(); ++i) {
            const auto& place = places[a * grid.distances.size() + i];
            file << std::format(
                "{:10.3f}{:10.3f}{:10.3f}{:10.4f}{:10.4f}{:15.6E}{:15.6E}"
                "{:15.6E}\n",
                grid.depths[depth] / 1e3, grid.distances[i], grid.azimuths[a],
                place.latitude, place.longitude,
                snapshots(source, time, depth, i, a, Up),
                snapshots(source, time, depth, i, a, Away),
                snapshots(source, time, depth, i, a, Transverse));
          }
        }
      }
      if (!file) throw std::runtime_error("cannot write " + path.string());
    }
  }
}

Snapshots RunSnapshots(const std::filesystem::path& path,
                       const Progress& progress) {
  const auto input = ReadSnapshotInput(path);
  auto snapshots = ComputeSnapshots(input, progress);
  WriteSnapshots(snapshots, input);
  return snapshots;
}

KernelResult ComputeKernels(const KernelInput& input,
                            const Progress& progress) {
  const auto model = Model::ReadDeck(input.model);
  const auto& calculation = input.calculation;
  auto seismogram =
      ComputeSeismograms(model, calculation, std::span{&input.source, 1},
                         std::span{&input.receiver, 1});
  auto measurement = Traveltime(model, calculation, input.source,
                                input.receiver, input.direction, input.window);
  auto kernels = ComputeRadialKernels(model, calculation, input.source,
                                      measurement, input.mesh, progress);
  auto named = Reparameterise(kernels, model, input.parameters, input.scaling);
  return {.seismogram = std::move(seismogram),
          .measurement = std::move(measurement),
          .kernels = std::move(kernels),
          .named = std::move(named)};
}

void WriteKernels(const KernelResult& result, const KernelInput& input) {
  const auto quantity = [&] {
    switch (input.calculation.quantity) {
      case Quantity::Displacement:
        return "displacement";
      case Quantity::Velocity:
        return "velocity";
      case Quantity::Acceleration:
        return "acceleration";
    }
    return "";
  }();
  const auto& window = input.window;
  const auto measured = std::format(
      "the cross-correlation traveltime (s) of {} {} at the receiver, in a "
      "window from {} s to {} s",
      input.component, quantity, window.start, window.end);

  auto open = [&](std::string_view suffix) {
    auto path = input.prefix;
    path += suffix;
    auto file = std::ofstream(path);
    if (!file) throw std::runtime_error("cannot write " + path.string());
    return std::pair{std::move(file), path};
  };

  {
    auto [file, path] = open(".measurement");
    file << "# yspec-kernels: " << measured << "\n";
    file << "# time_s seismogram window adjoint_source\n";
    const auto& d = input.direction;
    const auto& u = result.seismogram;
    const auto& j = result.measurement.adjointSource;
    for (auto n = std::size_t{0}; n < u.Samples(); ++n) {
      const auto t = static_cast<double>(n) * u.TimeStep();
      file << std::format("{:12.3f}{:15.6E}{:10.4f}{:15.6E}\n", t,
                          d.up * u(0, 0, Vertical, n) +
                              d.north * u(0, 0, North, n) +
                              d.east * u(0, 0, East, n),
                          window(t), j[n]);
    }
    if (!file) throw std::runtime_error("cannot write " + path.string());
  }

  auto [file, path] = open(".kernels");
  const auto& named = result.named;
  const auto relative = input.scaling == Scaling::Relative;
  file << "# yspec-kernels: radial kernels of " << measured << "\n";
  file << std::format(
      "# dJ = integral of sum over X of K_X d{}X r^2 dr, plus sum over X of "
      "P_X d{}X at the source; K in s m^-3{}, P in s{}\n",
      relative ? "ln " : "", relative ? "ln " : "",
      relative ? "" : " per unit of X", relative ? "" : " per unit of X");
  file << std::format("# the point parts P_X at the source, at {:.3f} km:\n",
                      result.kernels.SourceRadius() / 1e3);
  file << "# at_source";
  for (auto value : named.atSource) file << std::format(" {:.6E}", value);
  file << "\n# radius_km depth_km weight_m3";
  for (const auto& name : named.names) file << " " << name;
  file << "\n";
  const auto& kernels = result.kernels;
  const auto surface = kernels.Radii().empty() ? 0.0 : [&] {
    const auto model = Model::ReadDeck(input.model);
    return model.SurfaceRadius() * model.GetScales().Length();
  }();
  for (auto i = std::size_t{0}; i < kernels.Size(); ++i) {
    const auto r = kernels.Radii()[i];
    file << std::format("{:15.6f}{:15.6f}{:15.6E}", r / 1e3,
                        (surface - r) / 1e3, kernels.Weights()[i]);
    for (const auto& values : named.values) {
      file << std::format("{:15.6E}", values[i]);
    }
    file << "\n";
  }
  if (!file) throw std::runtime_error("cannot write " + path.string());
}

KernelResult RunKernels(const std::filesystem::path& path,
                        const Progress& progress) {
  const auto input = ReadKernelInput(path);
  auto result = ComputeKernels(input, progress);
  WriteKernels(result, input);
  return result;
}

}  // namespace YSpec
