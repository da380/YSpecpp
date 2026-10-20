#include <YSpec/Model.hpp>
#include <YSpec/Run.hpp>
#include <format>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace YSpec {

Seismograms ComputeSeismograms(const SeismogramInput& input,
                               const Progress& progress) {
  const auto model = Model::ReadDeck(input.model);
  return ComputeSeismograms(model, input.calculation, input.sources,
                            input.receivers, progress);
}

void WriteSeismograms(const Seismograms& seismograms,
                      const SeismogramInput& input) {
  for (auto source = std::size_t{0}; source < seismograms.Sources(); ++source) {
    for (auto receiver = std::size_t{0}; receiver < seismograms.Receivers();
         ++receiver) {
      auto path = input.prefix;
      path += seismograms.Sources() == 1
                  ? std::format(".{}", receiver + 1)
                  : std::format(".{}.{}", input.names.at(source), receiver + 1);
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

}  // namespace YSpec
