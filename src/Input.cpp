#include <YSpec/Input.hpp>

#include "InputSections.hpp"

namespace YSpec {

using namespace InputSections;

SeismogramInput ReadSeismogramInput(const std::filesystem::path& path) {
  const auto json = planetmodel::ReadJsonFile(path);
  auto members = SharedMembers();
  members.insert(members.end(), {"output", "receivers"});
  Allow(json, members);
  const auto directory = path.parent_path();

  auto shared = ReadShared(json, directory);
  auto input = SeismogramInput{.model = std::move(shared.model),
                               .calculation = shared.calculation,
                               .sources = std::move(shared.sources),
                               .names = std::move(shared.names),
                               .receivers = {},
                               .prefix = {}};

  const auto& output = json.At("output");
  Allow(output, {"prefix", "quantity", "potential_and_tilt_corrections"});
  input.prefix = directory / output.At("prefix").AsString();
  auto& calculation = input.calculation;
  calculation.quantity = ReadQuantity(output, calculation.quantity);
  if (const auto* corrections = output.Find("potential_and_tilt_corrections")) {
    calculation.corrections = corrections->AsBool();
  }

  input.receivers = ReadReceivers(json.At("receivers"), directory);
  return input;
}

SnapshotInput ReadSnapshotInput(const std::filesystem::path& path) {
  const auto json = planetmodel::ReadJsonFile(path);
  auto members = SharedMembers();
  members.insert(members.end(), {"output", "snapshots"});
  Allow(json, members);
  const auto directory = path.parent_path();

  auto shared = ReadShared(json, directory);
  auto input = SnapshotInput{.model = std::move(shared.model),
                             .calculation = shared.calculation,
                             .sources = std::move(shared.sources),
                             .names = std::move(shared.names),
                             .grid = {},
                             .times = {},
                             .prefix = {}};

  const auto& output = json.At("output");
  Allow(output, {"prefix", "quantity"});
  input.prefix = directory / output.At("prefix").AsString();
  input.calculation.quantity = ReadQuantity(output, input.calculation.quantity);

  const auto& snapshots = json.At("snapshots");
  Allow(snapshots, {"times_s", "depths_km", "distances_deg", "azimuths_deg"});
  const auto& times = snapshots.At("times_s");
  input.times = ReadValues(times);
  for (auto time : input.times) {
    if (!(time >= 0 && time <= input.calculation.sampling.duration)) {
      Fail(times, "should lie between zero and the length of the time series");
    }
  }
  const auto& depths = snapshots.At("depths_km");
  input.grid.depths = ReadValues(depths, 1e3);
  for (auto depth : input.grid.depths) {
    if (!(depth >= 0)) Fail(depths, "should not be negative");
  }
  const auto& distances = snapshots.At("distances_deg");
  input.grid.distances = ReadValues(distances);
  for (auto distance : input.grid.distances) {
    if (!(distance >= 0 && distance <= 180)) {
      Fail(distances, "should lie between 0 and 180");
    }
  }
  input.grid.azimuths = ReadValues(snapshots.At("azimuths_deg"));
  return input;
}

}  // namespace YSpec
