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

KernelInput ReadKernelInput(const std::filesystem::path& path) {
  const auto json = planetmodel::ReadJsonFile(path);
  auto members = SharedMembers();
  members.insert(members.end(), {"output", "measurement", "kernels"});
  Allow(json, members);
  const auto directory = path.parent_path();

  auto shared = ReadShared(json, directory);
  if (shared.sources.size() != 1) {
    Fail(json.At("sources"), "should hold one source for kernels");
  }
  auto input = KernelInput{.model = std::move(shared.model),
                           .calculation = shared.calculation,
                           .source = shared.sources.front(),
                           .name = shared.names.front(),
                           .receiver = {},
                           .component = {},
                           .direction = {},
                           .window = {},
                           .parameters = ParameterSet::Velocities,
                           .scaling = Scaling::Absolute,
                           .mesh = {},
                           .prefix = {}};

  const auto& output = json.At("output");
  Allow(output, {"prefix", "quantity"});
  input.prefix = directory / output.At("prefix").AsString();
  input.calculation.quantity = ReadQuantity(output, input.calculation.quantity);

  const auto& measurement = json.At("measurement");
  Allow(measurement, {"receiver", "component", "traveltime"});
  input.receiver = ReadReceiver(measurement.At("receiver"),
                                input.source.location, directory);
  const auto& component = measurement.At("component");
  input.component = component.AsString();
  const auto& from = input.source.location;
  const auto& at = input.receiver.location;
  input.direction =
      Choice<Direction>(component, {{"Z", Along(Vertical)},
                                    {"N", Along(North)},
                                    {"E", Along(East)},
                                    {"R", RadialDirection(from, at)},
                                    {"T", TransverseDirection(from, at)}});
  const auto& traveltime = measurement.At("traveltime");
  Allow(traveltime, {"window_s", "taper_s"});
  const auto& window = traveltime.At("window_s");
  if (window.Size() != 2) Fail(window, "should be [start, end]");
  input.window = {.start = window.At(0).AsNumber(),
                  .end = window.At(1).AsNumber(),
                  .taper = Number(traveltime, "taper_s", 0)};
  const auto& w = input.window;
  if (!(w.start >= 0 && w.start < w.end &&
        w.end <= input.calculation.sampling.duration)) {
    Fail(window,
         "should satisfy 0 <= start < end <= the length of the time series");
  }
  if (!(w.taper >= 0 && 2 * w.taper <= w.end - w.start)) {
    Fail(traveltime.At("taper_s"),
         "should be at least zero and at most half the window");
  }

  if (const auto* kernels = json.Find("kernels")) {
    Allow(*kernels, {"parameters", "scaling", "mesh"});
    if (const auto* parameters = kernels->Find("parameters")) {
      input.parameters = Choice<ParameterSet>(
          *parameters, {{"love", ParameterSet::Love},
                        {"moduli", ParameterSet::Moduli},
                        {"velocities", ParameterSet::Velocities},
                        {"anisotropic", ParameterSet::Anisotropic}});
    }
    if (const auto* scaling = kernels->Find("scaling")) {
      input.scaling = Choice<Scaling>(
          *scaling,
          {{"absolute", Scaling::Absolute}, {"relative", Scaling::Relative}});
    }
    if (const auto* mesh = kernels->Find("mesh")) {
      Allow(*mesh, {"order", "max_element_km"});
      if (const auto* order = mesh->Find("order")) {
        input.mesh.order = order->AsInt();
        if (input.mesh.order < 1) Fail(*order, "should be at least 1");
      }
      input.mesh.maxElement =
          Number(*mesh, "max_element_km", input.mesh.maxElement / 1e3) * 1e3;
      if (!(input.mesh.maxElement > 0)) {
        Fail(mesh->At("max_element_km"), "should be positive");
      }
    }
  }
  return input;
}

}  // namespace YSpec
