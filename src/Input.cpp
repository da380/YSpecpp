#include <YSpec/Input.hpp>
#include <algorithm>
#include <format>
#include <initializer_list>
#include <planetmodel/json.hpp>
#include <set>
#include <string_view>
#include <utility>

namespace YSpec {

namespace {

using planetmodel::Json;

[[noreturn]] void Fail(const Json& json, std::string_view message) {
  throw planetmodel::Error(json.Source(), json.Path(), message);
}

// Refuses any member of an object not among those named.
void Allow(const Json& json, std::initializer_list<std::string_view> names) {
  for (const auto& [name, value] : json.AsObject()) {
    if (std::ranges::find(names, name) == names.end()) {
      Fail(value, "is not a recognised member");
    }
  }
}

template <class T>
T Choice(const Json& json,
         std::initializer_list<std::pair<std::string_view, T>> choices) {
  const auto& word = json.AsString();
  for (const auto& [name, value] : choices) {
    if (word == name) return value;
  }
  auto names = std::string{};
  for (const auto& [name, value] : choices) {
    names += (names.empty() ? "\"" : ", \"") + std::string(name) + "\"";
  }
  Fail(json, std::format("is \"{}\", which is none of {}", word, names));
}

double Number(const Json& object, std::string_view key) {
  return object.At(key).AsNumber();
}

double Number(const Json& object, std::string_view key, double fallback) {
  const auto* found = object.Find(key);
  return found ? found->AsNumber() : fallback;
}

double Latitude(const Json& json) {
  const auto latitude = json.AsNumber();
  if (!(latitude >= -90 && latitude <= 90)) {
    Fail(json, "should be a latitude, between -90 and 90 degrees");
  }
  return latitude;
}

Location ReadLocation(const Json& json) {
  if (json.Size() != 2) Fail(json, "should be [latitude, longitude]");
  return {Latitude(json.At(0)), json.At(1).AsNumber()};
}

// A depth in km, returned in m.
double Depth(const Json& object) {
  const auto& json = object.At("depth_km");
  const auto depth = json.AsNumber();
  if (!(depth >= 0)) Fail(json, "should not be negative");
  return depth * 1e3;
}

Source ReadSource(const Json& json) {
  Allow(json,
        {"name", "depth_km", "latitude", "longitude", "moment_tensor_Nm"});
  const auto& moment = json.At("moment_tensor_Nm");
  Allow(moment, {"rr", "rt", "rp", "tt", "tp", "pp"});
  return {
      .location = {Latitude(json.At("latitude")), Number(json, "longitude")},
      .depth = Depth(json),
      .moment = {Number(moment, "rr"), Number(moment, "rt"),
                 Number(moment, "rp"), Number(moment, "tt"),
                 Number(moment, "tp"), Number(moment, "pp")}};
}

}  // namespace

Input ReadInput(const std::filesystem::path& path) {
  const auto json = planetmodel::ReadJsonFile(path);
  Allow(json, {"model", "attenuation", "gravity", "output", "degrees",
               "frequencies_mHz", "time_series", "filter_mHz", "sources",
               "receivers", "solver"});
  const auto directory = path.parent_path();
  auto input = Input{};
  auto& calculation = input.calculation;

  input.model = directory / json.At("model").AsString();
  if (const auto* attenuation = json.Find("attenuation")) {
    calculation.attenuation = attenuation->AsBool();
  }
  if (const auto* gravity = json.Find("gravity")) {
    calculation.gravity =
        Choice<Gravity>(*gravity, {{"none", Gravity::None},
                                   {"cowling", Gravity::Cowling},
                                   {"full", Gravity::Full}});
  }

  const auto& output = json.At("output");
  Allow(output, {"prefix", "quantity", "potential_and_tilt_corrections"});
  input.prefix = directory / output.At("prefix").AsString();
  if (const auto* quantity = output.Find("quantity")) {
    calculation.quantity =
        Choice<Quantity>(*quantity, {{"displacement", Quantity::Displacement},
                                     {"velocity", Quantity::Velocity},
                                     {"acceleration", Quantity::Acceleration}});
  }
  if (const auto* corrections = output.Find("potential_and_tilt_corrections")) {
    calculation.corrections = corrections->AsBool();
  }

  const auto& degrees = json.At("degrees");
  Allow(degrees, {"min", "max"});
  if (const auto* min = degrees.Find("min")) {
    calculation.degrees.min = min->AsInt();
  }
  calculation.degrees.max = degrees.At("max").AsInt();
  if (calculation.degrees.min < 0 ||
      calculation.degrees.max < calculation.degrees.min) {
    Fail(degrees, "should satisfy 0 <= min <= max");
  }

  const auto& band = json.At("frequencies_mHz");
  Allow(band, {"min", "max"});
  calculation.band = {Number(band, "min") * 1e-3, Number(band, "max") * 1e-3};
  if (!(calculation.band.min > 0 &&
        calculation.band.max > calculation.band.min)) {
    Fail(band, "should satisfy 0 < min < max");
  }

  const auto& series = json.At("time_series");
  Allow(series, {"length_min", "time_step_s"});
  calculation.sampling = {Number(series, "length_min") * 60,
                          Number(series, "time_step_s")};
  if (!(calculation.sampling.duration > 0 &&
        calculation.sampling.timeStep > 0)) {
    Fail(series, "should have a positive length and time step");
  }

  if (const auto* filter = json.Find("filter_mHz")) {
    if (filter->Size() != 4) Fail(*filter, "should have four corners");
    auto corner = std::array<double, 4>{};
    for (auto i : {0, 1, 2, 3}) {
      corner[i] = filter->At(static_cast<std::size_t>(i)).AsNumber() * 1e-3;
    }
    if (!(corner[0] >= 0 && corner[0] <= corner[1] && corner[1] < corner[2] &&
          corner[2] <= corner[3])) {
      Fail(*filter,
           "should be in increasing order, with a pass band between the "
           "second and third corners");
    }
    calculation.taper = Taper{corner[0], corner[1], corner[2], corner[3]};
  }

  const auto& sources = json.At("sources");
  for (const auto& source : sources.AsArray()) {
    input.sources.push_back(ReadSource(source));
    const auto* name = source.Find("name");
    input.names.push_back(name ? name->AsString()
                               : std::to_string(input.sources.size()));
    if (input.names.back().empty() ||
        input.names.back().find('/') != std::string::npos) {
      Fail(name ? *name : source, "should be a name without a '/'");
    }
  }
  if (input.sources.empty()) Fail(sources, "is empty");
  if (std::set(input.names.begin(), input.names.end()).size() !=
      input.names.size()) {
    Fail(sources, "should have different names");
  }

  const auto& receivers = json.At("receivers");
  Allow(receivers, {"depth_km", "locations"});
  const auto depth = Depth(receivers);
  for (const auto& location : receivers.At("locations").AsArray()) {
    input.receivers.push_back({ReadLocation(location), depth});
  }
  if (input.receivers.empty()) Fail(receivers.At("locations"), "is empty");

  if (const auto* solver = json.Find("solver")) {
    Allow(*solver, {"tolerance", "start_decay"});
    auto& options = calculation.solver;
    options.tolerance = Number(*solver, "tolerance", options.tolerance);
    options.startDecay = Number(*solver, "start_decay", options.startDecay);
    if (!(options.tolerance > 0 && options.tolerance < 1)) {
      Fail(solver->At("tolerance"), "should lie between 0 and 1");
    }
    if (!(options.startDecay > 0)) {
      Fail(solver->At("start_decay"), "should be positive");
    }
  }
  return input;
}

}  // namespace YSpec
