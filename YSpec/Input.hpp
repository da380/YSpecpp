#pragma once

/**
 * @file Input.hpp
 * @brief Reads a calculation from a JSON file.
 *
 * @details The file is strict JSON (RFC 8259). Its members, with the units in
 * their names, are:
 *
 * @code{.json}
 * {
 *   "model": "prem.200",
 *   "attenuation": true,
 *   "gravity": "full",
 *   "output": {
 *     "prefix": "yspec.out",
 *     "quantity": "velocity",
 *     "potential_and_tilt_corrections": false
 *   },
 *   "degrees": { "min": 0, "max": 600 },
 *   "frequencies_mHz": { "min": 0.2, "max": 40 },
 *   "time_series": { "length_min": 60, "time_step_s": 1.0 },
 *   "filter_mHz": [1.0, 5.0, 38, 40],
 *   "source": {
 *     "depth_km": 29.5, "latitude": 36.87, "longitude": 35.58,
 *     "moment_tensor_Nm": { "rr": 4.22e17, "rt": 3.08e17, "rp": 6.99e17,
 *                           "tt": -2.976e18, "tp": 6.54e17, "pp": 2.55e18 }
 *   },
 *   "receivers": { "depth_km": 3, "locations": [[0, 4], [0, 6]] },
 *   "solver": { "tolerance": 1e-6, "start_decay": 11 }
 * }
 * @endcode
 *
 * "gravity" is one of "none", "cowling" and "full"; "quantity" one of
 * "displacement", "velocity" and "acceleration". The moment tensor is in
 * spherical coordinates at the source (r up, theta south, phi east). Receiver
 * locations are [latitude, longitude] in geographic degrees. The filter is a
 * raised-cosine taper rising between its first two corners and falling
 * between its last two; it may be left out. "attenuation", "gravity", the
 * corrections, "degrees.min", the filter and "solver" have defaults (true,
 * "full", false, 0, none, and the solver's defaults). Relative paths are
 * relative to the directory of the input file. A member not listed here is an
 * error, so that a misspelt one is not silently ignored.
 */

#include <array>
#include <filesystem>
#include <format>
#include <initializer_list>
#include <planetmodel/json.hpp>
#include <stdexcept>
#include <string>
#include <string_view>

#include "Synthetics.hpp"

namespace YSpec {

/** @brief A calculation read from a file. */
struct Input {
  Settings settings;             ///< The calculation.
  std::filesystem::path model;   ///< The model deck.
  std::filesystem::path prefix;  ///< The prefix of the output files.
};

namespace InputDetails {

using planetmodel::Json;

[[noreturn]] inline void Fail(const Json& json, std::string_view message) {
  throw planetmodel::Error(json.Source(), json.Path(), message);
}

// Refuses any member of an object not among those named.
inline void Allow(const Json& json,
                  std::initializer_list<std::string_view> names) {
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

inline double Number(const Json& object, std::string_view key) {
  return object.At(key).AsNumber();
}

inline double Number(const Json& object, std::string_view key,
                     double fallback) {
  const auto* found = object.Find(key);
  return found ? found->AsNumber() : fallback;
}

inline double Latitude(const Json& json) {
  const auto latitude = json.AsNumber();
  if (!(latitude >= -90 && latitude <= 90)) {
    Fail(json, "should be a latitude, between -90 and 90 degrees");
  }
  return latitude;
}

inline Location ReadLocation(const Json& json) {
  if (json.Size() != 2) Fail(json, "should be [latitude, longitude]");
  return {Latitude(json.At(0)), json.At(1).AsNumber()};
}

inline double Depth(const Json& object) {
  const auto& json = object.At("depth_km");
  const auto depth = json.AsNumber();
  if (!(depth >= 0)) Fail(json, "should not be negative");
  return depth;
}

}  // namespace InputDetails

/**
 * @brief Reads a calculation from a JSON file.
 * @throws planetmodel::Error naming the file, the member and the fault.
 */
inline Input ReadInput(const std::filesystem::path& path) {
  using namespace InputDetails;
  const auto json = planetmodel::ReadJsonFile(path);
  Allow(json, {"model", "attenuation", "gravity", "output", "degrees",
               "frequencies_mHz", "time_series", "filter_mHz", "source",
               "receivers", "solver"});
  const auto directory = path.parent_path();
  auto input = Input{};
  auto& s = input.settings;

  input.model = directory / json.At("model").AsString();
  if (const auto* a = json.Find("attenuation")) s.attenuation = a->AsBool();
  if (const auto* g = json.Find("gravity")) {
    s.gravity = Choice<Gravity>(*g, {{"none", Gravity::None},
                                     {"cowling", Gravity::Cowling},
                                     {"full", Gravity::Full}});
  }

  const auto& output = json.At("output");
  Allow(output, {"prefix", "quantity", "potential_and_tilt_corrections"});
  input.prefix = directory / output.At("prefix").AsString();
  if (const auto* q = output.Find("quantity")) {
    s.quantity =
        Choice<Quantity>(*q, {{"displacement", Quantity::Displacement},
                              {"velocity", Quantity::Velocity},
                              {"acceleration", Quantity::Acceleration}});
  }
  if (const auto* c = output.Find("potential_and_tilt_corrections")) {
    s.corrections = c->AsBool();
  }

  const auto& degrees = json.At("degrees");
  Allow(degrees, {"min", "max"});
  if (const auto* l = degrees.Find("min")) s.lMin = l->AsInt();
  s.lMax = degrees.At("max").AsInt();
  if (s.lMin < 0 || s.lMax < s.lMin) {
    Fail(degrees, "should satisfy 0 <= min <= max");
  }

  const auto& frequencies = json.At("frequencies_mHz");
  Allow(frequencies, {"min", "max"});
  s.fMin = Number(frequencies, "min");
  s.fMax = Number(frequencies, "max");
  if (!(s.fMin > 0 && s.fMax > s.fMin)) {
    Fail(frequencies, "should satisfy 0 < min < max");
  }

  const auto& series = json.At("time_series");
  Allow(series, {"length_min", "time_step_s"});
  s.length = Number(series, "length_min");
  s.timeStep = Number(series, "time_step_s");
  if (!(s.length > 0 && s.timeStep > 0)) {
    Fail(series, "should have a positive length and time step");
  }

  if (const auto* filter = json.Find("filter_mHz")) {
    if (filter->Size() != 4) Fail(*filter, "should have four corners");
    for (auto i : {0, 1, 2, 3}) {
      s.filter[i] = filter->At(static_cast<std::size_t>(i)).AsNumber();
    }
    const auto none = s.filter == std::array<double, 4>{};
    if (!none && !(s.filter[0] >= 0 && s.filter[0] <= s.filter[1] &&
                   s.filter[1] < s.filter[2] && s.filter[2] <= s.filter[3])) {
      Fail(*filter,
           "should be in increasing order, with a pass band between the "
           "second and third corners (or all zero for no filter)");
    }
  }

  const auto& source = json.At("source");
  Allow(source, {"depth_km", "latitude", "longitude", "moment_tensor_Nm"});
  s.sourceDepth = Depth(source);
  s.source = {Latitude(source.At("latitude")), Number(source, "longitude")};
  const auto& moment = source.At("moment_tensor_Nm");
  Allow(moment, {"rr", "rt", "rp", "tt", "tp", "pp"});
  s.moment = {Number(moment, "rr"), Number(moment, "rt"), Number(moment, "rp"),
              Number(moment, "tt"), Number(moment, "tp"), Number(moment, "pp")};

  const auto& receivers = json.At("receivers");
  Allow(receivers, {"depth_km", "locations"});
  s.receiverDepth = Depth(receivers);
  for (const auto& location : receivers.At("locations").AsArray()) {
    s.receivers.push_back(ReadLocation(location));
  }
  if (s.receivers.empty()) Fail(receivers.At("locations"), "is empty");

  if (const auto* solver = json.Find("solver")) {
    Allow(*solver, {"tolerance", "start_decay"});
    s.solver.integration.tolerance =
        Number(*solver, "tolerance", s.solver.integration.tolerance);
    s.solver.startDecay = Number(*solver, "start_decay", s.solver.startDecay);
    if (!(s.solver.integration.tolerance > 0 &&
          s.solver.integration.tolerance < 1)) {
      Fail(solver->At("tolerance"), "should lie between 0 and 1");
    }
    if (!(s.solver.startDecay > 0)) {
      Fail(solver->At("start_decay"), "should be positive");
    }
  }
  return input;
}

}  // namespace YSpec
