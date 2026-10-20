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
 *   "sources": [
 *     { "name": "event", "depth_km": 29.5, "latitude": 36.87,
 *       "longitude": 35.58,
 *       "moment_tensor_Nm": { "rr": 4.22e17, "rt": 3.08e17, "rp": 6.99e17,
 *                             "tt": -2.976e18, "tp": 6.54e17, "pp": 2.55e18 }
 *     }
 *   ],
 *   "receivers": { "depth_km": 3, "locations": [[0, 4], [0, 6]] },
 *   "solver": { "tolerance": 1e-5, "start_decay": 11 }
 * }
 * @endcode
 *
 * "gravity" is one of "none", "cowling" and "full"; "quantity" one of
 * "displacement", "velocity" and "acceleration". The moment tensor is in
 * spherical coordinates at the source (r up, theta south, phi east). Receiver
 * locations are [latitude, longitude] in geographic degrees. The filter is a
 * raised-cosine taper rising between its first two corners and falling
 * between its last two; it may be left out. "attenuation", "gravity", the
 * corrections, "degrees.min", the filter, the sources' names and "solver"
 * have defaults (true, "full", false, 0, none, their numbers from 1, and the
 * solver's defaults). Relative paths are relative to the directory of the
 * input file. A member not listed here is an error, so that a misspelt one is
 * not silently ignored.
 */

#include <filesystem>
#include <string>
#include <vector>

#include "Calculation.hpp"
#include "Receiver.hpp"
#include "Source.hpp"

namespace YSpec {

/** @brief A calculation read from a file, in SI units. */
struct Input {
  std::filesystem::path model;      ///< The model deck.
  Calculation calculation;          ///< What to compute.
  std::vector<Source> sources;      ///< The sources.
  std::vector<std::string> names;   ///< The names of the sources.
  std::vector<Receiver> receivers;  ///< The receivers.
  std::filesystem::path prefix;     ///< The prefix of the output files.
};

/**
 * @brief Reads a calculation from a JSON file.
 * @throws std::runtime_error naming the file, the member and the fault.
 */
Input ReadInput(const std::filesystem::path& path);

}  // namespace YSpec
