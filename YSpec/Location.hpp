#pragma once

/**
 * @file Location.hpp
 * @brief A point on the surface of the earth.
 */

namespace YSpec {

/** @brief A point on the surface, in geographic degrees. */
struct Location {
  double latitude = 0;   ///< Geographic latitude, degrees north.
  double longitude = 0;  ///< Longitude, degrees east.
};

}  // namespace YSpec
