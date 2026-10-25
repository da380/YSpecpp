#pragma once

/**
 * @file Location.hpp
 * @brief A point on the surface of the earth, and lines of points along great
 * circles.
 *
 * @details Locations are in geographic latitude when an explicit
 * HydrostaticFigure is supplied. Model-less GreatCircle overloads interpret
 * locations as geocentric and use spherical geometry.
 */

#include <cstddef>
#include <vector>

#include "HydrostaticFigure.hpp"

namespace YSpec {

/** @brief A point on the surface. Latitude and longitude are in degrees;
 * latitude is geographic when a computed figure is supplied and geocentric for
 * model-less spherical helpers. */
struct Location {
  double latitude =
      0;  ///< Degrees north: geographic with a figure, otherwise geocentric.
  double longitude = 0;  ///< Longitude, degrees east.
};

/**
 * @brief Points equally spaced along the shorter great circle between two
 * points, both included.
 * @param from The first point.
 * @param to The last point.
 * @param count The number of points; one gives @p from alone.
 * @throws std::invalid_argument if the count is zero, or if the points are
 * antipodal, so that the great circle is undefined.
 */
std::vector<Location> GreatCircle(const Location& from, const Location& to,
                                  std::size_t count);
/** Great circle between geographic locations using the supplied model figure.
 */
std::vector<Location> GreatCircle(const Location& from, const Location& to,
                                  std::size_t count,
                                  const HydrostaticFigure& figure);

/**
 * @brief Points along the great circle leaving a point at an azimuth, at
 * distances equally spaced between two limits, both included.
 * @param from The point the distances are measured from.
 * @param azimuth The azimuth at @p from, degrees clockwise from north.
 * @param minDistance The distance of the first point (degrees).
 * @param maxDistance The distance of the last point (degrees).
 * @param count The number of points; one gives the point at @p minDistance.
 * @throws std::invalid_argument if the count is zero.
 */
std::vector<Location> GreatCircle(const Location& from, double azimuth,
                                  double minDistance, double maxDistance,
                                  std::size_t count);
/** Great circle from a geographic location using the supplied model figure. */
std::vector<Location> GreatCircle(const Location& from, double azimuth,
                                  double minDistance, double maxDistance,
                                  std::size_t count,
                                  const HydrostaticFigure& figure);

}  // namespace YSpec
