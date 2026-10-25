#pragma once

// Spherical geometry in geocentric coordinates. Model-derived conversion
// between geographic and geocentric latitude is owned by HydrostaticFigure.

#include <array>
#include <cmath>
#include <numbers>

namespace YSpec::Sphere {

inline constexpr auto Degree = std::numbers::pi / 180;

// Convert a geocentric latitude in degrees to colatitude in radians.
inline double Colatitude(double geocentricLatitude) {
  return std::numbers::pi / 2 - geocentricLatitude * Degree;
}

// Convert colatitude in radians to geocentric latitude in degrees.
inline double Latitude(double colatitude) {
  return 90 - colatitude / Degree;
}

using Vector = std::array<double, 3>;

// The unit vector at a colatitude and longitude (radians).
inline Vector Unit(double theta, double phi) {
  return {std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi),
          std::cos(theta)};
}

}  // namespace YSpec::Sphere
