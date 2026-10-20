#pragma once

/**
 * @file Receiver.hpp
 * @brief A receiver.
 */

#include "Location.hpp"

namespace YSpec {

/** @brief A receiver: where on the surface, and how deep (m). */
struct Receiver {
  Location location;  ///< The point on the surface above it.
  double depth = 0;   ///< The depth (m).
};

}  // namespace YSpec
