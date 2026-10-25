#pragma once

/**
 * @file HydrostaticFigure.hpp
 * @brief The first-order hydrostatic figure of a slowly rotating model.
 */

#include <span>
#include <vector>

namespace YSpec {

class Model;
struct Calculation;

/// IERS nominal mean angular velocity of the Earth (rad/s).
inline constexpr double EarthRotationRate = 7.292115e-5;

/** A sampled, first-order hydrostatic ellipticity profile and its conversion. */
class HydrostaticFigure {
 public:
  /** One radial sample and derivative of the first-order figure profile. */
  struct Node {
    /** Radius of the profile sample (m). */
    double radius;            ///< Radius (m).
    /** First-order ellipticity at this radius. */
    double ellipticity;       ///< First-order figure ellipticity at the radius.
    /** Radial derivative of ellipticity (m^-1). */
    double radialDerivative;  ///< d(ellipticity)/dr (m^-1).
  };

  /** A spherical figure: supplied coordinates are already geocentric. */
  static HydrostaticFigure Spherical(double rotationRate = 0.0);

  /**
   * @brief The computed figure from a complete spherical model.
   * @param model Density model with a regular centre.
   * @param rotationRate Uniform angular rotation rate (rad/s).
   * @param tolerance Relative adaptive integration tolerance.
   */
  static HydrostaticFigure FromModel(const Model& model,
                                     double rotationRate = EarthRotationRate,
                                     double tolerance = 1e-9);

  /** Whether a model-derived coordinate correction is enabled. */
  bool Enabled() const { return enabled_; }
  /** Uniform rotation rate used to compute this figure (rad/s). */
  double RotationRate() const { return rotationRate_; }
  /** Ellipticity at the model surface. */
  double SurfaceEllipticity() const { return surfaceEllipticity_; }
  /** Sampled radial ellipticity profile; a diagnostic, not a coordinate map. */
  std::span<const Node> Nodes() const { return nodes_; }

  /** Convert geographic latitude to geocentric latitude (degrees). */
  double GeocentricLatitude(double geographicLatitude) const;
  /** Convert geocentric latitude to geographic latitude (degrees). */
  double GeographicLatitude(double geocentricLatitude) const;
  /** The first-order figure ellipticity at a radius (m). */
  double Ellipticity(double radius) const;

 private:
  bool enabled_ = false;
  double rotationRate_ = 0;
  double surfaceEllipticity_ = 0;
  std::vector<Node> nodes_;
};

/** Resolve the calculation's model-aware coordinate convention. */
HydrostaticFigure FigureFor(const Model& model, const Calculation& calculation);

}  // namespace YSpec
