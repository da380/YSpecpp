#include <YSpec/HydrostaticFigure.hpp>
#include <YSpec/Calculation.hpp>

#include <YSpec/Model.hpp>
#include <YSpec/Units.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <ranges>
#include <stdexcept>

#include "RungeKutta.hpp"

namespace YSpec {

HydrostaticFigure HydrostaticFigure::Spherical(double rotationRate) {
  if (!(rotationRate >= 0) || !std::isfinite(rotationRate)) {
    throw std::invalid_argument("rotation rate must be finite and non-negative");
  }
  HydrostaticFigure figure;
  figure.rotationRate_ = rotationRate;
  figure.nodes_ = {{0, 0, 0}};
  return figure;
}

HydrostaticFigure HydrostaticFigure::FromModel(const Model& model,
                                                double rotationRate,
                                                double tolerance) {
  if (!(rotationRate >= 0) || !std::isfinite(rotationRate)) {
    throw std::invalid_argument("rotation rate must be finite and non-negative");
  }
  if (!(tolerance > 0 && tolerance < 1) || !std::isfinite(tolerance)) {
    throw std::invalid_argument("figure tolerance must lie between zero and one");
  }
  if (std::abs(model.InnerRadius()) > 1e-12 * model.SurfaceRadius()) {
    throw std::invalid_argument(
        "hydrostatic figure requires a model extending to a regular centre");
  }
  const auto& scales = model.GetScales();
  const auto radius = model.SurfaceRadius();
  if (!(radius > 0)) throw std::invalid_argument("model surface radius must be positive");

  HydrostaticFigure figure;
  figure.enabled_ = true;
  figure.rotationRate_ = rotationRate;
  const auto surfaceMaterial = model(model.NumberOfLayers() - 1, radius);
  const auto surfaceGravity = surfaceMaterial.g * scales.Acceleration();
  const auto surfaceRadius = radius * scales.Length();
  if (!(surfaceGravity > 0) || !std::isfinite(surfaceGravity)) {
    throw std::invalid_argument("hydrostatic figure requires positive surface gravity");
  }
  const auto m = rotationRate * rotationRate * surfaceRadius / surfaceGravity;
  if (!(m < 0.1)) {
    throw std::invalid_argument(
        "first-order hydrostatic figure requires slow rotation (m < 0.1)");
  }

  // Integrate in x=r/R. The regular central solution is epsilon=1, d/dx=0;
  // starting at 1e-6 R bounds the neglected centre-series term at O(1e-12)
  // for smooth density profiles.
  using State = ScaledVector<2>;
  auto state = State{};
  state.value = {1.0, 0.0};
  auto step = 1e-3;
  const auto start = 1e-6;
  auto ratio = [&](Int layer, double x) {
    if (x <= 0) return 1.0;
    const auto r = x * radius;
    const auto material = model(layer, r);
    const auto rho = material.rho * scales.Density();
    const auto gravity = material.g * scales.Acceleration();
    if (!(rho > 0) || !(gravity > 0)) {
      throw std::invalid_argument(
          "hydrostatic figure requires positive density and enclosed mass");
    }
    const auto rSI = r * scales.Length();
    const auto mean =
        3 * gravity /
        (4 * std::numbers::pi * GravitationalConstant * rSI);
    return rho / mean;
  };
  auto derivative = [&](Int layer, double x, const Vector<2>& y) {
    const auto q = ratio(layer, x);
    return Vector<2>{y[1], -6 * q / x * y[1] -
                               6 * (q - 1) / (x * x) * y[0]};
  };
  const auto integrator = RungeKutta<DormandPrince54>(
      IntegrationOptions{.adaptive = true, .tolerance = tolerance});
  figure.nodes_.push_back({0, 1, 0});
  auto x = start;
  const auto first = model(Int{0}, x * radius);
  if (!(first.rho > 0)) {
    throw std::invalid_argument("hydrostatic figure requires positive central density");
  }
  for (auto layer : std::views::iota(Int{0}, model.NumberOfLayers())) {
    const auto& material = model.GetLayer(layer);
    for (const auto knot : material.Knots()) {
      auto end = knot / radius;
      if (end <= x) continue;
      integrator.Integrate(
          [&](double at, const Vector<2>& y) { return derivative(layer, at, y); },
          x, end, state, step);
      x = end;
      const auto value = Rescale(state.value[0], state.logScale).real();
      const auto derivativeValue =
          Rescale(state.value[1], state.logScale).real() / surfaceRadius;
      const auto node =
          Node{knot * scales.Length(), value, derivativeValue};
      if (figure.nodes_.back().radius == node.radius) {
        figure.nodes_.back() = node;
      } else {
        figure.nodes_.push_back(node);
      }
    }
  }
  if (x < 1) {
    const auto last = model.NumberOfLayers() - 1;
    integrator.Integrate(
        [&](double at, const Vector<2>& y) { return derivative(last, at, y); },
        x, 1.0, state, step);
  }
  const auto e = Rescale(state.value[0], state.logScale).real();
  const auto de = Rescale(state.value[1], state.logScale).real();
  const auto denominator = de + 2 * e;
  if (!(denominator > 0) || !std::isfinite(denominator)) {
    throw std::runtime_error("hydrostatic figure surface normalization failed");
  }
  const auto normalization = 2.5 * m / denominator;
  figure.surfaceEllipticity_ = normalization * e;
  for (auto& node : figure.nodes_) {
    node.ellipticity *= normalization;
    node.radialDerivative *= normalization;
  }
  if (figure.nodes_.back().radius == surfaceRadius) {
    figure.nodes_.back().ellipticity = figure.surfaceEllipticity_;
    figure.nodes_.back().radialDerivative = normalization * de / surfaceRadius;
  } else {
    figure.nodes_.push_back(
        {surfaceRadius, figure.surfaceEllipticity_, normalization * de / surfaceRadius});
  }
  if (!(1 - 2 * figure.surfaceEllipticity_ > 0)) {
    throw std::invalid_argument(
        "first-order hydrostatic figure exceeds its positive latitude-map range");
  }
  return figure;
}

double HydrostaticFigure::GeocentricLatitude(double geographicLatitude) const {
  if (!enabled_) return geographicLatitude;
  const auto ratio = 1 - 2 * surfaceEllipticity_;
  const auto angle = geographicLatitude * std::numbers::pi / 180;
  return std::atan2(ratio * std::sin(angle), std::cos(angle)) * 180 /
         std::numbers::pi;
}

double HydrostaticFigure::GeographicLatitude(double geocentricLatitude) const {
  if (!enabled_) return geocentricLatitude;
  const auto ratio = 1 - 2 * surfaceEllipticity_;
  const auto angle = geocentricLatitude * std::numbers::pi / 180;
  return std::atan2(std::sin(angle), ratio * std::cos(angle)) * 180 /
         std::numbers::pi;
}

double HydrostaticFigure::Ellipticity(double radius) const {
  if (!std::isfinite(radius)) {
    throw std::invalid_argument("figure radius must be finite");
  }
  if (!enabled_ || nodes_.empty()) return 0;
  if (radius <= nodes_.front().radius) return nodes_.front().ellipticity;
  if (radius >= nodes_.back().radius) return nodes_.back().ellipticity;
  const auto upper = std::ranges::lower_bound(nodes_, radius, {}, &Node::radius);
  const auto& b = *upper;
  const auto& a = *(upper - 1);
  const auto width = b.radius - a.radius;
  const auto t = (radius - a.radius) / width;
  const auto t2 = t * t;
  const auto t3 = t2 * t;
  return (2 * t3 - 3 * t2 + 1) * a.ellipticity +
         (t3 - 2 * t2 + t) * width * a.radialDerivative +
         (-2 * t3 + 3 * t2) * b.ellipticity +
         (t3 - t2) * width * b.radialDerivative;
}

HydrostaticFigure FigureFor(const Model& model,
                            const Calculation& calculation) {
  if (!calculation.latitudeCorrection) {
    return HydrostaticFigure::Spherical(calculation.rotationRate);
  }
  return HydrostaticFigure::FromModel(model, calculation.rotationRate);
}

}  // namespace YSpec
