#pragma once

/**
 * @file Model.hpp
 * @brief A spherically symmetric, transversely isotropic, anelastic earth
 * model, read from a deck and made non-dimensional.
 *
 * @details The model is a stack of layers separated by discontinuities. Within
 * a layer each property is a natural cubic spline through the deck's knots,
 * which is what Minos and the original YSpec use. All the properties of a
 * layer share its knots, so a single search finds the segment and every
 * property is then a cubic in the distance from its lower knot.
 *
 * Gravity is computed from the density splines exactly: the mass below a
 * radius is the integral of a cubic times @f$r^2@f$, a polynomial integrated
 * in closed form.
 *
 * Consecutive layers of the same phase form a region. The equations change
 * form only between regions, where solid meets fluid.
 */

#include <Interpolation/CubicSplineSystem.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <format>
#include <numbers>
#include <planetmodel/deck.hpp>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Units.hpp"

namespace YSpec {

/// The signed index type, as in GSHTrans.
using Int = std::ptrdiff_t;

/** @brief Whether a layer can support shear. */
enum class Phase { Solid, Fluid };

/** @brief Which of two layers meeting at a discontinuity is meant. */
enum class Side { Below, Above };

/**
 * @brief The properties of the model at one radius, non-dimensional.
 *
 * @details The anelastic quantities are inverse quality factors, zero where
 * the deck gives no attenuation. Following Minos, @f$A@f$ and @f$C@f$ are
 * attenuated with the effective inverse quality factor @c qA, @f$F@f$ with
 * @c qF, and @f$L@f$ and @f$N@f$ with @c qMu; for an isotropic material this
 * is the same as attenuating the bulk and shear moduli.
 */
struct Material {
  double rho;     ///< Density.
  double A;       ///< Love modulus A.
  double C;       ///< Love modulus C.
  double F;       ///< Love modulus F.
  double L;       ///< Love modulus L.
  double N;       ///< Love modulus N.
  double kappa;   ///< Voigt average bulk modulus.
  double mu;      ///< Voigt average shear modulus.
  double qKappa;  ///< Inverse quality factor in bulk.
  double qMu;     ///< Inverse quality factor in shear.
  double qA;      ///< Effective inverse quality factor of A and C.
  double qF;      ///< Effective inverse quality factor of F.
  double g;       ///< Gravitational acceleration.

  /** @brief The isotropic P-wave speed. */
  auto Alpha() const { return std::sqrt((kappa + 4 * mu / 3) / rho); }
  /** @brief The isotropic S-wave speed. */
  auto Beta() const { return std::sqrt(mu / rho); }
};

namespace ModelDetails {

/// The number of splined properties: Material less gravity.
inline constexpr auto NumberOfFields = std::size_t{12};
using Fields = std::array<double, NumberOfFields>;

inline Material MakeMaterial(const Fields& f, double g) {
  return {.rho = f[0],
          .A = f[1],
          .C = f[2],
          .F = f[3],
          .L = f[4],
          .N = f[5],
          .kappa = f[6],
          .mu = f[7],
          .qKappa = f[8],
          .qMu = f[9],
          .qA = f[10],
          .qF = f[11],
          .g = g};
}

}  // namespace ModelDetails

/**
 * @brief One layer of the model: the region between two discontinuities.
 */
class Layer {
 public:
  /// The splined properties at a knot, in Material's order less gravity.
  using Fields = ModelDetails::Fields;

  /**
   * @brief Builds a layer from its knots and the values of each property.
   * @param radii The knots, strictly increasing, at least two.
   * @param values For each knot, the properties in Material's order, less
   * gravity.
   * @param phase Whether the layer is solid or fluid.
   * @param massBelow The integral of @f$\rho r^2@f$ from the centre to the
   * layer's lower radius.
   */
  Layer(std::vector<double> radii, const std::vector<Fields>& values,
        Phase phase, double massBelow)
      : r_{std::move(radii)}, phase_{phase} {
    if (r_.size() < 2 || values.size() != r_.size()) {
      throw std::invalid_argument(
          "Layer: need at least two knots, with a value at each");
    }
    const auto n = r_.size();
    auto system = Interpolation::CubicSplineSystem{std::views::all(r_)};
    auto y = std::vector<double>(n);
    auto curvature = std::vector<double>(n);
    coefficients_.assign(n - 1, {});
    for (auto f : std::views::iota(std::size_t{0}, Fields{}.size())) {
      std::ranges::transform(values, y.begin(),
                             [f](const auto& v) { return v[f]; });
      system.Solve(y, std::span{curvature});
      for (auto k : std::views::iota(std::size_t{0}, n - 1)) {
        const auto h = r_[k + 1] - r_[k];
        auto& c = coefficients_[k];
        c[0][f] = y[k];
        c[1][f] = (y[k + 1] - y[k]) / h -
                  h * (2 * curvature[k] + curvature[k + 1]) / 6;
        c[2][f] = curvature[k] / 2;
        c[3][f] = (curvature[k + 1] - curvature[k]) / (6 * h);
      }
    }
    // The mass below each knot, from the exact integral of the density.
    massBelow_.reserve(n);
    massBelow_.push_back(massBelow);
    for (auto k : std::views::iota(std::size_t{0}, n - 1)) {
      massBelow_.push_back(massBelow_.back() +
                           MassInSegment(k, r_[k + 1] - r_[k]));
    }
  }

  /** @brief The lower radius. */
  auto Lower() const { return r_.front(); }
  /** @brief The upper radius. */
  auto Upper() const { return r_.back(); }
  /** @brief Whether the layer is solid or fluid. */
  auto GetPhase() const { return phase_; }
  /** @brief Whether the layer is fluid. */
  auto IsFluid() const { return phase_ == Phase::Fluid; }
  /** @brief The knots. */
  std::span<const double> Knots() const { return r_; }
  /** @brief The integral of @f$\rho r^2@f$ from the centre to the top. */
  auto MassBelowTop() const { return massBelow_.back(); }

  /**
   * @brief The properties at a radius.
   * @details A radius outside the layer is given the value of the end cubic
   * continued, which is what an integration step overshooting the layer by
   * rounding needs.
   */
  Material operator()(double r) const {
    const auto k = Segment(r);
    const auto t = r - r_[k];
    const auto& c = coefficients_[k];
    auto f = Fields{};
    for (auto i : std::views::iota(std::size_t{0}, f.size())) {
      f[i] = ((c[3][i] * t + c[2][i]) * t + c[1][i]) * t + c[0][i];
    }
    return ModelDetails::MakeMaterial(f, Gravity(k, t, r));
  }

  /** @brief The gravitational acceleration at a radius. */
  double Gravity(double r) const {
    const auto k = Segment(r);
    return Gravity(k, r - r_[k], r);
  }

 private:
  // Per segment, the four coefficients of the cubic in the distance from the
  // segment's lower knot, for each property.
  using Cubic = std::array<Fields, 4>;

  std::vector<double> r_;
  Phase phase_;
  std::vector<Cubic> coefficients_;
  std::vector<double> massBelow_;

  std::size_t Segment(double r) const {
    const auto it = std::ranges::upper_bound(r_, r);
    const auto k = std::distance(r_.begin(), it) - 1;
    return static_cast<std::size_t>(
        std::clamp<std::ptrdiff_t>(k, 0, std::ssize(r_) - 2));
  }

  // The integral of rho r^2 over [r_k, r_k + t], for the density cubic
  // rho = a0 + a1 t + a2 t^2 + a3 t^3 and r^2 = r_k^2 + 2 r_k t + t^2.
  double MassInSegment(std::size_t k, double t) const {
    const auto& c = coefficients_[k];
    const auto r0 = r_[k];
    auto sum = 0.0;
    auto tPower = t;  // t^(i+1)
    for (auto i : std::views::iota(0, 4)) {
      const auto a = c[i][0];
      sum += a * (r0 * r0 * tPower / (i + 1) + 2 * r0 * tPower * t / (i + 2) +
                  tPower * t * t / (i + 3));
      tPower *= t;
    }
    return sum;
  }

  double Gravity(std::size_t k, double t, double r) const {
    if (r <= 0) return 0;
    return 4 * PiG * (massBelow_[k] + MassInSegment(k, t)) / (r * r);
  }
};

/**
 * @brief A run of consecutive layers of the same phase.
 */
struct Region {
  Int first;     ///< The index of the lowest layer.
  Int last;      ///< The index of the highest layer.
  Phase phase;   ///< The phase of every layer in it.
  double lower;  ///< The lower radius.
  double upper;  ///< The upper radius.

  /** @brief Whether the region is fluid. */
  auto IsFluid() const { return phase == Phase::Fluid; }
  /** @brief The layer indices, lowest first. */
  auto Layers() const { return std::views::iota(first, last + 1); }
};

/**
 * @brief The earth model.
 */
class Model {
 public:
  /**
   * @brief Reads a model in the deck format of Minos.
   * @param path The deck.
   * @param scales The scales with which to make it non-dimensional.
   * @throws planetmodel::Error if the deck cannot be read, and
   * std::invalid_argument if it does not describe a usable model.
   */
  static Model ReadDeck(const std::filesystem::path& path,
                        Scales scales = Scales{}) {
    const auto deck = planetmodel::ReadDeck(path, planetmodel::MINEOS());
    const auto header =
        planetmodel::ReadMineosHeader(deck.Header(), deck.Source());
    return Model(deck, header.tref, scales);
  }

  /**
   * @brief Builds the model from a deck already read.
   * @param deck The deck, in SI units, in either the six or the nine column
   * format.
   * @param referencePeriod The period in seconds at which the deck's
   * velocities apply.
   * @param scales The scales with which to make it non-dimensional.
   */
  Model(const planetmodel::Deck& deck, double referencePeriod, Scales scales)
      : scales_{scales},
        referenceFrequency_{2 * std::numbers::pi / referencePeriod /
                            scales.Frequency()} {
    const auto anisotropic = deck.HasColumn("vpv");
    const auto& radius = deck.Radius();
    const auto& rho = deck.Column("rho");
    const auto& vpv = deck.Column(anisotropic ? "vpv" : "vp");
    const auto& vsv = deck.Column(anisotropic ? "vsv" : "vs");
    const auto& vph = deck.Column(anisotropic ? "vph" : "vp");
    const auto& vsh = deck.Column(anisotropic ? "vsh" : "vs");
    const auto& qKappa = deck.Column("qkappa");
    const auto& qMu = deck.Column("qmu");

    const auto velocity = scales.Velocity();
    auto inverse = [](double q) { return q > 0 ? 1 / q : 0.0; };
    auto massBelow = 0.0;
    for (auto [begin, end] : deck.Layers()) {
      auto r = std::vector<double>{};
      auto values = std::vector<Layer::Fields>{};
      auto shear = 0.0;
      for (auto i : std::views::iota(begin, end)) {
        const auto density = rho[i] / scales.Density();
        auto modulus = [&](double v) {
          const auto u = v / velocity;
          return density * u * u;
        };
        const auto A = modulus(vph[i]);
        const auto C = modulus(vpv[i]);
        const auto L = modulus(vsv[i]);
        const auto N = modulus(vsh[i]);
        const auto eta = anisotropic ? deck.Column("eta")[i] : 1.0;
        const auto F = eta * (A - 2 * L);
        const auto kappa = (C + 4 * A - 4 * N + 4 * F) / 9;
        const auto mu = (C + A + 6 * L + 5 * N - 2 * F) / 15;
        const auto qk = inverse(qKappa[i]);
        const auto qm = inverse(qMu[i]);
        // Minos's effective inverse quality factors for A, C and F.
        const auto ratio = 4 * mu / (3 * kappa + 4 * mu);
        const auto qA = (1 - ratio) * qk + ratio * qm;
        const auto qF = ((1 - ratio) * qk - ratio * qm / 2) / (1 - 1.5 * ratio);
        r.push_back(radius[i] / scales.Length());
        values.push_back({density, A, C, F, L, N, kappa, mu, qk, qm, qA, qF});
        shear = std::max({shear, L, N});
      }
      const auto phase = shear > 0 ? Phase::Solid : Phase::Fluid;
      layers_.emplace_back(std::move(r), values, phase, massBelow);
      massBelow = layers_.back().MassBelowTop();
    }
    if (layers_.empty()) throw std::invalid_argument("Model: empty deck");

    // Group the layers into regions of one phase.
    for (auto i : std::views::iota(Int{0}, std::ssize(layers_))) {
      const auto& layer = layers_[i];
      if (regions_.empty() || regions_.back().phase != layer.GetPhase()) {
        regions_.push_back({.first = i,
                            .last = i,
                            .phase = layer.GetPhase(),
                            .lower = layer.Lower(),
                            .upper = layer.Upper()});
      } else {
        regions_.back().last = i;
        regions_.back().upper = layer.Upper();
      }
    }
  }

  /** @brief The scales with which the model was made non-dimensional. */
  const Scales& GetScales() const { return scales_; }

  /** @brief The frequency at which the moduli are given, non-dimensional. */
  auto ReferenceFrequency() const { return referenceFrequency_; }

  /** @brief The number of layers. */
  auto NumberOfLayers() const { return std::ssize(layers_); }
  /** @brief A layer, numbered from the centre. */
  const Layer& GetLayer(Int i) const { return layers_[i]; }
  /** @brief The layers, from the centre. */
  std::span<const Layer> Layers() const { return layers_; }

  /** @brief The regions of one phase, from the centre. */
  std::span<const Region> Regions() const { return regions_; }

  /** @brief The radius of the centre (zero unless the model is hollow). */
  auto InnerRadius() const { return layers_.front().Lower(); }
  /** @brief The radius of the surface. */
  auto SurfaceRadius() const { return layers_.back().Upper(); }

  /**
   * @brief The layer containing a radius.
   * @param r The radius.
   * @param side At a discontinuity, which of the two layers to return.
   * @throws std::out_of_range if the radius is outside the model.
   */
  Int LayerIndex(double r, Side side = Side::Below) const {
    if (r < InnerRadius() || r > SurfaceRadius()) {
      throw std::out_of_range(std::format("radius {} is outside the model", r));
    }
    auto it = side == Side::Below
                  ? std::ranges::lower_bound(layers_, r, {}, &Layer::Upper)
                  : std::ranges::upper_bound(layers_, r, {}, &Layer::Upper);
    if (it == layers_.end()) --it;
    return std::distance(layers_.begin(), it);
  }

  /** @brief The index of the region containing a layer. */
  Int RegionIndex(Int layer) const {
    const auto it = std::ranges::find_if(
        regions_, [layer](const auto& region) { return layer <= region.last; });
    return std::distance(regions_.begin(), it);
  }

  /** @brief The properties at a radius, in the layer given. */
  Material operator()(Int layer, double r) const { return layers_[layer](r); }

 private:
  Scales scales_;
  double referenceFrequency_;
  std::vector<Layer> layers_;
  std::vector<Region> regions_;
};

}  // namespace YSpec
