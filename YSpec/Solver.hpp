#pragma once

/**
 * @file Solver.hpp
 * @brief The response at a receiver to unit jumps at a source, for one degree
 * and one complex frequency.
 *
 * @details The source and the receiver lie in the same solid region. For each
 * system a solution regular at the bottom (integrated upwards from a starting
 * level) and one satisfying the conditions at the top (integrated downwards
 * from the surface) are combined with the jump at the source, in the manner of
 * Friederich & Dalkolmo (1995) and Al-Attar & Woodhouse (2008).
 *
 * The integration starts at the deepest point at which the solution can
 * matter: below the turning point of the slowest wave at the given degree and
 * frequency the solution decays, and the WKB estimate of that decay is
 * accumulated downwards until it reaches a given number of e-folds.
 *
 * Integrations restart at every discontinuity of the model, and the systems
 * change from solid to fluid form and back wherever the phase does, so a fluid
 * inner core, an ocean, or any other sequence of solid and fluid regions needs
 * no special treatment.
 *
 * A solver holds no state that a solve changes, so one solver can be shared
 * by any number of threads.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <format>
#include <limits>
#include <stdexcept>
#include <vector>

#include "Equations.hpp"
#include "Model.hpp"
#include "RungeKutta.hpp"

namespace YSpec {

/** @brief A radius in a specified layer. */
struct Point {
  Int layer;  ///< The layer, which settles which side of a discontinuity.
  double r;   ///< The radius.
};

/** @brief Settings of the radial solver. */
struct SolverOptions {
  /// The integration settings. The default tolerance makes seismograms
  /// accurate to a few parts in 10^4 of their peak, about ten times better
  /// than the original YSpec; 1e-4 matches it at less cost.
  IntegrationOptions integration = {.tolerance = 1e-5};
  /// The decay, in e-folds, below a turning point at which to start.
  double startDecay = 11.0;
  /// The first step in a layer is this many times @f$\alpha/\omega@f$, for the
  /// largest P-wave speed in the layer (with fixed steps, every step is).
  double stepFactor = 0.25;
};

/** @brief The response of the toroidal system: W for the two unit jumps. */
using ToroidalResponse = std::array<Complex, 2>;

/** @brief The response of the radial system: U for the two unit jumps. */
using RadialResponse = std::array<Complex, 2>;

/**
 * @brief The response of the spheroidal system: U, V, the potential and the
 * gravity-flux variable (the last two zero unless self-gravitating), each for
 * the four unit jumps.
 */
using SpheroidalResponse = std::array<std::array<Complex, 4>, 4>;

namespace SolverDetails {

// The spheroidal equations without self-gravitation: 4x4 in a solid, 2x2 in a
// fluid.
template <bool Cowling>
struct Equations4 {
  static constexpr std::size_t SolidSize = 5;
  static constexpr std::size_t FluidSize = 2;
  static constexpr std::size_t BSize = 4;
  // With the receiver above the source the solution is +Solution/Delta.
  static constexpr double Sign = 1;

  static Block4 Solid(const Context& c, const ComplexMaterial& m) {
    return SolidBlock4(c, m, Cowling);
  }
  static Block2 Fluid(const Context& c, const ComplexMaterial& m) {
    return FluidBlock2(c, m, Cowling);
  }
  static Vector<5> SolidStart(const Context& c, const ComplexMaterial& m) {
    return SolidStart4(c, m);
  }
  static Vector<2> FluidStart(const Context& c, const ComplexMaterial& m) {
    return FluidStart2(c, m);
  }
  static Vector<5> SolidMinors(const Block4& b, const Vector<5>& m) {
    return MinorDerivative(b, m);
  }
  static Vector<2> FluidMinors(const Block2& b, const Vector<2>& m) {
    return Derivative(b, m);
  }
  // The jumps in (U, V, P, S).
  static Vector<4> Jump(int j) {
    auto s = Vector<4>{};
    s[j] = 1;
    return s;
  }
  static Vector<4> B(const Vector<5>& m, const Vector<4>& s) {
    return BStart(m, s);
  }
  // U, V, the potential and the gravity-flux variable.
  static std::array<Complex, 4> Components(const Vector<4>& y) {
    return {y[0], y[1], Complex{}, Complex{}};
  }
};

// The self-gravitating spheroidal equations: 6x6 in a solid, 4x4 in a fluid.
struct Equations6 {
  static constexpr std::size_t SolidSize = 14;
  static constexpr std::size_t FluidSize = 5;
  static constexpr std::size_t BSize = 14;
  static constexpr double Sign = -1;

  static Block6 Solid(const Context& c, const ComplexMaterial& m) {
    return SolidBlock6(c, m);
  }
  static Block4 Fluid(const Context& c, const ComplexMaterial& m) {
    return FluidBlock4(c, m);
  }
  static Vector<14> SolidStart(const Context& c, const ComplexMaterial& m) {
    return SolidStart6(c, m);
  }
  static Vector<5> FluidStart(const Context& c, const ComplexMaterial& m) {
    return FluidStart4(c, m);
  }
  static Vector<14> SolidMinors(const Block6& b, const Vector<14>& m) {
    return MinorDerivative(b, m);
  }
  static Vector<5> FluidMinors(const Block4& b, const Vector<5>& m) {
    return MinorDerivative(b, m);
  }
  // The jumps in (U, V, P, S), placed in (U, V, phi, P, S, psi): the
  // potential and its flux are continuous at the source.
  static Vector<6> Jump(int j) {
    constexpr auto position = std::array{0, 1, 3, 4};
    auto s = Vector<6>{};
    s[position[j]] = 1;
    return s;
  }
  static Vector<14> B(const Vector<14>& m, const Vector<6>& s) {
    return BStart(m, s);
  }
  static std::array<Complex, 4> Components(const Vector<6>& y) {
    return {y[0], y[1], y[2], y[5]};
  }
};

// Minors integrated through regions of either phase.
template <class E>
struct MinorState {
  ScaledVector<E::SolidSize> solid;
  ScaledVector<E::FluidSize> fluid;
  bool inFluid = false;

  void ChangePhase(bool toFluid) {
    if (toFluid == inFluid) return;
    if (toFluid) {
      fluid = {SolidToFluid(solid.value), solid.logScale};
      fluid.Normalise();
    } else {
      solid = {FluidToSolid(fluid.value), fluid.logScale};
      solid.Normalise();
    }
    inFluid = toFluid;
  }
};

}  // namespace SolverDetails

/**
 * @brief Solves the radial equations for a fixed source and receiver radius.
 * @tparam Tableau The Runge-Kutta method.
 */
template <ButcherTableau Tableau = DormandPrince54>
class RadialSolver {
 public:
  /**
   * @param model The earth model.
   * @param sourceRadius The radius of the source, non-dimensional.
   * @param receiverRadius The radius of the receivers, non-dimensional.
   * @param gravity How gravity enters.
   * @param attenuation Whether to include attenuation.
   * @param options Integration settings.
   * @throws std::invalid_argument unless the source and receiver lie in the
   * same solid region, or if the options are out of range. A source or
   * receiver on a boundary between solid and fluid is taken to be on the
   * solid side.
   */
  RadialSolver(const Model& model, double sourceRadius, double receiverRadius,
               Gravity gravity, bool attenuation = true,
               SolverOptions options = {})
      : model_{&model},
        gravity_{gravity},
        attenuation_{attenuation},
        options_{options},
        integrator_{options.integration},
        source_{SolidSide(model, sourceRadius)},
        receiver_{SolidSide(model, receiverRadius)} {
    if (!(options.integration.tolerance > 0 &&
          options.integration.tolerance < 1)) {
      throw std::invalid_argument("the tolerance must lie between 0 and 1");
    }
    if (!(options.startDecay > 0)) {
      throw std::invalid_argument("the start decay must be positive");
    }
    if (!(options.stepFactor > 0)) {
      throw std::invalid_argument("the step factor must be positive");
    }
    region_ = model.RegionIndex(source_.layer);
    const auto& region = model.Regions()[region_];
    if (region.IsFluid()) {
      throw std::invalid_argument(
          "the source is in a fluid region of the model (an ocean?); it must "
          "be in a solid");
    }
    if (model.GetLayer(receiver_.layer).IsFluid()) {
      throw std::invalid_argument(
          "the receivers are in a fluid region of the model (an ocean?); they "
          "must be in a solid");
    }
    if (model.RegionIndex(receiver_.layer) != region_) {
      throw std::invalid_argument(
          "the source and the receivers must lie in the same solid region of "
          "the model");
    }
    for (const auto& layer : model.Layers()) {
      auto alpha = 0.0;
      for (auto r : layer.Knots()) alpha = std::max(alpha, layer(r).Alpha());
      maxAlpha_.push_back(alpha);
    }
  }

  /** @brief The model. */
  const Model& GetModel() const { return *model_; }
  /** @brief The source. */
  Point Source() const { return source_; }
  /** @brief The receiver. */
  Point Receiver() const { return receiver_; }
  /** @brief How gravity enters. */
  Gravity GetGravity() const { return gravity_; }
  /** @brief Whether attenuation is included. */
  bool Attenuation() const { return attenuation_; }

  /**
   * @brief The toroidal response for degree @p l >= 1.
   */
  ToroidalResponse Toroidal(int l, Complex omega) const {
    const auto c = Context(*model_, l, omega, attenuation_);
    const auto& region = model_->Regions()[region_];
    auto rhs = [&c](Int layer) {
      return [&c, layer](double r, const Vector<2>& y) {
        return Derivative(ToroidalBlock(c, c(layer, r)), y);
      };
    };
    const auto lower = Lower();
    const auto start = StartPoint(l, omega.real(), lower, region_);
    auto h = 0.0;
    auto low = ScaledVector<2>{};
    auto from = start;
    if (start.r == region.lower && region.lower > 0) {
      // Fluid below: the traction vanishes.
      low.value = {1.0, 0.0};
    } else {
      from = StartRadius(start, omega, h);
      low.value = ToroidalStart(c, c(from.layer, from.r));
      low.Normalise();
    }
    Walk(from, lower, low, h, omega, rhs);
    auto up = ScaledVector<2>{{1.0, 0.0}};
    auto hUp = 0.0;
    Walk(Point{region.last, region.upper}, Upper(), up, hUp, omega, rhs);
    return TwoByTwo(low, up, h, omega, rhs);
  }

  /**
   * @brief The radial (l = 0) response.
   */
  RadialResponse Radial(Complex omega) const {
    const auto c = Context(*model_, 0, omega, attenuation_);
    const auto gravity = gravity_ != Gravity::None;
    auto rhs = [&c, gravity](Int layer) {
      return [&c, layer, gravity](double r, const Vector<2>& y) {
        return Derivative(RadialBlock(c, c(layer, r), gravity), y);
      };
    };
    auto h = 0.0;
    const auto start = StartRadius(Point{0, model_->InnerRadius()}, omega, h);
    auto low =
        ScaledVector<2>{RadialStart(c, c(start.layer, start.r), gravity)};
    low.Normalise();
    Walk(start, Lower(), low, h, omega, rhs);
    auto up = ScaledVector<2>{{1.0, 0.0}};
    auto hUp = 0.0;
    Walk(Point{model_->NumberOfLayers() - 1, model_->SurfaceRadius()}, Upper(),
         up, hUp, omega, rhs);
    return TwoByTwo(low, up, h, omega, rhs);
  }

  /**
   * @brief The spheroidal response for degree @p l >= 1.
   */
  SpheroidalResponse Spheroidal(int l, Complex omega) const {
    switch (gravity_) {
      case Gravity::None:
        return SpheroidalImpl<SolverDetails::Equations4<false>>(l, omega);
      case Gravity::Cowling:
        return SpheroidalImpl<SolverDetails::Equations4<true>>(l, omega);
      case Gravity::Full:
        return SpheroidalImpl<SolverDetails::Equations6>(l, omega);
    }
    return {};
  }

  /**
   * @brief Where the upward integration starts for a degree and frequency.
   * @param l The degree.
   * @param omega The real part of the frequency.
   * @param top The point the integration goes to.
   * @param lowestRegion The deepest region to look in.
   * @return A knot below @p top, or the bottom of @p lowestRegion if the
   * solution does not decay enough above it.
   */
  Point StartPoint(int l, double omega, Point top, Int lowestRegion) const {
    const auto& lowest = model_->Regions()[lowestRegion];
    const auto bottom = Point{lowest.first, lowest.lower};
    if (l == 0 || omega <= 0) return bottom;
    const auto p2 = static_cast<double>(l) * (l + 1) / (omega * omega);
    // The square of the vertical slowness of the slowest wave, negative
    // where the solution is evanescent.
    auto q2 = [&](Int layer, double r) {
      if (r <= 0) return -std::numeric_limits<double>::infinity();
      const auto& L = model_->GetLayer(layer);
      const auto m = L(r);
      const auto v = L.IsFluid() ? m.Alpha() : m.Beta();
      return 1 / (v * v) - p2 / (r * r);
    };
    auto sum = 0.0;
    auto region = model_->RegionIndex(top.layer);
    for (auto layer = top.layer; layer >= lowest.first; --layer) {
      if (model_->RegionIndex(layer) != region) {
        region = model_->RegionIndex(layer);
        sum = 0;
      }
      const auto knots = model_->GetLayer(layer).Knots();
      for (auto k = std::ssize(knots) - 1; k > 0; --k) {
        const auto upper = knots[k];
        const auto lower = knots[k - 1];
        if (lower >= top.r) continue;
        const auto qu = q2(layer, upper);
        const auto ql = q2(layer, lower);
        if (qu < 0 && ql < 0) {
          const auto su = std::isinf(qu) ? std::sqrt(-ql) : std::sqrt(-qu);
          const auto sl = std::isinf(ql) ? su : std::sqrt(-ql);
          sum += 0.5 * omega * (su + sl) * (std::min(upper, top.r) - lower);
        } else {
          sum = 0;
        }
        if (sum > options_.startDecay) return Point{layer, lower};
      }
    }
    return bottom;
  }

 private:
  const Model* model_;
  Gravity gravity_;
  bool attenuation_;
  SolverOptions options_;
  RungeKutta<Tableau> integrator_;
  Point source_;
  Point receiver_;
  Int region_;
  std::vector<double> maxAlpha_;

  Point Lower() const { return source_.r <= receiver_.r ? source_ : receiver_; }
  Point Upper() const { return source_.r <= receiver_.r ? receiver_ : source_; }
  bool ReceiverAbove() const { return receiver_.r >= source_.r; }

  // The first step in a layer.
  double FirstStep(Int layer, Complex omega) const {
    return options_.stepFactor * maxAlpha_[layer] / std::abs(omega);
  }

  // A radius in its layer, taking the solid side of a solid-fluid boundary.
  static Point SolidSide(const Model& model, double r) {
    const auto below = model.LayerIndex(r, Side::Below);
    const auto above = model.LayerIndex(r, Side::Above);
    if (model.GetLayer(below).IsFluid() && !model.GetLayer(above).IsFluid()) {
      return {above, r};
    }
    return {below, r};
  }

  // A start at the centre is moved off it, since the starting solutions
  // divide by the radius.
  Point StartRadius(Point start, Complex omega, double& h) const {
    h = FirstStep(start.layer, omega);
    if (start.r > 0) return start;
    const auto knots = model_->GetLayer(start.layer).Knots();
    return {start.layer, std::min(0.1 * h, 0.1 * knots[1])};
  }

  // Integrates a system of one form from one point to another through the
  // layers between them. makeRhs(layer) gives the right-hand side in a layer.
  template <std::size_t N, class MakeRhs>
  void Walk(Point from, Point to, ScaledVector<N>& y, double& h, Complex omega,
            MakeRhs&& makeRhs) const {
    WalkLayers(from, to, h, omega,
               [&](Int layer, double r0, double r1, double& step) {
                 integrator_.Integrate(makeRhs(layer), r0, r1, y, step);
               });
  }

  // Calls integrate(layer, r0, r1, h) for each layer from one point to
  // another, in order, with the interval of the layer between them.
  template <class IntegrateLayer>
  void WalkLayers(Point from, Point to, double& h, Complex omega,
                  IntegrateLayer&& integrate) const {
    const auto up = to.layer >= from.layer;
    for (auto layer = from.layer;; layer += up ? 1 : -1) {
      const auto& L = model_->GetLayer(layer);
      const auto r0 =
          layer == from.layer ? from.r : (up ? L.Lower() : L.Upper());
      const auto r1 = layer == to.layer ? to.r : (up ? L.Upper() : L.Lower());
      if (!options_.integration.adaptive || h <= 0) {
        h = FirstStep(layer, omega);
      }
      integrate(layer, r0, r1, h);
      if (layer == to.layer) break;
    }
  }

  // Combines the lower and upper solutions of a 2x2 system: each is at the
  // nearer of the source and receiver on entry. h is the step in use.
  template <class MakeRhs>
  std::array<Complex, 2> TwoByTwo(ScaledVector<2>& low, ScaledVector<2>& up,
                                  double& h, Complex omega,
                                  MakeRhs&& makeRhs) const {
    // The solution at the source is continued to the receiver, where the
    // Wronskian is formed (it is constant, so any radius would do).
    auto& atSource = ReceiverAbove() ? low : up;
    const auto& atReceiver = ReceiverAbove() ? up : low;
    auto continued = atSource;
    Walk(source_, receiver_, continued, h, omega, makeRhs);
    const auto& lowR = ReceiverAbove() ? continued : atReceiver;
    const auto& upR = ReceiverAbove() ? atReceiver : continued;
    const auto delta = ScaledComplex{
        upR.value[1] * lowR.value[0] - upR.value[0] * lowR.value[1],
        upR.logScale + lowR.logScale};
    // The jumps (1, 0) and (0, 1).
    const auto b = std::array{
        atSource[1], ScaledComplex{-atSource.value[0], atSource.logScale}};
    auto result = std::array<Complex, 2>{};
    for (auto i : {0, 1}) {
      const auto y = atReceiver[0] * b[i] / delta;
      result[i] = -y.Get();
    }
    return result;
  }

  template <class E>
  SpheroidalResponse SpheroidalImpl(int l, Complex omega) const {
    using State = SolverDetails::MinorState<E>;
    const auto c = Context(*model_, l, omega, attenuation_);

    // Minors, through layers of either phase.
    auto walkMinors = [&](Point from, Point to, State& s, double& h) {
      WalkLayers(from, to, h, omega,
                 [&](Int layer, double r0, double r1, double& step) {
                   s.ChangePhase(model_->GetLayer(layer).IsFluid());
                   if (s.inFluid) {
                     integrator_.Integrate(
                         [&c, layer](double r, const Vector<E::FluidSize>& m) {
                           return E::FluidMinors(E::Fluid(c, c(layer, r)), m);
                         },
                         r0, r1, s.fluid, step);
                   } else {
                     integrator_.Integrate(
                         [&c, layer](double r, const Vector<E::SolidSize>& m) {
                           return E::SolidMinors(E::Solid(c, c(layer, r)), m);
                         },
                         r0, r1, s.solid, step);
                   }
                 });
    };

    // Upwards from the starting level to the lower of source and receiver.
    const auto lower = Lower();
    const auto upper = Upper();
    const auto start = StartPoint(l, omega.real(), lower, 0);
    auto h = 0.0;
    const auto point = StartRadius(start, omega, h);
    auto s1 = State{};
    {
      const auto m = c(point.layer, point.r);
      s1.inFluid = model_->GetLayer(point.layer).IsFluid();
      if (s1.inFluid) {
        s1.fluid.value = E::FluidStart(c, m);
        s1.fluid.Normalise();
      } else {
        s1.solid.value = E::SolidStart(c, m);
        s1.solid.Normalise();
      }
    }
    walkMinors(point, lower, s1, h);
    const auto ms1 = s1.solid;

    // Continued to the upper of the two.
    auto s3 = s1;
    walkMinors(lower, upper, s3, h);
    const auto& ms3 = s3.solid;

    // Downwards from the surface, where the traction vanishes.
    auto s2 = State{};
    {
      const auto top = model_->NumberOfLayers() - 1;
      s2.inFluid = model_->GetLayer(top).IsFluid();
      s2.fluid.value[0] = 1;
      s2.solid.value[0] = 1;
      auto hDown = 0.0;
      walkMinors(Point{top, model_->SurfaceRadius()}, upper, s2, hDown);
    }
    const auto& ms2 = s2.solid;

    const auto delta = ScaledComplex{SecularDeterminant(ms3.value, ms2.value),
                                     ms3.logScale + ms2.logScale};

    // The b-vectors for the four jumps, carried together from the source to
    // the receiver, since their coefficients are the same.
    auto b = std::array<ScaledVector<E::BSize>, 4>{};
    for (auto j : {0, 1, 2, 3}) {
      if (ReceiverAbove()) {
        b[j].value = E::B(ms1.value, E::Jump(j));
        b[j].logScale = ms1.logScale;
      } else {
        b[j].value = E::B(ms2.value, E::Jump(j));
        const auto scale = ScaledComplex{-1.0, ms2.logScale} / delta;
        for (auto& v : b[j].value) v *= scale.value;
        b[j].logScale = scale.logScale;
      }
      b[j].Normalise();
    }
    auto hb = h;
    WalkLayers(source_, receiver_, hb, omega,
               [&](Int layer, double r0, double r1, double& step) {
                 integrator_.IntegrateTogether(
                     [&c, layer](double r, const Vectors<E::BSize, 4>& v) {
                       const auto block = E::Solid(c, c(layer, r));
                       auto d = Vectors<E::BSize, 4>{};
                       for (auto k = std::size_t{0}; k < 4; ++k) {
                         d[k] = BDerivative(block, v[k]);
                       }
                       return d;
                     },
                     r0, r1, b, step);
               });

    auto result = SpheroidalResponse{};
    for (auto j : {0, 1, 2, 3}) {
      auto y = std::array<Complex, 4>{};
      auto logScale = 0.0;
      if (ReceiverAbove()) {
        y = E::Components(Solution(ms2.value, b[j].value));
        const auto scale =
            ScaledComplex{E::Sign, ms2.logScale + b[j].logScale} / delta;
        for (auto& v : y) v *= scale.value;
        logScale = scale.logScale;
      } else {
        y = E::Components(Solution(ms1.value, b[j].value));
        logScale = ms1.logScale + b[j].logScale;
      }
      for (auto i : {0, 1, 2, 3}) result[i][j] = Rescale(y[i], logScale);
    }
    return result;
  }
};

}  // namespace YSpec
