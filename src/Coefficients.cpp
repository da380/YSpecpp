#include "Coefficients.hpp"

#include <format>
#include <map>
#include <ranges>
#include <stdexcept>

#include "SourceWeights.hpp"

namespace YSpec {

using namespace std::complex_literals;

FrequencyGrid MakeGrid(const Calculation& calculation, const Scales& scales) {
  const auto time = scales.Time();
  const auto& sampling = calculation.sampling;
  return FrequencyGrid(calculation.band.min * time, calculation.band.max * time,
                       sampling.timeStep / time, sampling.duration / time);
}

std::optional<Taper> MakeTaper(const Calculation& calculation,
                               const Scales& scales) {
  auto taper = calculation.taper;
  if (taper) {
    for (auto* f : {&taper->f1, &taper->f2, &taper->f3, &taper->f4}) {
      *f *= scales.Time();
    }
  }
  return taper;
}

// The SI unit of a quantity, in the model's scales.
double Unit(const Scales& scales, Quantity quantity) {
  switch (quantity) {
    case Quantity::Displacement:
      return scales.Length();
    case Quantity::Velocity:
      return scales.Velocity();
    case Quantity::Acceleration:
      return scales.Acceleration();
  }
  return scales.Length();
}

// The number of time derivatives of displacement a quantity is.
int Derivatives(Quantity quantity) {
  switch (quantity) {
    case Quantity::Displacement:
      return 0;
    case Quantity::Velocity:
      return 1;
    case Quantity::Acceleration:
      return 2;
  }
  return 0;
}

namespace {

// Groups depths (m), as radii.
std::vector<DepthGroup> GroupByDepth(std::span<const double> depths,
                                     const Model& model, const char* what) {
  const auto surface = model.SurfaceRadius();
  const auto length = model.GetScales().Length();
  auto groups = std::map<double, std::vector<std::size_t>>{};
  for (auto i : std::views::iota(std::size_t{0}, depths.size())) {
    groups[depths[i]].push_back(i);
  }
  auto result = std::vector<DepthGroup>{};
  for (auto& [depth, members] : groups) {
    const auto r = surface - depth / length;
    if (depth < 0 || r < model.InnerRadius()) {
      throw std::invalid_argument(std::format(
          "a {} depth of {} km is outside the model, whose radius is {} km",
          what, depth / 1e3, surface * length / 1e3));
    }
    result.push_back({r, std::move(members)});
  }
  return result;
}

std::vector<double> Radii(const std::vector<DepthGroup>& groups) {
  auto radii = std::vector<double>{};
  for (const auto& group : groups) radii.push_back(group.radius);
  return radii;
}

std::vector<double> SourceDepths(std::span<const Source> sources) {
  auto depths = std::vector<double>{};
  for (const auto& source : sources) depths.push_back(source.depth);
  return depths;
}

}  // namespace

DepthSolves::DepthSolves(const Model& model, const Calculation& calculation,
                         std::span<const Source> sources,
                         std::span<const double> depths)
    : sourceDepths{GroupByDepth(SourceDepths(sources), model, "source")},
      receiverDepths{GroupByDepth(depths, model, "receiver")},
      solver{model,
             Radii(sourceDepths),
             Radii(receiverDepths),
             calculation.gravity,
             calculation.attenuation,
             calculation.solver} {
  for (auto i : std::views::iota(std::size_t{0}, solver.Receivers())) {
    const auto receiver = solver.Receiver(i);
    const auto at = model(receiver.layer, receiver.r);
    gravity.push_back({.r = receiver.r,
                       .g = at.g,
                       .dg = 4 * PiG * at.rho - 2 * at.g / receiver.r});
  }
}

Responses Record(const DepthSolves& solves, int l, Complex omega,
                 const Calculation& calculation) {
  const auto& solver = solves.solver;
  const auto& gravity = solves.gravity;
  const auto derivative = [&](Complex y) {
    switch (calculation.quantity) {
      case Quantity::Displacement:
        return y;
      case Quantity::Velocity:
        return 1i * omega * y;
      case Quantity::Acceleration:
        return -omega * omega * y;
    }
    return y;
  };
  const auto corrections =
      calculation.corrections && calculation.quantity == Quantity::Acceleration;
  const auto full = calculation.gravity == Gravity::Full;
  auto result = Responses{};
  if (l == 0) {
    result.radial = solver.Radial(omega);
  } else {
    result.toroidal = solver.Toroidal(l, omega);
    result.spheroidal = solver.Spheroidal(l, omega);
  }
  for (auto source : std::views::iota(std::size_t{0}, solver.Sources())) {
    for (auto receiver : std::views::iota(std::size_t{0}, solver.Receivers())) {
      const auto& at = gravity[receiver];
      if (l == 0) {
        for (auto& u : result.radial(source, receiver).U) {
          // For l = 0 the change in gravity reduces to the free-air term.
          u = derivative(u) - (corrections ? 2 * at.g * u / at.r : 0.0);
        }
        continue;
      }
      for (auto& w : result.toroidal(source, receiver).W) w = derivative(w);
      auto& spheroidal = result.spheroidal(source, receiver);
      for (auto j : {0, 1, 2, 3}) {
        const auto U = spheroidal.U[j];
        const auto V = spheroidal.V[j];
        spheroidal.U[j] = derivative(U);
        spheroidal.V[j] = derivative(V);
        if (!corrections) continue;
        if (full) {
          // The change in gravity at the moving instrument, phi' + U g',
          // and the tilt of its axis, (phi + g U) / r.
          spheroidal.U[j] += spheroidal.DPhi[j] + U * at.dg;
          spheroidal.V[j] += (spheroidal.Phi[j] + at.g * U) / at.r;
        } else {
          // Without the potential, the free-air term and the tilt.
          spheroidal.U[j] -= 2 * at.g * U / at.r;
          spheroidal.V[j] += at.g * U / at.r;
        }
      }
    }
  }
  return result;
}

Coefficients Weight(const Responses& responses, std::size_t sourceDepth,
                    std::size_t receiverDepth, const Context& c,
                    const ComplexMaterial& atSource,
                    const MomentTensor& moment) {
  auto result = Coefficients{};
  const auto l = c.Degree();
  if (l == 0) {
    const auto& u = responses.radial(sourceDepth, receiverDepth).U;
    const auto weights = Radial(c, atSource, moment);
    result.u[0] = weights[0] * u[0] + weights[1] * u[1];
    return result;
  }
  const auto& toroidal = responses.toroidal(sourceDepth, receiverDepth);
  const auto& spheroidal = responses.spheroidal(sourceDepth, receiverDepth);
  const auto toroidalWeights = Toroidal(c, atSource, moment);
  const auto spheroidalWeights = Spheroidal(c, atSource, moment);
  for (auto m : OrderArray<Complex>::Orders(l)) {
    for (auto j : {0, 1}) {
      result.w[m] += toroidalWeights[m][j] * toroidal.W[j];
    }
    for (auto j : {0, 1, 2, 3}) {
      result.u[m] += spheroidalWeights[m][j] * spheroidal.U[j];
      result.v[m] += spheroidalWeights[m][j] * spheroidal.V[j];
    }
  }
  return result;
}

}  // namespace YSpec
