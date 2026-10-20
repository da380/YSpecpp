#pragma once

/**
 * @file Synthetics.hpp
 * @brief Synthetic seismograms: the radial solutions for every degree and
 * frequency, combined with the source, summed over degree and order at the
 * receivers, and transformed to time series.
 *
 * @details The work is a loop over frequency, threaded, with a loop over
 * degree inside it: each thread adds into the spectra at its own frequency,
 * so nothing is shared but read-only data. The receivers all lie at one
 * depth, and the source is rotated to the north pole, where it excites only
 * the orders @f$|m| \le 2@f$.
 *
 * The output is displacement, velocity or acceleration at each receiver, in
 * SI units, in the vertical (up), north and east directions. For acceleration
 * the change in gravity seen by a moving instrument, and the tilt of its
 * axis, can be included ("potential and tilt corrections"); for degree zero
 * the correction reduces to the free-air term @f$-2gU/r@f$.
 */

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstddef>
#include <format>
#include <functional>
#include <numbers>
#include <ranges>
#include <stdexcept>
#include <vector>

#include "Equations.hpp"
#include "Model.hpp"
#include "Receivers.hpp"
#include "Solver.hpp"
#include "Source.hpp"
#include "TimeSeries.hpp"

namespace YSpec {

/** @brief The quantity the seismograms record. */
enum class Quantity { Displacement, Velocity, Acceleration };

/**
 * @brief Everything that defines a calculation, in the units of the input
 * file: km, degrees, mHz, minutes, seconds and N m.
 */
struct Settings {
  Gravity gravity = Gravity::Full;  ///< How gravity enters.
  bool attenuation = true;          ///< Whether to include attenuation.
  Quantity quantity = Quantity::Displacement;  ///< What is recorded.
  /// Whether to correct acceleration for the potential and for tilt.
  bool corrections = false;
  int lMin = 0;         ///< The lowest degree.
  int lMax = 0;         ///< The highest degree.
  double fMin = 0;      ///< The lowest frequency (mHz).
  double fMax = 0;      ///< The highest frequency (mHz).
  double length = 0;    ///< The length of the seismograms (minutes).
  double timeStep = 0;  ///< The time step (s).
  /// The corners of the cosine taper (mHz); all zero for none.
  std::array<double, 4> filter = {};
  double sourceDepth = 0;           ///< The depth of the source (km).
  Location source = {};             ///< The epicentre.
  MomentTensor moment = {};         ///< The moment tensor (N m).
  double receiverDepth = 0;         ///< The depth of the receivers (km).
  std::vector<Location> receivers;  ///< The receivers.
  SolverOptions solver = {};        ///< Settings of the radial integration.
};

/**
 * @brief Seismograms at a set of receivers.
 */
class Seismograms {
 public:
  /** @brief Zero seismograms of a given size. */
  Seismograms(std::size_t receivers, std::size_t samples, double timeStep)
      : receivers_{receivers},
        samples_{samples},
        timeStep_{timeStep},
        data_(receivers * 3 * samples) {}

  /** @brief The number of receivers. */
  std::size_t Receivers() const { return receivers_; }
  /** @brief The number of samples. */
  std::size_t Samples() const { return samples_; }
  /** @brief The time step (s). */
  double TimeStep() const { return timeStep_; }

  /**
   * @brief A sample.
   * @param receiver The receiver.
   * @param component 0 for vertical (up), 1 for north, 2 for east.
   * @param i The sample, at time @f$i\,\Delta t@f$ after the origin time.
   */
  double& operator()(std::size_t receiver, int component, std::size_t i) {
    return data_[(receiver * 3 + component) * samples_ + i];
  }
  /** @brief A sample, read-only. */
  double operator()(std::size_t receiver, int component, std::size_t i) const {
    return data_[(receiver * 3 + component) * samples_ + i];
  }

 private:
  std::size_t receivers_;
  std::size_t samples_;
  double timeStep_;
  std::vector<double> data_;
};

/**
 * @brief Called as frequencies are completed, with the number done and the
 * number in all; called from one thread at a time.
 */
using Progress = std::function<void(std::size_t, std::size_t)>;

namespace SyntheticsDetails {

// Applies the time derivatives of the quantity recorded.
inline Complex Derivative(Quantity quantity, Complex omega, Complex y) {
  switch (quantity) {
    case Quantity::Displacement:
      return y;
    case Quantity::Velocity:
      return Complex(0, 1) * omega * y;
    case Quantity::Acceleration:
      return -omega * omega * y;
  }
  return y;
}

}  // namespace SyntheticsDetails

/**
 * @brief Computes the spectra at the receivers.
 * @param model The earth model.
 * @param settings The calculation.
 * @param grid The frequency grid.
 * @param progress Called as frequencies are completed, if given.
 * @return For each receiver and each of the components r, theta and phi in the
 * frame with the source at the pole, the spectrum at the grid's frequencies,
 * non-dimensional.
 */
inline std::vector<Complex> ComputeSpectra(const Model& model,
                                           const Settings& settings,
                                           const FrequencyGrid& grid,
                                           const Progress& progress = {}) {
  using namespace std::complex_literals;
  const auto& scales = model.GetScales();
  const auto surface = model.SurfaceRadius();
  const auto rs = surface - settings.sourceDepth * 1e3 / scales.Length();
  const auto rr = surface - settings.receiverDepth * 1e3 / scales.Length();
  const auto moment = settings.moment / scales.Moment();
  const auto lMin = settings.lMin;
  const auto lMax = settings.lMax;
  if (lMin < 0 || lMax < lMin) {
    throw std::invalid_argument("the degrees must satisfy 0 <= min <= max");
  }

  auto inside = [&](double r, double depth, const char* what) {
    if (r < model.InnerRadius() || r > surface) {
      throw std::invalid_argument(std::format(
          "the {} depth of {} km is outside the model, whose radius is {} km",
          what, depth, surface * scales.Length() / 1e3));
    }
  };
  inside(rs, settings.sourceDepth, "source");
  inside(rr, settings.receiverDepth, "receiver");
  const auto solver = RadialSolver<>(model, rs, rr, settings.gravity,
                                     settings.attenuation, settings.solver);
  const auto source = solver.Source();
  const auto receiver = solver.Receiver();
  const auto g = model(receiver.layer, rr).g;
  const auto quantity = settings.quantity;
  const auto corrections =
      settings.corrections && quantity == Quantity::Acceleration;

  // Each receiver's place in the frame with the source at the pole.
  const auto nr = settings.receivers.size();
  auto distance = std::vector<double>(nr);
  auto phase = std::vector<std::array<Complex, 5>>(nr);
  for (auto k : std::views::iota(std::size_t{0}, nr)) {
    const auto geometry = Geometry(settings.receivers[k], settings.source);
    distance[k] = geometry.distance;
    for (auto m : std::views::iota(-2, 3)) {
      phase[k][m + 2] = std::exp(1i * (m * geometry.phi));
    }
  }
  const auto angular = AngularFunctions(lMax, distance);

  const auto nf = grid.Count();
  auto spectra = std::vector<Complex>(nr * 3 * nf);

  auto done = std::atomic<std::size_t>{0};
  auto error = std::exception_ptr{};
  const auto count = static_cast<long>(nf);

#pragma omp parallel for schedule(dynamic)
  for (long jj = 0; jj < count; ++jj) {
    try {
      const auto j = static_cast<std::size_t>(jj);
      const auto omega = grid.Omega(j);
      auto derivative = [&](Complex y) {
        return SyntheticsDetails::Derivative(quantity, omega, y);
      };
      // The sums at this frequency, gathered locally so that threads at
      // neighbouring frequencies do not write to the same cache lines.
      auto local = std::vector<Complex>(nr * 3);
      auto at = [&local](std::size_t k, int component) -> Complex& {
        return local[k * 3 + component];
      };
      for (auto l = lMin; l <= lMax; ++l) {
        const auto c = Context(model, l, omega, settings.attenuation);
        const auto atSource = c(source.layer, source.r);
        if (l == 0) {
          const auto y = solver.Radial(omega);
          const auto w = Radial(c, atSource, moment);
          const auto u0 = (w[0] * y[0] + w[1] * y[1]) / rr;
          auto u = derivative(u0);
          if (corrections) u -= 2 * g * u0 / rr;
          for (auto k : std::views::iota(std::size_t{0}, nr)) {
            at(k, 0) += u * angular(k, 0, 0).x;
          }
          continue;
        }
        const auto zeta = c.Zeta();
        const auto mMax = std::min(l, 2);

        // Toroidal.
        {
          const auto y = solver.Toroidal(l, omega);
          const auto w = Toroidal(c, atSource, moment);
          auto coefficient = std::array<Complex, 5>{};
          for (auto m : std::views::iota(-mMax, mMax + 1)) {
            const auto& wm = w[m + 2];
            coefficient[m + 2] =
                derivative((wm[0] * y[0] + wm[1] * y[1]) / (zeta * rr));
          }
          for (auto k : std::views::iota(std::size_t{0}, nr)) {
            auto ut = Complex{};
            auto up = Complex{};
            for (auto m : std::views::iota(-mMax, mMax + 1)) {
              const auto& a = angular(k, l, m);
              const auto term = coefficient[m + 2] * phase[k][m + 2];
              ut += 1i * a.mx * term;
              up -= a.dx * term;
            }
            at(k, 1) += ut;
            at(k, 2) += up;
          }
        }

        // Spheroidal.
        {
          auto y = solver.Spheroidal(l, omega);
          for (auto i : {0, 1, 2, 3}) {
            const auto U = y[0][i];
            const auto V = y[1][i];
            const auto P = y[2][i];
            const auto Psi = y[3][i];
            y[0][i] = derivative(U);
            y[1][i] = derivative(V);
            if (corrections) {
              // The variables are r U, zeta r V, r phi and the flux
              // psi = r (phi' + (l+1) phi / r + 4 pi G rho U) / (4 pi G). The
              // radial correction is r (phi' + U g'), with g' = 4 pi G rho
              // - 2 g / r; the psi term vanishes at a free surface, where it
              // reduces to the original YSpec's.
              y[0][i] += 4 * PiG * Psi - ((l + 1.0) * P + 2 * g * U) / rr;
              y[1][i] += zeta * (P + g * U) / rr;
            }
          }
          const auto w = Spheroidal(c, atSource, moment);
          auto u = std::array<Complex, 5>{};
          auto v = std::array<Complex, 5>{};
          for (auto m : std::views::iota(-mMax, mMax + 1)) {
            for (auto i : {0, 1, 2, 3}) {
              u[m + 2] += y[0][i] * w[m + 2][i];
              v[m + 2] += y[1][i] * w[m + 2][i];
            }
            u[m + 2] /= rr;
            v[m + 2] /= zeta * rr;
          }
          for (auto k : std::views::iota(std::size_t{0}, nr)) {
            auto ur = Complex{};
            auto ut = Complex{};
            auto up = Complex{};
            for (auto m : std::views::iota(-mMax, mMax + 1)) {
              const auto& a = angular(k, l, m);
              const auto e = phase[k][m + 2];
              ur += u[m + 2] * a.x * e;
              ut += v[m + 2] * a.dx * e;
              up += 1i * v[m + 2] * a.mx * e;
            }
            at(k, 0) += ur;
            at(k, 1) += ut;
            at(k, 2) += up;
          }
        }
      }
      for (auto k : std::views::iota(std::size_t{0}, nr)) {
        for (auto c : {0, 1, 2}) spectra[(k * 3 + c) * nf + j] = at(k, c);
      }
      const auto finished = ++done;
      if (progress) {
#pragma omp critical(YSpecProgress)
        progress(finished, nf);
      }
    } catch (...) {
#pragma omp critical(YSpecError)
      if (!error) error = std::current_exception();
    }
  }
  if (error) std::rethrow_exception(error);
  return spectra;
}

/**
 * @brief Computes synthetic seismograms.
 * @param model The earth model.
 * @param settings The calculation.
 * @param progress Called as frequencies are completed, if given.
 */
inline Seismograms ComputeSeismograms(const Model& model,
                                      const Settings& settings,
                                      const Progress& progress = {}) {
  const auto& scales = model.GetScales();
  const auto time = scales.Time();
  // Frequencies in mHz to cycles per unit time.
  auto frequency = [time](double mHz) { return mHz / 1000 * time; };
  const auto grid =
      FrequencyGrid(frequency(settings.fMin), frequency(settings.fMax),
                    settings.timeStep / time, settings.length * 60 / time);
  const auto filter = CosineFilter{
      frequency(settings.filter[0]), frequency(settings.filter[1]),
      frequency(settings.filter[2]), frequency(settings.filter[3])};

  const auto spectra = ComputeSpectra(model, settings, grid, progress);
  const auto nr = settings.receivers.size();
  const auto series = InverseTransform(grid, filter, spectra, nr * 3);

  // The samples up to the length asked for.
  const auto samples = std::min<std::size_t>(
      grid.Samples(), static_cast<std::size_t>(
                          settings.length * 60 / settings.timeStep + 1e-9) +
                          1);
  auto unit = scales.Length();
  if (settings.quantity == Quantity::Velocity) unit = scales.Velocity();
  if (settings.quantity == Quantity::Acceleration) unit = scales.Acceleration();

  auto seismograms = Seismograms(nr, samples, settings.timeStep);
  const auto nt = grid.Samples();
  for (auto k : std::views::iota(std::size_t{0}, nr)) {
    const auto rotation =
        Geometry(settings.receivers[k], settings.source).rotation;
    const auto c = std::cos(rotation);
    const auto s = std::sin(rotation);
    for (auto i : std::views::iota(std::size_t{0}, samples)) {
      const auto ur = series[(k * 3 + 0) * nt + i] * unit;
      const auto ut = series[(k * 3 + 1) * nt + i] * unit;
      const auto up = series[(k * 3 + 2) * nt + i] * unit;
      seismograms(k, 0, i) = ur;
      seismograms(k, 1, i) = -ut * c + up * s;
      seismograms(k, 2, i) = ut * s + up * c;
    }
  }
  return seismograms;
}

}  // namespace YSpec
