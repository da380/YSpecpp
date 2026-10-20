#include <YSpec/RadialSolver.hpp>
#include <YSpec/Synthetics.hpp>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <format>
#include <map>
#include <ranges>
#include <stdexcept>
#include <vector>

#include "AngularFunctions.hpp"
#include "Equations.hpp"
#include "Geometry.hpp"
#include "OrderArray.hpp"
#include "SourceWeights.hpp"
#include "TimeSeries.hpp"

namespace YSpec {

namespace {

using namespace std::complex_literals;

// The members of a set that share a depth, and the radius of that depth.
struct DepthGroup {
  double radius;
  std::vector<std::size_t> members;
};

template <class Points>
std::vector<DepthGroup> GroupByDepth(const Points& points, const Model& model,
                                     const char* what) {
  const auto surface = model.SurfaceRadius();
  const auto length = model.GetScales().Length();
  auto groups = std::map<double, std::vector<std::size_t>>{};
  for (auto i : std::views::iota(std::size_t{0}, std::size(points))) {
    groups[points[i].depth].push_back(i);
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

// The responses at one degree and frequency for every source depth and
// receiver depth, converted to the quantity recorded and corrected, if asked,
// for the potential and tilt.
struct Responses {
  ResponseTable<ToroidalResponse> toroidal{0, 0};
  ResponseTable<SpheroidalResponse> spheroidal{0, 0};
  ResponseTable<RadialResponse> radial{0, 0};
};

// What the corrections need at a receiver depth.
struct ReceiverGravity {
  double r;   // The radius.
  double g;   // Gravity.
  double dg;  // Its radial derivative.
};

Responses Record(const RadialSolver& solver, int l, Complex omega,
                 const Calculation& calculation,
                 std::span<const ReceiverGravity> gravity) {
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

// The displacement coefficients for one source, order by order.
struct Coefficients {
  OrderArray<Complex> u;  // radial
  OrderArray<Complex> v;  // consoidal
  OrderArray<Complex> w;  // toroidal
};

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

// Where each receiver is seen from each source, and the angular functions
// there. These depend only on where the receiver is on the surface, so
// receivers at one place and several depths, as in a snapshot, share them.
class Views {
 public:
  Views(std::span<const Source> sources, std::span<const Receiver> receivers,
        int maxDegree)
      : place_(receivers.size()) {
    auto places = std::vector<Location>{};
    for (auto k : std::views::iota(std::size_t{0}, receivers.size())) {
      const auto& location = receivers[k].location;
      const auto found = std::ranges::find_if(places, [&](const Location& p) {
        return p.latitude == location.latitude &&
               p.longitude == location.longitude;
      });
      place_[k] =
          static_cast<std::size_t>(std::distance(places.begin(), found));
      if (found == places.end()) places.push_back(location);
    }
    places_ = places.size();
    auto distances = std::vector<double>{};
    for (const auto& source : sources) {
      for (const auto& place : places) {
        auto view =
            View{.geometry = Geometry(place, source.location), .phase = {}};
        for (auto m : OrderArray<Complex>::Orders(2)) {
          view.phase[m] = std::exp(1i * (m * view.geometry.phi));
        }
        views_.push_back(view);
        distances.push_back(view.geometry.distance);
      }
    }
    angular_ = AngularFunctions(maxDegree, distances);
  }

  // A receiver as seen from a source.
  struct View {
    ReceiverGeometry geometry;
    OrderArray<Complex> phase;  // exp(i m phi)
  };

  const View& Of(std::size_t source, std::size_t receiver) const {
    return views_[Index(source, receiver)];
  }
  const Angular& Functions(std::size_t source, std::size_t receiver, int l,
                           int m) const {
    return angular_(static_cast<Int>(Index(source, receiver)), l, m);
  }

 private:
  std::vector<std::size_t> place_;  // the place of each receiver
  std::size_t places_ = 0;
  std::vector<View> views_;
  AngularFunctions angular_{0, {}};

  std::size_t Index(std::size_t source, std::size_t receiver) const {
    return source * places_ + place_[receiver];
  }
};

// The displacement at a receiver, in the frame with the source at the pole,
// from the coefficients at one degree.
std::array<Complex, 3> Sum(const Coefficients& coefficients, int l,
                           const Views& views, std::size_t source,
                           std::size_t receiver) {
  const auto& phase = views.Of(source, receiver).phase;
  auto uR = Complex{};
  auto uTheta = Complex{};
  auto uPhi = Complex{};
  for (auto m : OrderArray<Complex>::Orders(l)) {
    const auto& f = views.Functions(source, receiver, l, m);
    const auto e = phase[m];
    const auto u = coefficients.u[m];
    const auto v = coefficients.v[m];
    const auto w = coefficients.w[m];
    uR += u * f.value * e;
    uTheta += (v * f.dTheta + 1i * w * f.mOverSin) * e;
    uPhi += (1i * v * f.mOverSin - w * f.dTheta) * e;
  }
  return {uR, uTheta, uPhi};
}

// The spectra in the frame with the source at the pole, for every source,
// receiver and component (r, theta, phi), at the grid's frequencies.
class Spectra {
 public:
  Spectra(std::size_t sources, std::size_t receivers, std::size_t frequencies)
      : receivers_{receivers},
        frequencies_{frequencies},
        values_(sources * receivers * 3 * frequencies) {}

  Complex& operator()(std::size_t source, std::size_t receiver, int component,
                      std::size_t frequency) {
    return values_[((source * receivers_ + receiver) * 3 + component) *
                       frequencies_ +
                   frequency];
  }
  const std::vector<Complex>& Values() const { return values_; }

 private:
  std::size_t receivers_;
  std::size_t frequencies_;
  std::vector<Complex> values_;
};

}  // namespace

Seismograms ComputeSeismograms(const Model& model,
                               const Calculation& calculation,
                               std::span<const Source> sources,
                               std::span<const Receiver> receivers,
                               const Progress& progress) {
  const auto& scales = model.GetScales();
  const auto& degrees = calculation.degrees;
  if (degrees.min < 0 || degrees.max < degrees.min) {
    throw std::invalid_argument("the degrees must satisfy 0 <= min <= max");
  }
  if (sources.empty() || receivers.empty()) {
    throw std::invalid_argument("there must be a source and a receiver");
  }

  // The frequency grid and the taper, non-dimensional.
  const auto time = scales.Time();
  const auto& sampling = calculation.sampling;
  const auto grid =
      FrequencyGrid(calculation.band.min * time, calculation.band.max * time,
                    sampling.timeStep / time, sampling.duration / time);
  auto taper = calculation.taper;
  if (taper) {
    for (auto* f : {&taper->f1, &taper->f2, &taper->f3, &taper->f4}) {
      *f *= time;
    }
  }

  // One solver for every source depth and every receiver depth.
  const auto sourceDepths = GroupByDepth(sources, model, "source");
  const auto receiverDepths = GroupByDepth(receivers, model, "receiver");
  auto radii = [](const std::vector<DepthGroup>& groups) {
    auto r = std::vector<double>{};
    for (const auto& group : groups) r.push_back(group.radius);
    return r;
  };
  const auto solver = RadialSolver(model, radii(sourceDepths),
                                   radii(receiverDepths), calculation.gravity,
                                   calculation.attenuation, calculation.solver);
  auto gravity = std::vector<ReceiverGravity>{};
  for (auto i : std::views::iota(std::size_t{0}, receiverDepths.size())) {
    const auto receiver = solver.Receiver(i);
    const auto at = model(receiver.layer, receiver.r);
    gravity.push_back({.r = receiver.r,
                       .g = at.g,
                       .dg = 4 * PiG * at.rho - 2 * at.g / receiver.r});
  }

  const auto views = Views(sources, receivers, degrees.max);
  const auto receiverCount = receivers.size();
  const auto frequencyCount = grid.Count();
  auto spectra = Spectra(sources.size(), receiverCount, frequencyCount);
  auto done = std::atomic<std::size_t>{0};
  auto error = std::exception_ptr{};

#pragma omp parallel for schedule(dynamic)
  for (long index = 0; index < static_cast<long>(frequencyCount); ++index) {
    try {
      const auto frequency = static_cast<std::size_t>(index);
      const auto omega = grid.Omega(frequency);
      // The sums at this frequency, gathered locally so that threads at
      // neighbouring frequencies do not write to the same cache lines.
      auto sums =
          std::vector<std::array<Complex, 3>>(sources.size() * receiverCount);
      for (auto l = degrees.min; l <= degrees.max; ++l) {
        const auto responses = Record(solver, l, omega, calculation, gravity);
        const auto c = Context(model, l, omega, calculation.attenuation);
        for (auto sourceDepth :
             std::views::iota(std::size_t{0}, sourceDepths.size())) {
          const auto point = solver.Source(sourceDepth);
          const auto atSource = c(point.layer, point.r);
          for (auto source : sourceDepths[sourceDepth].members) {
            const auto moment = sources[source].moment / scales.Moment();
            for (auto receiverDepth :
                 std::views::iota(std::size_t{0}, receiverDepths.size())) {
              const auto coefficients = Weight(
                  responses, sourceDepth, receiverDepth, c, atSource, moment);
              for (auto receiver : receiverDepths[receiverDepth].members) {
                const auto u = Sum(coefficients, l, views, source, receiver);
                auto& sum = sums[source * receiverCount + receiver];
                for (auto i : {0, 1, 2}) sum[i] += u[i];
              }
            }
          }
        }
      }
      for (auto source : std::views::iota(std::size_t{0}, sources.size())) {
        for (auto receiver : std::views::iota(std::size_t{0}, receiverCount)) {
          const auto& sum = sums[source * receiverCount + receiver];
          for (auto i : {0, 1, 2}) {
            spectra(source, receiver, i, frequency) = sum[i];
          }
        }
      }
      const auto finished = ++done;
      if (progress) {
#pragma omp critical(YSpecProgress)
        progress(finished, frequencyCount);
      }
    } catch (...) {
#pragma omp critical(YSpecError)
      if (!error) error = std::current_exception();
    }
  }
  if (error) std::rethrow_exception(error);

  // To time series, rotated to north and east, in SI units.
  const auto series = InverseTransform(grid, taper, spectra.Values(),
                                       sources.size() * receiverCount * 3);
  const auto samples = std::min<std::size_t>(
      grid.Samples(),
      static_cast<std::size_t>(sampling.duration / sampling.timeStep + 1e-9) +
          1);
  auto unit = scales.Length();
  if (calculation.quantity == Quantity::Velocity) unit = scales.Velocity();
  if (calculation.quantity == Quantity::Acceleration) {
    unit = scales.Acceleration();
  }
  auto seismograms =
      Seismograms(sources.size(), receiverCount, samples, sampling.timeStep);
  for (auto source : std::views::iota(std::size_t{0}, sources.size())) {
    for (auto receiver : std::views::iota(std::size_t{0}, receiverCount)) {
      const auto rotation = views.Of(source, receiver).geometry.rotation;
      const auto cosine = std::cos(rotation);
      const auto sine = std::sin(rotation);
      auto component = [&](int i) {
        const auto first =
            ((source * receiverCount + receiver) * 3 + i) * grid.Samples();
        return std::span{series}.subspan(first, samples);
      };
      const auto uR = component(0);
      const auto uTheta = component(1);
      const auto uPhi = component(2);
      for (auto i : std::views::iota(std::size_t{0}, samples)) {
        seismograms(source, receiver, Vertical, i) = uR[i] * unit;
        seismograms(source, receiver, North, i) =
            (-uTheta[i] * cosine + uPhi[i] * sine) * unit;
        seismograms(source, receiver, East, i) =
            (uTheta[i] * sine + uPhi[i] * cosine) * unit;
      }
    }
  }
  return seismograms;
}

}  // namespace YSpec
