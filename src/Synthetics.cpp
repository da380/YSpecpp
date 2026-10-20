#include <YSpec/Synthetics.hpp>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <ranges>
#include <stdexcept>
#include <vector>

#include "AngularFunctions.hpp"
#include "Coefficients.hpp"
#include "Geometry.hpp"
#include "OrderArray.hpp"
#include "TimeSeries.hpp"

namespace YSpec {

namespace {

using namespace std::complex_literals;

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
    uTheta += (v * f.dTheta + TimesI(w) * f.mOverSin) * e;
    uPhi += (TimesI(v) * f.mOverSin - w * f.dTheta) * e;
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

// The response of an instrument to the quantity computed, at the grid's
// frequencies (non-dimensional, with the unit of time given): its response to
// the motion it responds to, converted by a power of i omega.
std::vector<Complex> InstrumentResponse(const Instrument& instrument,
                                        const FrequencyGrid& grid, double time,
                                        Quantity computed) {
  const auto power = Derivatives(instrument.input) - Derivatives(computed);
  auto response = std::vector<Complex>{};
  for (auto j : std::views::iota(std::size_t{0}, grid.Count())) {
    const auto omega = grid.Omega(j) / time;
    response.push_back(instrument.Response(omega) *
                       std::pow(1i * omega, power));
  }
  return response;
}

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

  const auto grid = MakeGrid(calculation, scales);
  const auto taper = MakeTaper(calculation, scales);
  const auto& sampling = calculation.sampling;

  // One solver for every source depth and every receiver depth.
  auto depths = std::vector<double>{};
  for (const auto& receiver : receivers) depths.push_back(receiver.depth);
  const auto solves = DepthSolves(model, calculation, sources, depths);

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
        ForEachSource(solves, model, calculation, sources, l, omega,
                      [&](std::size_t source, std::size_t receiverDepth,
                          const Coefficients& coefficients) {
                        for (auto receiver :
                             solves.receiverDepths[receiverDepth].members) {
                          const auto u =
                              Sum(coefficients, l, views, source, receiver);
                          auto& sum = sums[source * receiverCount + receiver];
                          for (auto i : {0, 1, 2}) sum[i] += u[i];
                        }
                      });
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

  // Rotated to vertical, north and east, and through the instruments, if
  // any: both are linear and the same at every time, so they act on the
  // spectra. The instruments are evaluated at the same complex frequencies,
  // which is exact: the damping passes through a convolution unchanged.
  for (auto receiver : std::views::iota(std::size_t{0}, receiverCount)) {
    const auto& instrument = receivers[receiver].instrument;
    const auto response =
        instrument ? InstrumentResponse(*instrument, grid, scales.Time(),
                                        calculation.quantity)
                   : std::vector<Complex>(frequencyCount, 1.0);
    for (auto source : std::views::iota(std::size_t{0}, sources.size())) {
      const auto rotation = views.Of(source, receiver).geometry.rotation;
      const auto cosine = std::cos(rotation);
      const auto sine = std::sin(rotation);
      for (auto frequency : std::views::iota(std::size_t{0}, frequencyCount)) {
        auto& uR = spectra(source, receiver, 0, frequency);
        auto& uTheta = spectra(source, receiver, 1, frequency);
        auto& uPhi = spectra(source, receiver, 2, frequency);
        const auto vertical = uR;
        const auto north = -uTheta * cosine + uPhi * sine;
        const auto east = uTheta * sine + uPhi * cosine;
        uR = vertical * response[frequency];
        uTheta = north * response[frequency];
        uPhi = east * response[frequency];
      }
    }
  }

  // To time series, in SI units or those of the instruments.
  const auto series = InverseTransform(grid, taper, spectra.Values(),
                                       sources.size() * receiverCount * 3);
  const auto samples = std::min<std::size_t>(
      grid.Samples(),
      static_cast<std::size_t>(sampling.duration / sampling.timeStep + 1e-9) +
          1);
  const auto unit = Unit(scales, calculation.quantity);
  auto seismograms =
      Seismograms(sources.size(), receiverCount, samples, sampling.timeStep);
  for (auto source : std::views::iota(std::size_t{0}, sources.size())) {
    for (auto receiver : std::views::iota(std::size_t{0}, receiverCount)) {
      for (auto component : {Vertical, North, East}) {
        const auto first =
            ((source * receiverCount + receiver) * 3 + component) *
            grid.Samples();
        const auto trace = std::span{series}.subspan(first, samples);
        for (auto i : std::views::iota(std::size_t{0}, samples)) {
          seismograms(source, receiver, component, i) = trace[i] * unit;
        }
      }
    }
  }
  return seismograms;
}

}  // namespace YSpec
