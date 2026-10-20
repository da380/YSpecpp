#include <GaussQuad/All>
#include <YSpec/Kernels.hpp>
#include <YSpec/RadialSolver.hpp>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <numbers>
#include <ranges>
#include <stdexcept>
#include <string>
#include <utility>

#include "AngularFunctions.hpp"
#include "Coefficients.hpp"
#include "Geometry.hpp"
#include "KernelParts.hpp"
#include "OrderArray.hpp"
#include "SourceWeights.hpp"
#include "Spectra.hpp"
#include "TimeSeries.hpp"

namespace YSpec {

double Window::operator()(double t) const {
  if (t < start || t > end) return 0;
  auto rise = [this](double s) {
    return s >= taper ? 1.0
                      : 0.5 * (1 - std::cos(std::numbers::pi * s / taper));
  };
  return std::min(rise(t - start), rise(end - t));
}

Direction Along(Component component) {
  switch (component) {
    case Vertical:
      return {.up = 1, .north = 0, .east = 0};
    case North:
      return {.up = 0, .north = 1, .east = 0};
    case East:
      return {.up = 0, .north = 0, .east = 1};
  }
  return {};
}

// The azimuth of the source seen from the receiver is 2 pi less the angle
// that turns the components there; radial is opposite it, and transverse a
// quarter turn anticlockwise from radial.
Direction RadialDirection(const Location& source, const Location& receiver) {
  const auto back = 2 * std::numbers::pi - Geometry(receiver, source).rotation;
  return {.up = 0, .north = -std::cos(back), .east = -std::sin(back)};
}

Direction TransverseDirection(const Location& source,
                              const Location& receiver) {
  const auto back = 2 * std::numbers::pi - Geometry(receiver, source).rotation;
  return {.up = 0, .north = -std::sin(back), .east = std::cos(back)};
}

std::vector<Measurement> Traveltimes(const Model& model,
                                     const Calculation& calculation,
                                     const Source& source,
                                     const Receiver& receiver,
                                     std::span<const Direction> directions,
                                     std::span<const Window> windows) {
  if (directions.size() != windows.size()) {
    throw std::invalid_argument(
        "there must be a window for each direction of a traveltime");
  }
  const auto& scales = model.GetScales();
  const auto spectra = SeismogramSpectra(
      model, calculation, std::span{&source, 1}, std::span{&receiver, 1});
  const auto grid = MakeGrid(calculation, scales);
  const auto unit = Unit(scales, calculation.quantity) / scales.Time();
  const auto dt = calculation.sampling.timeStep;
  const auto samples = SampleCount(grid, calculation.sampling);
  auto measurements = std::vector<Measurement>{};
  for (auto i : std::views::iota(std::size_t{0}, directions.size())) {
    const auto& direction = directions[i];
    const auto& window = windows[i];
    // The rate of the seismogram, whose spectrum is i omega times its own at
    // the complex frequency: the damping is undone after the transform.
    auto rate = std::vector<Complex>(grid.Count());
    for (auto k : std::views::iota(std::size_t{0}, grid.Count())) {
      const auto measured = direction.up * spectra(0, 0, Vertical, k) +
                            direction.north * spectra(0, 0, North, k) +
                            direction.east * spectra(0, 0, East, k);
      rate[k] = TimesI(grid.Omega(k)) * measured;
    }
    const auto series =
        InverseTransform(grid, MakeTaper(calculation, scales), rate, 1);
    auto& measurement = measurements.emplace_back(Measurement{
        .receiver = receiver, .direction = direction, .adjointSource = {}});
    auto& j = measurement.adjointSource;
    auto norm = 0.0;
    for (auto n : std::views::iota(std::size_t{0}, samples)) {
      const auto w = window(static_cast<double>(n) * dt);
      const auto value = series[n] * unit;
      j.push_back(-w * value);
      norm += w * value * value * dt;
    }
    if (!(norm > 0)) {
      throw std::invalid_argument(
          "the synthetic is zero in the window of a traveltime");
    }
    for (auto& value : j) value /= norm;
  }
  return measurements;
}

Measurement Traveltime(const Model& model, const Calculation& calculation,
                       const Source& source, const Receiver& receiver,
                       const Direction& direction, const Window& window) {
  return std::move(Traveltimes(model, calculation, source, receiver,
                               std::span{&direction, 1}, std::span{&window, 1})
                       .front());
}

RadialKernels::RadialKernels(std::vector<double> radii,
                             std::vector<double> weights, double sourceRadius)
    : radii_{std::move(radii)},
      weights_{std::move(weights)},
      sourceRadius_{sourceRadius} {
  if (weights_.size() != radii_.size()) {
    throw std::invalid_argument(
        "RadialKernels: there must be a weight for every radius");
  }
  for (auto& values : values_) values.assign(radii_.size(), 0.0);
}

namespace {

// The mesh in radius, non-dimensional: elements between the discontinuities
// of the model and the given breaks, divided to be no larger than asked, each
// with
// the nodes of a Gauss-Legendre rule. Elements never cross a discontinuity,
// a source or a receiver, where the fields or their derivatives jump.
class Mesh {
 public:
  Mesh(const Model& model, const RadialMesh& options,
       std::span<const double> breaks)
      : order_{options.order} {
    if (order_ < 1 || !(options.maxElement > 0)) {
      throw std::invalid_argument(
          "the radial mesh needs at least one node in an element, and "
          "elements of positive size");
    }
    const auto rule = GaussQuad::GaussLegendreQuadrature1D<double>(order_);
    const auto& x = rule.Points();
    const auto& w = rule.Weights();
    const auto maxElement = options.maxElement / model.GetScales().Length();
    // Edges closer than this are one.
    const auto close = 1e-9 * model.SurfaceRadius();
    for (auto layer : std::views::iota(Int{0}, model.NumberOfLayers())) {
      const auto& here = model.GetLayer(layer);
      auto edges = std::vector{here.Lower(), here.Upper()};
      for (auto r : breaks) {
        if (r > edges.front() && r < edges.back()) edges.push_back(r);
      }
      std::ranges::sort(edges);
      auto lower = edges.front();
      for (auto upper : edges | std::views::drop(1)) {
        if (upper - lower < close) continue;
        const auto pieces = std::max(
            1,
            static_cast<int>(std::ceil((upper - lower) / maxElement - 1e-9)));
        const auto half = (upper - lower) / (2 * pieces);
        for (auto piece : std::views::iota(0, pieces)) {
          const auto centre = lower + (2 * piece + 1) * half;
          halves_.push_back(half);
          for (auto j : std::views::iota(0, order_)) {
            r_.push_back(centre + half * x[j]);
            weights_.push_back(half * w[j]);
            layers_.push_back(layer);
          }
        }
        lower = upper;
      }
    }
    // The integrals from each node to the top of its element of the
    // Lagrange polynomials through the nodes, on [-1, 1], by the same rule
    // mapped onto [x_i, 1], which is exact for them.
    auto basis = [&](int j, double t) {
      auto value = 1.0;
      for (auto k : std::views::iota(0, order_)) {
        if (k != j) value *= (t - x[k]) / (x[j] - x[k]);
      }
      return value;
    };
    above_.assign(static_cast<std::size_t>(order_ * order_), 0.0);
    for (auto i : std::views::iota(0, order_)) {
      const auto half = (1 - x[i]) / 2;
      for (auto j : std::views::iota(0, order_)) {
        auto sum = 0.0;
        for (auto q : std::views::iota(0, order_)) {
          sum += w[q] * basis(j, x[i] + half * (x[q] + 1));
        }
        above_[static_cast<std::size_t>(i * order_ + j)] = half * sum;
      }
    }
  }

  std::size_t Size() const { return r_.size(); }
  const std::vector<double>& Radii() const { return r_; }
  double Radius(std::size_t i) const { return r_[i]; }
  // The weight of a node in the integral over radius.
  double Weight(std::size_t i) const { return weights_[i]; }
  Int Layer(std::size_t i) const { return layers_[i]; }

  // The integral of f from each node to the top of the model.
  std::vector<double> IntegralAbove(std::span<const double> f) const {
    const auto order = static_cast<std::size_t>(order_);
    auto result = std::vector<double>(Size());
    auto tail = 0.0;
    for (auto element = halves_.size(); element-- > 0;) {
      const auto first = element * order;
      for (auto i : std::views::iota(std::size_t{0}, order)) {
        auto partial = 0.0;
        for (auto j : std::views::iota(std::size_t{0}, order)) {
          partial += above_[i * order + j] * f[first + j];
        }
        result[first + i] = tail + halves_[element] * partial;
      }
      for (auto j : std::views::iota(std::size_t{0}, order)) {
        tail += weights_[first + j] * f[first + j];
      }
    }
    return result;
  }

 private:
  int order_;
  std::vector<double> r_;
  std::vector<double> weights_;
  std::vector<Int> layers_;
  std::vector<double> halves_;  // half the width of each element
  std::vector<double> above_;
};

// The sum over orders of (-1)^m a[m] b[-m]. For fields of spin weights N and
// -N, expanded in X^N_lm exp(i m phi) and X^-N_lm exp(i m phi), this is the
// integral of their product over the unit sphere but for a factor (-1)^N,
// which the products below fold in; for the displacement, with u expanded as
// U Y r + V grad1 Y - W r x grad1 Y, the integral of u . u° is
// Pair(U, U°) + zeta^2 (Pair(V, V°) + Pair(W, W°)).
Complex Pair(const OrderArray<Complex>& a, const OrderArray<Complex>& b,
             int l) {
  auto sum = Complex{};
  for (auto m : OrderArray<Complex>::Orders(l)) {
    const auto term = a[m] * b[-m];
    sum += m % 2 == 0 ? term : -term;
  }
  return sum;
}

// The integrals over a sphere of the products of the forward and adjoint
// fields that the kernels need, summed over degree. With h the horizontal
// trace of the strain,
//
//   A: h h°,   C: e_rr e°_rr,   F: e_rr h° + h e°_rr,
//   L: 4 e_rh . e°_rh,   N: 2 e_hh : e°_hh - 2 h h°,
//
// the coefficients of the moduli in the elastic energy density, and
//
//   kinetic: s . s°,   radial: s_r s°_r,
//   coupling: r (s_h . grad_h s°_r + s°_h . grad_h s_r),
//   potential: s° . grad phi + s . grad phi°.
//
// The radial products are also kept apart at the degrees with a potential,
// where the gravitational terms include the density through Poisson's
// equation.
struct Products {
  Complex A, C, F, L, N;
  Complex kinetic;
  Complex radial;
  Complex radialWithPotential;
  Complex coupling;
  Complex potential;

  // They are linear in the adjoint field.
  Products& operator+=(const Products& p) {
    const auto mine = Members();
    const auto theirs = p.Members();
    for (auto i : std::views::iota(std::size_t{0}, mine.size())) {
      *mine[i] += *theirs[i];
    }
    return *this;
  }
  friend Products operator*(Complex factor, Products p) {
    for (auto* member : p.Members()) *member *= factor;
    return p;
  }

 private:
  std::array<Complex*, 10> Members() {
    return {&A,        &C,        &F,      &L,
            &N,        &kinetic,  &radial, &radialWithPotential,
            &coupling, &potential};
  }
  std::array<const Complex*, 10> Members() const {
    return {&A,        &C,        &F,      &L,
            &N,        &kinetic,  &radial, &radialWithPotential,
            &coupling, &potential};
  }
};

// The coefficients of a function of a field's coefficients, order by order.
template <class Function>
OrderArray<Complex> Each(int l, Function f) {
  auto result = OrderArray<Complex>{};
  for (auto m : OrderArray<Complex>::Orders(l)) result[m] = f(m);
  return result;
}

// The horizontal trace of the strain.
OrderArray<Complex> Trace(const Coefficients& c, int l, double r) {
  const auto zeta2 = l * (l + 1.0);
  return Each(l, [&](int m) { return (2.0 * c.u[m] - zeta2 * c.v[m]) / r; });
}

// The consoidal and toroidal parts of 2 e_rh, the shear strain on a
// horizontal plane.
OrderArray<Complex> Shear(const Coefficients& c, int l, double r) {
  return Each(l, [&](int m) { return c.dv[m] + (c.u[m] - c.v[m]) / r; });
}
OrderArray<Complex> Twist(const Coefficients& c, int l, double r) {
  return Each(l, [&](int m) { return c.dw[m] - c.w[m] / r; });
}

// Adds the products at one degree, at radius r. In a fluid the shear terms
// are left zero: the radial derivative of the tangential displacement is not
// given there, and the fluid has no shear modulus to perturb.
void Add(Products& p, const Coefficients& a, const Coefficients& b, int l,
         double r, bool solid, bool potential) {
  const auto zeta2 = l * (l + 1.0);
  auto pair = [l](const OrderArray<Complex>& x, const OrderArray<Complex>& y) {
    return Pair(x, y, l);
  };
  const auto ha = Trace(a, l, r);
  const auto hb = Trace(b, l, r);
  const auto hh = pair(ha, hb);
  const auto uu = pair(a.u, b.u);
  const auto tangential = pair(a.v, b.v) + pair(a.w, b.w);
  p.A += hh;
  p.C += pair(a.du, b.du);
  p.F += pair(a.du, hb) + pair(ha, b.du);
  if (solid) {
    p.L += zeta2 * (pair(Shear(a, l, r), Shear(b, l, r)) +
                    pair(Twist(a, l, r), Twist(b, l, r)));
    p.N += zeta2 * (zeta2 - 2) / (r * r) * tangential - hh;
  }
  p.kinetic += uu + zeta2 * tangential;
  p.radial += uu;
  p.coupling += zeta2 * (pair(a.v, b.u) + pair(a.u, b.v));
  if (potential) {
    p.radialWithPotential += uu;
    p.potential += pair(b.u, a.dphi) + pair(a.u, b.dphi) +
                   zeta2 * (pair(b.v, a.phi) + pair(a.v, b.phi)) / r;
  }
}

// The jumps in the forward field's displacement at the source, at radius r,
// from the weights of the unit jumps in rU, zeta r V and zeta r W.
Coefficients Jumps(const SourceWeights& weights, int l, double r) {
  auto jumps = Coefficients{};
  if (l == 0) {
    jumps.u[0] = weights.radial[0] / r;
    return jumps;
  }
  const auto zeta = std::sqrt(l * (l + 1.0));
  for (auto m : OrderArray<Complex>::Orders(l)) {
    jumps.u[m] = weights.spheroidal[m][0] / r;
    jumps.v[m] = weights.spheroidal[m][1] / (zeta * r);
    jumps.w[m] = weights.toroidal[m][0] / (zeta * r);
  }
  return jumps;
}

// Adds the products at one degree of the part of the forward field's strain
// concentrated on the source's radius r, the jumps in U, V and W times a delta
// function in r, with the adjoint field there.
void AddAtSource(Products& p, const Coefficients& jumps, const Coefficients& b,
                 int l, double r) {
  const auto zeta2 = l * (l + 1.0);
  p.C += Pair(jumps.u, b.du, l);
  p.F += Pair(jumps.u, Trace(b, l, r), l);
  p.L += zeta2 *
         (Pair(jumps.v, Shear(b, l, r), l) + Pair(jumps.w, Twist(b, l, r), l));
}

// What is accumulated at each node: the kernels of the parameters, and the
// two angular integrals from which the adjoint potential follows, M_rr and
// the integrand of its radial integral.
constexpr auto Parameters = static_cast<std::size_t>(RadialParameter::Count);
constexpr auto RadialMoment = Parameters;
constexpr auto Flux = Parameters + 1;
using Accumulator = RadialAccumulator;

// Adds the products at one frequency, with its weight in the sum over
// frequency, into the kernels: the sum is of the real parts, and each modulus
// is dispersed by its quality factor. The gravitational terms are those of the
// weight, through g, and where there is a potential those through its
// gradient and, in the gradient of gravity, the density.
void Accumulate(Accumulator& sum, const Products& p, const Material& m,
                double r, double weight, Complex dispersion, Complex omega2,
                bool gravity) {
  auto add = [&](std::size_t slot, Complex value) {
    sum[slot] += weight * value.real();
  };
  auto at = [](RadialParameter q) { return static_cast<std::size_t>(q); };
  add(at(RadialParameter::A), -(1.0 + m.qA * dispersion) * p.A);
  add(at(RadialParameter::C), -(1.0 + m.qA * dispersion) * p.C);
  add(at(RadialParameter::F), -(1.0 + m.qF * dispersion) * p.F);
  add(at(RadialParameter::L), -(1.0 + m.qMu * dispersion) * p.L);
  add(at(RadialParameter::N), -(1.0 + m.qMu * dispersion) * p.N);
  add(at(RadialParameter::QA), -dispersion * (m.A * p.A + m.C * p.C));
  add(at(RadialParameter::QF), -dispersion * m.F * p.F);
  add(at(RadialParameter::QMu), -dispersion * (m.L * p.L + m.N * p.N));
  auto weightTerms = Complex{};
  if (gravity) {
    const auto flux = (p.coupling - 4.0 * p.radial) / r;
    weightTerms = m.g * flux + 4 * PiG * m.rho * p.radialWithPotential;
    add(RadialMoment, m.rho * p.radialWithPotential);
    add(Flux, m.rho * flux);
  }
  add(at(RadialParameter::Rho), omega2 * p.kinetic - weightTerms - p.potential);
}

}  // namespace

std::vector<RadialAccumulation> AccumulateRadial(
    const Model& model, const Calculation& calculation, const Source& source,
    std::span<const Measurement> measurements, std::span<const Point> points,
    const Progress& progress) {
  const auto& scales = model.GetScales();
  const auto& degrees = calculation.degrees;
  const auto radii =
      CheckKernelInputs(model, calculation, source, measurements);
  const auto rr = radii.receiver;
  const auto grid = MakeGrid(calculation, scales);
  const auto nodes = points.size();
  const auto solver = KernelSolver(
      model, calculation, radii,
      std::vector<Point>(points.begin(), points.end()), radii.source);
  const auto receiverIndex = AdjointIndex(radii);

  // Where the receiver is in the frame with the source at the pole, the
  // functions there that place the adjoint field's force, and the forces of
  // the basis of the measurements' directions in that frame.
  const auto geometry =
      Geometry(measurements.front().receiver.location, source.location);
  const auto spin = SpinFunctions(degrees.max, {geometry.distance});
  const auto basis = MakeBasis(measurements);
  auto forces = std::vector<Force>{};
  for (const auto& direction : basis.directions) {
    forces.push_back(InSourceFrame(direction, geometry.rotation));
  }
  const auto bases = forces.size();
  const auto count = measurements.size();

  // The spectra of the measurements' forces: their weights in the sum over
  // frequency are the same.
  auto adjoint = std::vector<std::vector<Complex>>{};
  auto frequencyWeight = std::vector<double>{};
  for (const auto& measurement : measurements) {
    auto spectrum = MakeAdjointSpectrum(model, calculation, measurement);
    adjoint.push_back(std::move(spectrum.force));
    frequencyWeight = std::move(spectrum.weight);
  }

  // The model at the nodes and at the source.
  auto materials = std::vector<Material>{};
  auto solid = std::vector<bool>{};
  for (const auto& point : points) {
    materials.push_back(model(point.layer, point.r));
    solid.push_back(!model.GetLayer(point.layer).IsFluid());
  }
  const auto at = solver.Source(0);
  const auto atSource = model(at.layer, at.r);
  const auto gravity = calculation.gravity != Gravity::None;
  const auto full = calculation.gravity == Gravity::Full;
  const auto moment = source.moment / scales.Moment();
  const auto force = source.force / scales.Force();

  auto total = std::vector<std::vector<Accumulator>>(
      count, std::vector<Accumulator>(nodes));
  auto totalAtSource = std::vector<Accumulator>(count);
  auto done = std::atomic<std::size_t>{0};
  auto error = std::exception_ptr{};

#pragma omp parallel
  {
    auto local = std::vector<std::vector<Accumulator>>(
        count, std::vector<Accumulator>(nodes));
    auto localAtSource = std::vector<Accumulator>(count);
    // The products with the adjoint field of each force of the basis, with a
    // spectrum of one.
    auto products =
        std::vector<std::vector<Products>>(bases, std::vector<Products>(nodes));
    auto productsAtSource = std::vector<Products>(bases);
#pragma omp for schedule(dynamic)
    for (long index = 0; index < static_cast<long>(grid.Count()); ++index) {
      try {
        const auto k = static_cast<std::size_t>(index);
        const auto wanted = std::ranges::any_of(
            adjoint, [k](const auto& a) { return a[k] != Complex{}; });
        if (wanted) {
          const auto omega = grid.Omega(k);
          for (auto& row : products) std::ranges::fill(row, Products{});
          std::ranges::fill(productsAtSource, Products{});
          for (auto l = degrees.min; l <= degrees.max; ++l) {
            const auto c = Context(model, l, omega, calculation.attenuation);
            auto responses = Responses{};
            if (l == 0) {
              responses.radial = solver.Radial(omega);
            } else {
              responses.toroidal = solver.Toroidal(l, omega);
              responses.spheroidal = solver.Spheroidal(l, omega);
            }
            const auto forward =
                PoleWeights(c, c(at.layer, at.r), moment, force);
            auto x = OrderArray<std::array<double, 5>>{};
            for (auto m : OrderArray<Complex>::Orders(l)) x[m] = spin(0, l, m);
            const auto potential = full && l > 0;
            for (auto b : std::views::iota(std::size_t{0}, bases)) {
              const auto backward =
                  ForceWeights(c, rr, forces[b], x, geometry.phi);
              for (auto i : std::views::iota(std::size_t{0}, nodes)) {
                Add(products[b][i], Combine(responses, 0, i, l, forward),
                    Combine(responses, receiverIndex, i, l, backward), l,
                    points[i].r, solid[i], potential);
              }
              AddAtSource(productsAtSource[b], Jumps(forward, l, at.r),
                          Combine(responses, receiverIndex, nodes, l, backward),
                          l, at.r);
            }
          }
          const auto weight = frequencyWeight[k];
          const auto dispersion =
              Context(model, 0, omega, calculation.attenuation).Dispersion();
          const auto omega2 = omega * omega;
          for (auto j : std::views::iota(std::size_t{0}, count)) {
            if (adjoint[j][k] == Complex{}) continue;
            // This measurement's force, of its spectrum along its direction.
            auto by = std::vector<Complex>{};
            for (auto b : std::views::iota(std::size_t{0}, bases)) {
              by.push_back(adjoint[j][k] * basis.coordinates[j][b]);
            }
            auto combined = [&](auto&& part) {
              auto p = Products{};
              for (auto b : std::views::iota(std::size_t{0}, bases)) {
                p += by[b] * part(b);
              }
              return p;
            };
            for (auto i : std::views::iota(std::size_t{0}, nodes)) {
              Accumulate(local[j][i], combined([&](std::size_t b) {
                           return products[b][i];
                         }),
                         materials[i], points[i].r, weight, dispersion, omega2,
                         gravity);
            }
            Accumulate(localAtSource[j], combined([&](std::size_t b) {
                         return productsAtSource[b];
                       }),
                       atSource, at.r, weight, dispersion, omega2, false);
          }
        }
        const auto finished = ++done;
        if (progress) {
#pragma omp critical(YSpecProgress)
          progress(finished, grid.Count());
        }
      } catch (...) {
#pragma omp critical(YSpecError)
        if (!error) error = std::current_exception();
      }
    }
#pragma omp critical(YSpecKernels)
    for (auto j : std::views::iota(std::size_t{0}, count)) {
      for (auto i : std::views::iota(std::size_t{0}, nodes)) {
        for (auto slot : std::views::iota(std::size_t{0}, total[j][i].size())) {
          total[j][i][slot] += local[j][i][slot];
        }
      }
      for (auto slot : std::views::iota(std::size_t{0}, Parameters + 2)) {
        totalAtSource[j][slot] += localAtSource[j][slot];
      }
    }
  }
  if (error) std::rethrow_exception(error);
  auto result = std::vector<RadialAccumulation>{};
  for (auto j : std::views::iota(std::size_t{0}, count)) {
    result.push_back({.points = std::move(total[j]),
                      .atSource = totalAtSource[j],
                      .source = at});
  }
  return result;
}

std::vector<RadialKernels> ComputeRadialKernels(
    const Model& model, const Calculation& calculation, const Source& source,
    std::span<const Measurement> measurements, const RadialMesh& meshOptions,
    const Progress& progress) {
  const auto& scales = model.GetScales();
  const auto [rs, rr] =
      CheckKernelInputs(model, calculation, source, measurements);
  const auto length = scales.Length();

  // The kernels at every node of the mesh.
  const auto mesh = Mesh(model, meshOptions, std::array{rs, rr});
  const auto nodes = mesh.Size();
  auto points = std::vector<Point>{};
  for (auto i : std::views::iota(std::size_t{0}, nodes)) {
    points.push_back({mesh.Layer(i), mesh.Radius(i)});
  }
  const auto accumulations = AccumulateRadial(model, calculation, source,
                                              measurements, points, progress);
  const auto gravity = calculation.gravity != Gravity::None;

  auto result = std::vector<RadialKernels>{};
  for (const auto& [total, totalAtSource, at] : accumulations) {
    // The density's part through the equilibrium gravity: the adjoint
    // potential, which for spherically symmetric perturbations is
    // -4 pi G [M_rr(r) + integral from r to the surface of the flux].
    auto potential = std::vector<double>(nodes);
    if (gravity) {
      auto flux = std::vector<double>(nodes);
      for (auto i : std::views::iota(std::size_t{0}, nodes)) {
        flux[i] = total[i][Flux];
      }
      const auto above = mesh.IntegralAbove(flux);
      for (auto i : std::views::iota(std::size_t{0}, nodes)) {
        potential[i] = -4 * PiG * (total[i][RadialMoment] + above[i]);
      }
    }

    // In SI units. The point part at the source is the integral over radius
    // of the delta function there, with its r^2.
    auto radiiSI = std::vector<double>{};
    auto weights = std::vector<double>{};
    const auto volume = length * length * length;
    for (auto i : std::views::iota(std::size_t{0}, nodes)) {
      const auto r = mesh.Radius(i);
      radiiSI.push_back(r * length);
      weights.push_back(mesh.Weight(i) * r * r * volume);
    }
    auto& kernels = result.emplace_back(std::move(radiiSI), std::move(weights),
                                        at.r * length);
    for (auto parameter : std::views::iota(std::size_t{0}, Parameters)) {
      const auto p = static_cast<RadialParameter>(parameter);
      const auto scale = p == RadialParameter::Rho  ? scales.Density()
                         : p >= RadialParameter::QA ? 1.0
                                                    : scales.Modulus();
      auto values = kernels[p];
      for (auto i : std::views::iota(std::size_t{0}, nodes)) {
        const auto extra = p == RadialParameter::Rho ? potential[i] : 0.0;
        values[i] = (total[i][parameter] + extra) / (scale * volume);
      }
      kernels.AtSource(p) = at.r * at.r * totalAtSource[parameter] / scale;
    }
  }
  return result;
}

RadialKernels ComputeRadialKernels(const Model& model,
                                   const Calculation& calculation,
                                   const Source& source,
                                   const Measurement& measurement,
                                   const RadialMesh& mesh,
                                   const Progress& progress) {
  return std::move(ComputeRadialKernels(model, calculation, source,
                                        std::span{&measurement, 1}, mesh,
                                        progress)
                       .front());
}

NamedKernels Reparameterise(const RadialKernels& kernels, const Model& model,
                            ParameterSet set, Scaling scaling) {
  constexpr auto parameters = static_cast<std::size_t>(RadialParameter::Count);
  auto result = NamedKernels{.names = Names(set), .values = {}, .atSource = {}};
  const auto count = result.names.size();
  result.values.assign(count, std::vector<double>(kernels.Size()));
  auto transformed = [&](auto kernel, double r) {
    auto k = std::array<double, parameters>{};
    for (auto p : std::views::iota(std::size_t{0}, parameters)) {
      k[p] = kernel(static_cast<RadialParameter>(p));
    }
    const auto x = r / model.GetScales().Length();
    auto [values, parameter] =
        Transform(k, ValuesAt(model, model.LayerIndex(x), x), set);
    if (scaling == Scaling::Relative) {
      for (auto j : std::views::iota(std::size_t{0}, count)) {
        values[j] *= parameter[j];
      }
    }
    return values;
  };
  for (auto i : std::views::iota(std::size_t{0}, kernels.Size())) {
    const auto values = transformed(
        [&](RadialParameter p) { return kernels[p][i]; }, kernels.Radii()[i]);
    for (auto j : std::views::iota(std::size_t{0}, count)) {
      result.values[j][i] = values[j];
    }
  }
  result.atSource =
      transformed([&](RadialParameter p) { return kernels.AtSource(p); },
                  kernels.SourceRadius());
  return result;
}

KernelRadii CheckKernelInputs(const Model& model,
                              const Calculation& calculation,
                              const Source& source,
                              const Measurement& measurement) {
  const auto& degrees = calculation.degrees;
  if (degrees.min < 0 || degrees.max < degrees.min) {
    throw std::invalid_argument("the degrees must satisfy 0 <= min <= max");
  }
  if (calculation.corrections &&
      calculation.quantity == Quantity::Acceleration) {
    throw std::invalid_argument(
        "kernels of seismograms with the potential and tilt corrections are "
        "not available");
  }
  const auto& scales = model.GetScales();
  const auto samples =
      SampleCount(MakeGrid(calculation, scales), calculation.sampling);
  if (measurement.adjointSource.size() != samples) {
    throw std::invalid_argument(
        "the adjoint source must have a value at each of the " +
        std::to_string(samples) + " samples of the seismogram");
  }
  const auto length = scales.Length();
  const auto surface = model.SurfaceRadius();
  const auto radii =
      KernelRadii{.source = surface - source.depth / length,
                  .receiver = surface - measurement.receiver.depth / length};
  for (auto r : {radii.source, radii.receiver}) {
    if (!(r >= model.InnerRadius() && r <= surface)) {
      throw std::invalid_argument(
          "the source or receiver is outside the model");
    }
  }
  return radii;
}

KernelRadii CheckKernelInputs(const Model& model,
                              const Calculation& calculation,
                              const Source& source,
                              std::span<const Measurement> measurements) {
  if (measurements.empty()) {
    throw std::invalid_argument("there must be a measurement");
  }
  const auto& first = measurements.front().receiver;
  for (const auto& measurement : measurements) {
    const auto& receiver = measurement.receiver;
    if (receiver.location.latitude != first.location.latitude ||
        receiver.location.longitude != first.location.longitude ||
        receiver.depth != first.depth) {
      throw std::invalid_argument(
          "the measurements whose kernels are computed together must be at "
          "one receiver");
    }
  }
  auto radii = KernelRadii{};
  for (const auto& measurement : measurements) {
    radii = CheckKernelInputs(model, calculation, source, measurement);
  }
  return radii;
}

DirectionBasis MakeBasis(std::span<const Measurement> measurements) {
  auto basis = DirectionBasis{};
  auto dot = [](const Direction& a, const Direction& b) {
    return a.up * b.up + a.north * b.north + a.east * b.east;
  };
  // Gram-Schmidt, keeping what is left of each direction unless it is
  // negligible.
  for (const auto& measurement : measurements) {
    auto d = measurement.direction;
    const auto size = std::sqrt(dot(d, d));
    for (const auto& e : basis.directions) {
      const auto along = dot(d, e);
      d = {.up = d.up - along * e.up,
           .north = d.north - along * e.north,
           .east = d.east - along * e.east};
    }
    const auto left = std::sqrt(dot(d, d));
    if (left > 1e-9 * size) {
      basis.directions.push_back(
          {.up = d.up / left, .north = d.north / left, .east = d.east / left});
    }
  }
  for (const auto& measurement : measurements) {
    auto& coordinates = basis.coordinates.emplace_back();
    for (const auto& e : basis.directions) {
      coordinates.push_back(dot(measurement.direction, e));
    }
  }
  return basis;
}

RadialSolver KernelSolver(const Model& model, const Calculation& calculation,
                          const KernelRadii& radii, std::vector<Point> points,
                          std::optional<double> atSource) {
  auto sourceRadii = std::vector{radii.source};
  if (radii.receiver != radii.source) sourceRadii.push_back(radii.receiver);
  if (atSource) {
    // On the solid side of a discontinuity, as the source itself is.
    const auto r = *atSource;
    const auto below = model.LayerIndex(r, Side::Below);
    const auto above = model.LayerIndex(r, Side::Above);
    points.push_back({model.GetLayer(below).IsFluid() ? above : below, r});
  }
  try {
    return RadialSolver(model, std::move(sourceRadii), std::move(points),
                        calculation.gravity, calculation.attenuation,
                        calculation.solver);
  } catch (const std::invalid_argument& error) {
    throw std::invalid_argument(
        std::string("the adjoint field has its source at the receiver, "
                    "which must be in the solid region of the source: ") +
        error.what());
  }
}

AdjointSpectrum MakeAdjointSpectrum(const Model& model,
                                    const Calculation& calculation,
                                    const Measurement& measurement) {
  const auto& scales = model.GetScales();
  const auto grid = MakeGrid(calculation, scales);
  const auto samples = SampleCount(grid, calculation.sampling);
  const auto& j = measurement.adjointSource;
  // The spectrum of the adjoint field's force, A_k = w_k R_k sum_n j(t_n)
  // exp(epsilon t_n) exp(i omega_k t_n) dt, with w the taper and R the
  // derivatives of the quantity recorded and the instrument, and the weight of
  // each frequency in the sum that replaces the integral over frequency:
  // 2 df, but df at the Nyquist frequency, and nothing above it, as the
  // inverse transform has it. The adjoint source is made non-dimensional so
  // that the measurement keeps its units.
  const auto time = scales.Time();
  const auto toNonDimensional = Unit(scales, calculation.quantity) * time;
  const auto dt = grid.TimeStep();
  const auto taper = MakeTaper(calculation, scales);
  const auto& instrument = measurement.receiver.instrument;
  const auto response = instrument ? InstrumentResponse(*instrument, grid, time,
                                                        calculation.quantity)
                                   : std::vector<Complex>(grid.Count(), 1.0);
  const auto nyquist = grid.Samples() / 2;
  auto adjoint = std::vector<Complex>(grid.Count());
  auto frequencyWeight = std::vector<double>(grid.Count());
  for (auto k : std::views::iota(std::size_t{0}, grid.Count())) {
    const auto index = grid.First() + k;
    if (index > nyquist) continue;
    const auto omega = grid.Omega(k);
    auto sum = Complex{};
    for (auto n : std::views::iota(std::size_t{0}, samples)) {
      if (j[n] == 0) continue;
      const auto t = static_cast<double>(n) * dt;
      const auto phase = 2 * std::numbers::pi *
                         static_cast<double>((index * n) % grid.Samples()) /
                         static_cast<double>(grid.Samples());
      sum += j[n] * toNonDimensional * std::exp(grid.Damping() * t) *
             std::polar(1.0, phase);
    }
    const auto w = taper ? TaperValue(*taper, grid.Frequency(k)) : 1.0;
    adjoint[k] = w *
                 std::pow(TimesI(omega), Derivatives(calculation.quantity)) *
                 response[k] * sum * dt;
    frequencyWeight[k] = (index == nyquist ? 1.0 : 2.0) * grid.Spacing();
  }

  return {.force = std::move(adjoint), .weight = std::move(frequencyWeight)};
}

void Scale(SourceWeights& weights, Complex factor) {
  for (auto m : std::views::iota(-OrderArray<Complex>::MaxOrder,
                                 OrderArray<Complex>::MaxOrder + 1)) {
    for (auto& w : weights.toroidal[m]) w *= factor;
    for (auto& w : weights.spheroidal[m]) w *= factor;
  }
  for (auto& w : weights.radial) w *= factor;
}

// A force along a direction at the receiver, in (r, theta, phi) of the frame
// with the source at the pole, where the angle given turns (theta, phi) to
// (south, east).
Force InSourceFrame(const Direction& direction, double rotation) {
  const auto c = std::cos(rotation);
  const auto s = std::sin(rotation);
  return {.r = direction.up,
          .t = -c * direction.north + s * direction.east,
          .p = s * direction.north + c * direction.east};
}

ModelValues ValuesAt(const Model& model, Int layer, double x) {
  const auto& scales = model.GetScales();
  const auto m = model(layer, x);
  const auto modulus = scales.Modulus();
  return {.rho = m.rho * scales.Density(),
          .A = m.A * modulus,
          .C = m.C * modulus,
          .F = m.F * modulus,
          .L = m.L * modulus,
          .N = m.N * modulus,
          .kappa = m.kappa * modulus,
          .mu = m.mu * modulus,
          .qA = m.qA,
          .qF = m.qF,
          .qMu = m.qMu,
          .qKappa = m.qKappa};
}

std::vector<std::string> Names(ParameterSet set) {
  switch (set) {
    case ParameterSet::Love:
      return {"rho", "A", "C", "F", "L", "N", "qA", "qF", "qmu"};
    case ParameterSet::Moduli:
      return {"rho", "kappa", "mu", "qkappa", "qmu"};
    case ParameterSet::Velocities:
      return {"rho", "vp", "vs", "qkappa", "qmu"};
    case ParameterSet::Anisotropic:
      return {"rho", "vpv", "vph", "vsv", "vsh", "eta", "qkappa", "qmu"};
  }
  return {};
}

std::pair<std::vector<double>, std::vector<double>> Transform(
    const std::array<double, static_cast<std::size_t>(RadialParameter::Count)>&
        k,
    const ModelValues& v, ParameterSet set) {
  using P = RadialParameter;
  auto K = [&](P p) { return k[static_cast<std::size_t>(p)]; };
  const auto rho = K(P::Rho);
  const auto A = K(P::A);
  const auto C = K(P::C);
  const auto F = K(P::F);
  const auto L = K(P::L);
  const auto N = K(P::N);
  // Minos's effective inverse quality factors, qA = (1 - x) qkappa + x qmu
  // and qF = ((1 - x) qkappa - x qmu / 2) / (1 - 3x/2), with x = 4 mu /
  // (3 kappa + 4 mu), at fixed moduli.
  const auto x = 4 * v.mu / (3 * v.kappa + 4 * v.mu);
  const auto qKappa = (1 - x) * K(P::QA) + (1 - x) / (1 - 1.5 * x) * K(P::QF);
  const auto qMu =
      K(P::QMu) + x * K(P::QA) - 0.5 * x / (1 - 1.5 * x) * K(P::QF);
  // The Voigt averages perturbed isotropically.
  const auto kappa = A + C + F;
  const auto mu = 4.0 / 3 * (A + C) - 2.0 / 3 * F + L + N;
  switch (set) {
    case ParameterSet::Love:
      return {{rho, A, C, F, L, N, K(P::QA), K(P::QF), K(P::QMu)},
              {v.rho, v.A, v.C, v.F, v.L, v.N, v.qA, v.qF, v.qMu}};
    case ParameterSet::Moduli:
      return {{rho, kappa, mu, qKappa, qMu},
              {v.rho, v.kappa, v.mu, v.qKappa, v.qMu}};
    case ParameterSet::Velocities: {
      const auto vp = std::sqrt((v.kappa + 4 * v.mu / 3) / v.rho);
      const auto vs = std::sqrt(v.mu / v.rho);
      return {
          {rho + (v.kappa * kappa + v.mu * mu) / v.rho, 2 * v.rho * vp * kappa,
           2 * v.rho * vs * (mu - 4.0 / 3 * kappa), qKappa, qMu},
          {v.rho, vp, vs, v.qKappa, v.qMu}};
    }
    case ParameterSet::Anisotropic: {
      const auto vpv = std::sqrt(v.C / v.rho);
      const auto vph = std::sqrt(v.A / v.rho);
      const auto vsv = std::sqrt(v.L / v.rho);
      const auto vsh = std::sqrt(v.N / v.rho);
      const auto eta = v.F / (v.A - 2 * v.L);
      return {{rho + (v.A * A + v.C * C + v.F * F + v.L * L + v.N * N) / v.rho,
               2 * v.rho * vpv * C, 2 * v.rho * vph * (A + eta * F),
               2 * v.rho * vsv * (L - 2 * eta * F), 2 * v.rho * vsh * N,
               (v.A - 2 * v.L) * F, qKappa, qMu},
              {v.rho, vpv, vph, vsv, vsh, eta, v.qKappa, v.qMu}};
    }
  }
  return {};
}

}  // namespace YSpec
