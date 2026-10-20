#include <gtest/gtest.h>

#include <YSpec/Kernels.hpp>
#include <YSpec/Synthetics.hpp>
#include <array>
#include <cmath>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace {

using namespace YSpec;
using P = RadialParameter;

const auto kPrem = std::string(YSPEC_EXAMPLES_DIR) + "/prem.200";

Calculation Small() {
  auto c = Calculation{};
  c.degrees = {0, 30};
  c.band = {0.2e-3, 4e-3};
  c.sampling = {3600, 60};
  c.taper = Taper{0.3e-3, 0.6e-3, 3e-3, 4e-3};
  c.quantity = Quantity::Velocity;
  c.solver.tolerance = 1e-10;
  c.solver.startDecay = 25;
  return c;
}

// The model with a set of properties scaled by 1 + h bump(r) at its knots.
// The splines through the knots are linear in the values, so the
// perturbation between them is the spline through h bump times the
// property, and gravity follows from the perturbed density.
Model Perturbed(const Model& model, const std::vector<Property>& properties,
                const std::function<double(double)>& bump, double h) {
  auto layers = std::vector<Layer>{};
  auto massBelow = 0.0;
  for (const auto& layer : model.Layers()) {
    const auto knots = layer.Knots();
    auto radii = std::vector<double>(knots.begin(), knots.end());
    auto values = std::vector<Properties>{};
    for (auto r : radii) {
      const auto m = layer(r);
      auto p = Properties{};
      p[Property::Rho] = m.rho;
      p[Property::A] = m.A;
      p[Property::C] = m.C;
      p[Property::F] = m.F;
      p[Property::L] = m.L;
      p[Property::N] = m.N;
      p[Property::Kappa] = m.kappa;
      p[Property::Mu] = m.mu;
      p[Property::QKappa] = m.qKappa;
      p[Property::QMu] = m.qMu;
      p[Property::QA] = m.qA;
      p[Property::QF] = m.qF;
      for (auto property : properties) p[property] *= 1 + h * bump(r);
      values.push_back(p);
    }
    layers.emplace_back(std::move(radii), values, layer.GetPhase(), massBelow);
    massBelow = layers.back().MassBelowTop();
  }
  return Model(model.GetScales(), model.ReferenceFrequency(), layers);
}

// The value of a property in a material.
double Value(const Material& m, Property property) {
  switch (property) {
    case Property::Rho:
      return m.rho;
    case Property::A:
      return m.A;
    case Property::C:
      return m.C;
    case Property::F:
      return m.F;
    case Property::L:
      return m.L;
    case Property::N:
      return m.N;
    case Property::QA:
      return m.qA;
    case Property::QF:
      return m.qF;
    case Property::QMu:
      return m.qMu;
    default:
      return 0;
  }
}

// The kernel of each property, as a radial parameter, and the scale of its
// SI unit.
P ParameterOf(Property property) {
  switch (property) {
    case Property::Rho:
      return P::Rho;
    case Property::A:
      return P::A;
    case Property::C:
      return P::C;
    case Property::F:
      return P::F;
    case Property::L:
      return P::L;
    case Property::N:
      return P::N;
    case Property::QA:
      return P::QA;
    case Property::QF:
      return P::QF;
    default:
      return P::QMu;
  }
}

double UnitOf(const Scales& scales, Property property) {
  switch (property) {
    case Property::Rho:
      return scales.Density();
    case Property::QA:
    case Property::QF:
    case Property::QMu:
      return 1;
    default:
      return scales.Modulus();
  }
}

struct Case {
  std::string name;
  std::vector<Property> properties;
  double centre;  // m
  double width;   // m
};

Source TheSource() {
  return {{10, 20},
          30e3,
          MomentTensor{.rr = 1e19,
                       .rt = 2e19,
                       .rp = -1e19,
                       .tt = -3e18,
                       .tp = 5e18,
                       .pp = -7e18},
          Force{}};
}

Receiver TheReceiver() { return {{-15, 50}, 3e3}; }

// A traveltime along a direction, with its kernels, and the finite
// differences that check them.
class Check {
 public:
  Check(Gravity gravity, Direction direction)
      : model_{Model::ReadDeck(kPrem)}, calculation_{Small()} {
    calculation_.gravity = gravity;
    const auto window = Window{.start = 700, .end = 1500, .taper = 200};
    measurement_ = Traveltime(model_, calculation_, TheSource(), TheReceiver(),
                              direction, window);
    kernels_ =
        ComputeRadialKernels(model_, calculation_, TheSource(), measurement_,
                             {.order = 5, .maxElement = 20e3});
  }

  // The change predicted by the kernels, and that by central differences.
  std::pair<double, double> operator()(const Case& c) const {
    const auto& scales = model_.GetScales();
    const auto length = scales.Length();
    auto bump = [&](double r) {
      const auto x = (r * length - c.centre) / c.width;
      return std::exp(-x * x);
    };
    const auto h = 1e-4;
    const auto plus = Perturbed(model_, c.properties, bump, h);
    const auto minus = Perturbed(model_, c.properties, bump, -h);
    const auto difference = (Measure(plus) - Measure(minus)) / (2 * h);
    // The perturbation of the properties at a radius, from the perturbed
    // splines.
    auto change = [&](Property property, double r) {
      const auto layer = model_.LayerIndex(r);
      return (Value(plus(layer, r), property) -
              Value(minus(layer, r), property)) /
             (2 * h) * UnitOf(scales, property);
    };
    auto predicted = 0.0;
    const auto& k = *kernels_;
    for (auto i = std::size_t{0}; i < k.Size(); ++i) {
      for (auto property : c.properties) {
        predicted += k.Weights()[i] * k[ParameterOf(property)][i] *
                     change(property, k.Radii()[i] / length);
      }
    }
    // And the point parts at the source.
    for (auto property : c.properties) {
      predicted += k.AtSource(ParameterOf(property)) *
                   change(property, k.SourceRadius() / length);
    }
    return {predicted, difference};
  }

 private:
  Model model_;
  Calculation calculation_;
  Measurement measurement_;
  std::optional<RadialKernels> kernels_;

  // The measurement, linearised: the sum of its adjoint source times the
  // seismogram.
  double Measure(const Model& model) const {
    const auto source = TheSource();
    const auto receiver = TheReceiver();
    const auto u = ComputeSeismograms(
        model, calculation_, std::span{&source, 1}, std::span{&receiver, 1});
    const auto& d = measurement_.direction;
    auto sum = 0.0;
    for (auto n = std::size_t{0}; n < u.Samples(); ++n) {
      const auto measured = d.up * u(0, 0, Vertical, n) +
                            d.north * u(0, 0, North, n) +
                            d.east * u(0, 0, East, n);
      sum += measurement_.adjointSource[n] * measured;
    }
    return sum * u.TimeStep();
  }
};

void ExpectAgreement(const Check& check, const std::vector<Case>& cases) {
  for (const auto& c : cases) {
    const auto [predicted, difference] = check(c);
    EXPECT_NEAR(predicted, difference, 1e-4 * std::abs(difference))
        << c.name << ": predicted " << predicted << ", by differences "
        << difference;
  }
}

TEST(RadialKernels, PredictFiniteDifferences) {
  ExpectAgreement(Check(Gravity::Full, Along(Vertical)),
                  {
                      {"rho, mantle", {Property::Rho}, 5800e3, 300e3},
                      {"A, mantle", {Property::A}, 5800e3, 300e3},
                      {"C, mantle", {Property::C}, 5800e3, 300e3},
                      {"F, mantle", {Property::F}, 5800e3, 300e3},
                      {"L, mantle", {Property::L}, 5800e3, 300e3},
                      {"N, mantle", {Property::N}, 5800e3, 300e3},
                      {"qA, mantle", {Property::QA}, 5800e3, 300e3},
                      {"qF, mantle", {Property::QF}, 5800e3, 300e3},
                      {"qMu, mantle", {Property::QMu}, 5800e3, 300e3},
                      // Around the source, where the kernels of C, F and L have
                      // a point part, and clear of the ocean, where only
                      // isotropic perturbations have meaning.
                      {"rho, crust", {Property::Rho}, 6340e3, 10e3},
                      {"C, crust", {Property::C}, 6340e3, 10e3},
                      {"F, crust", {Property::F}, 6340e3, 10e3},
                      {"L, crust", {Property::L}, 6340e3, 10e3},
                      {"rho, outer core", {Property::Rho}, 2500e3, 500e3},
                      {"kappa, outer core",
                       {Property::A, Property::C, Property::F},
                       2500e3,
                       500e3},
                      // Reached only through the gravity above it.
                      {"rho, inner core", {Property::Rho}, 600e3, 300e3},
                  });
}

TEST(RadialKernels, PredictFiniteDifferencesInTheCowlingApproximation) {
  ExpectAgreement(
      Check(Gravity::Cowling,
            RadialDirection(TheSource().location, TheReceiver().location)),
      {
          {"rho, mantle", {Property::Rho}, 5800e3, 300e3},
          {"L, mantle", {Property::L}, 5800e3, 300e3},
          {"N, mantle", {Property::N}, 5800e3, 300e3},
          {"rho, outer core", {Property::Rho}, 2500e3, 500e3},
          {"rho, inner core", {Property::Rho}, 600e3, 300e3},
      });
}

TEST(RadialKernels, PredictFiniteDifferencesWithoutGravity) {
  ExpectAgreement(
      Check(Gravity::None,
            TransverseDirection(TheSource().location, TheReceiver().location)),
      {
          {"rho, mantle", {Property::Rho}, 5800e3, 300e3},
          {"L, crust", {Property::L}, 6340e3, 10e3},
          {"N, mantle", {Property::N}, 5800e3, 300e3},
          {"rho, outer core", {Property::Rho}, 2500e3, 500e3},
      });
}

// The kernels of each set of parameters predict, for a perturbation of those
// parameters, what the model's own predict for the same perturbation: here
// with made-up kernels, at radii in the anisotropic upper mantle, the
// isotropic lower mantle and the fluid outer core.
TEST(RadialKernels, ReparameteriseConsistently) {
  const auto model = Model::ReadDeck(kPrem);
  const auto& scales = model.GetScales();
  auto radii = std::vector<double>{};
  for (auto r = 3000e3; r < 6300e3; r += 150e3) radii.push_back(r);
  auto kernels =
      RadialKernels(radii, std::vector<double>(radii.size(), 1e17), 6000e3);
  constexpr auto count = static_cast<std::size_t>(P::Count);
  for (auto p = std::size_t{0}; p < count; ++p) {
    const auto parameter = static_cast<P>(p);
    for (auto i = std::size_t{0}; i < radii.size(); ++i) {
      kernels[parameter][i] = std::sin(1.0 + p + 0.7 * i) * 1e-20;
    }
    kernels.AtSource(parameter) = std::cos(2.0 + p) * 1e-11;
  }
  // Small relative perturbations of each derived parameter at a radius.
  auto relative = [](std::size_t j, double r) {
    return 1e-3 * std::sin(0.3 * j + r / 1e6);
  };
  // The change they predict, through the model's own parameters and through
  // the set's.
  auto check = [&](ParameterSet set, auto toLove) {
    const auto named = Reparameterise(kernels, model, set, Scaling::Relative);
    auto direct = 0.0;
    auto derived = 0.0;
    auto add = [&](double r, double weight, auto kernel, auto derivedKernel) {
      const auto x = r / scales.Length();
      const auto m = model(model.LayerIndex(x), x);
      const auto love =
          toLove(m, [&](std::size_t j) { return relative(j, r); });
      for (auto p = std::size_t{0}; p < count; ++p) {
        direct += weight * kernel(static_cast<P>(p)) * love[p];
      }
      for (auto j = std::size_t{0}; j < named.names.size(); ++j) {
        derived += weight * derivedKernel(j) * relative(j, r);
      }
    };
    for (auto i = std::size_t{0}; i < radii.size(); ++i) {
      add(
          radii[i], kernels.Weights()[i], [&](P p) { return kernels[p][i]; },
          [&](std::size_t j) { return named.values[j][i]; });
    }
    add(
        kernels.SourceRadius(), 1.0, [&](P p) { return kernels.AtSource(p); },
        [&](std::size_t j) { return named.atSource[j]; });
    EXPECT_NEAR(derived, direct, 1e-10 * std::abs(direct));
  };
  const auto density = scales.Density();
  const auto modulus = scales.Modulus();
  // The changes of rho, A, C, F, L, N, qA, qF and qMu, in SI units.
  using Love = std::array<double, count>;
  // Isotropic changes of kappa and mu, and of qkappa and qmu, by Minos's
  // effective quality factors.
  auto isotropic = [&](const Material& m, double dRho, double dKappa,
                       double dMu, double dqKappa, double dqMu) {
    const auto x = 4 * m.mu / (3 * m.kappa + 4 * m.mu);
    return Love{dRho,
                dKappa + 4 * dMu / 3,
                dKappa + 4 * dMu / 3,
                dKappa - 2 * dMu / 3,
                dMu,
                dMu,
                (1 - x) * dqKappa + x * dqMu,
                ((1 - x) * dqKappa - x * dqMu / 2) / (1 - 1.5 * x),
                dqMu};
  };
  check(ParameterSet::Love, [&](const Material& m, auto d) {
    return Love{d(0) * m.rho * density,
                d(1) * m.A * modulus,
                d(2) * m.C * modulus,
                d(3) * m.F * modulus,
                d(4) * m.L * modulus,
                d(5) * m.N * modulus,
                d(6) * m.qA,
                d(7) * m.qF,
                d(8) * m.qMu};
  });
  check(ParameterSet::Moduli, [&](const Material& m, auto d) {
    return isotropic(m, d(0) * m.rho * density, d(1) * m.kappa * modulus,
                     d(2) * m.mu * modulus, d(3) * m.qKappa, d(4) * m.qMu);
  });
  check(ParameterSet::Velocities, [&](const Material& m, auto d) {
    // rho vp^2 = kappa + 4 mu / 3 and rho vs^2 = mu, to first order.
    const auto dMu = m.mu * (d(0) + 2 * d(2));
    const auto dP = (m.kappa + 4 * m.mu / 3) * (d(0) + 2 * d(1));
    return isotropic(m, d(0) * m.rho * density, (dP - 4 * dMu / 3) * modulus,
                     dMu * modulus, d(3) * m.qKappa, d(4) * m.qMu);
  });
  check(ParameterSet::Anisotropic, [&](const Material& m, auto d) {
    // rho, vpv, vph, vsv, vsh, eta.
    const auto dA = m.A * (d(0) + 2 * d(2));
    const auto dC = m.C * (d(0) + 2 * d(1));
    const auto dL = m.L * (d(0) + 2 * d(3));
    const auto dN = m.N * (d(0) + 2 * d(4));
    const auto eta = m.F / (m.A - 2 * m.L);
    const auto dF = eta * (dA - 2 * dL) + d(5) * eta * (m.A - 2 * m.L);
    auto love = isotropic(m, d(0) * m.rho * density, 0, 0, d(6) * m.qKappa,
                          d(7) * m.qMu);
    love[1] = dA * modulus;
    love[2] = dC * modulus;
    love[3] = dF * modulus;
    love[4] = dL * modulus;
    love[5] = dN * modulus;
    return love;
  });
}

}  // namespace
