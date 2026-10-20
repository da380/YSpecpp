// The responses of the radial solvers to unit jumps, compared with those of
// the Fortran YSpec (tests/data/solver, made by
// tests/fortran/make_solver_reference.sh with the Fortran's step density
// raised so that it is accurate to about 1e-8).

#include <gtest/gtest.h>

#include <YSpec/Solver.hpp>
#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <map>
#include <numbers>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace {

using namespace YSpec;

struct Reference {
  std::string model;
  bool attenuation = true;
  Gravity gravity = Gravity::Full;
  double sourceDepth = 0;    // km
  double receiverDepth = 0;  // km
  double length = 0;         // min
  // Keyed by kind ('T', 'S', 'R'), degree and frequency (mHz); the values
  // in the order the driver prints them.
  std::map<std::tuple<char, int, double>, std::vector<Complex>> values;
};

Reference Read(const std::string& name) {
  auto file = std::ifstream(std::string(YSPEC_TEST_DATA_DIR) + "/solver/" +
                            name + ".txt");
  if (!file) throw std::runtime_error("cannot open reference " + name);
  auto reference = Reference{};
  auto line = std::string{};
  while (std::getline(file, line)) {
    auto in = std::istringstream(line);
    auto tag = std::string{};
    in >> tag;
    if (tag == "#") {
      auto key = std::string{};
      in >> key;
      if (key == "model") in >> reference.model;
      if (key == "attenuation") {
        auto a = 0;
        in >> a;
        reference.attenuation = a == 1;
      }
      if (key == "gravity") {
        auto g = 0;
        in >> g;
        reference.gravity = g == 0   ? Gravity::None
                            : g == 1 ? Gravity::Cowling
                                     : Gravity::Full;
      }
      if (key == "source_depth_km") in >> reference.sourceDepth;
      if (key == "receiver_depth_km") in >> reference.receiverDepth;
      if (key == "time_series_min") in >> reference.length;
      continue;
    }
    auto l = 0;
    auto f = 0.0;
    in >> l >> f;
    if (tag == "S") {
      auto j = 0;
      in >> j;
    }
    auto& values = reference.values[{tag[0], l, f}];
    for (auto re = 0.0, im = 0.0; in >> re >> im;) values.emplace_back(re, im);
  }
  return reference;
}

class SolverAgainstFortran : public testing::TestWithParam<std::string> {};

TEST_P(SolverAgainstFortran, AgreesWithTheFortran) {
  const auto reference = Read(GetParam());
  const auto model =
      Model::ReadDeck(std::string(YSPEC_EXAMPLES_DIR) + "/" + reference.model);
  const auto& scales = model.GetScales();
  const auto surface = model.SurfaceRadius();
  const auto rs = surface - reference.sourceDepth * 1e3 / scales.Length();
  const auto rr = surface - reference.receiverDepth * 1e3 / scales.Length();
  // The imaginary part of the frequency, as in the original YSpec.
  const auto epsilon = 5 / (reference.length * 60 / scales.Time());
  auto options = SolverOptions{};
  options.integration.tolerance = 1e-10;
  const auto solver = RadialSolver<>(model, rs, rr, reference.gravity,
                                     reference.attenuation, options);

  for (const auto& [key, expected] : reference.values) {
    const auto [kind, l, f] = key;
    const auto omega =
        Complex(2 * std::numbers::pi * f / 1000 * scales.Time(), -epsilon);
    auto computed = std::vector<Complex>{};
    if (kind == 'R') {
      const auto y = solver.Radial(omega);
      computed.assign(y.begin(), y.end());
    } else if (kind == 'T') {
      const auto y = solver.Toroidal(l, omega);
      computed.assign(y.begin(), y.end());
    } else {
      const auto y = solver.Spheroidal(l, omega);
      // The driver prints, for each jump, U, V and the potential.
      for (auto j : {0, 1, 2, 3}) {
        for (auto i : {0, 1, 2}) computed.push_back(y[i][j]);
      }
    }
    ASSERT_EQ(computed.size(), expected.size());
    auto scale = 0.0;
    for (const auto& v : expected) scale = std::max(scale, std::abs(v));
    for (auto i = std::size_t{0}; i < expected.size(); ++i) {
      EXPECT_LT(std::abs(computed[i] - expected[i]), 1e-5 * scale)
          << kind << " l = " << l << " f = " << f << " mHz, entry " << i << ": "
          << computed[i] << " against " << expected[i];
    }
  }
}

INSTANTIATE_TEST_SUITE_P(Cases, SolverAgainstFortran,
                         testing::Values("full", "full_below", "cowling",
                                         "cowling_below", "none", "none_below",
                                         "fic_full", "fic_cowling", "fic_none",
                                         "no_ocean_full", "no_ocean_none"));

TEST(Solver, RefusesASourceInTheOcean) {
  const auto model =
      Model::ReadDeck(std::string(YSPEC_EXAMPLES_DIR) + "/prem.200");
  const auto surface = model.SurfaceRadius();
  const auto oceanFloor = model.Regions().back().lower;
  EXPECT_THROW(
      RadialSolver<>(model, (surface + oceanFloor) / 2, surface, Gravity::Full),
      std::invalid_argument);
}

TEST(Solver, PutsASourceOnTheCoreMantleBoundaryInTheMantle) {
  const auto model =
      Model::ReadDeck(std::string(YSPEC_EXAMPLES_DIR) + "/prem.200");
  const auto cmb = model.Regions()[1].upper;
  const auto solver =
      RadialSolver<>(model, cmb, model.Regions()[2].upper, Gravity::Full);
  EXPECT_EQ(solver.Source().layer, model.Regions()[2].first);
  const auto y = solver.Spheroidal(10, Complex(20.0, -1.0));
  EXPECT_TRUE(std::isfinite(std::abs(y[0][0])));
}

TEST(Solver, HasNoGravityFluxAtAFreeSurface) {
  // The flux variable psi vanishes at the free surface, which is what makes
  // the potential correction there the original YSpec's.
  const auto model =
      Model::ReadDeck(std::string(YSPEC_EXAMPLES_DIR) + "/prem.200.no");
  const auto surface = model.SurfaceRadius();
  const auto solver =
      RadialSolver<>(model, surface - 0.005, surface, Gravity::Full);
  for (auto l : {2, 20}) {
    const auto y = solver.Spheroidal(l, Complex(20.0, -0.5));
    for (auto j : {0, 1, 2, 3}) {
      EXPECT_LT(std::abs(y[3][j]), 1e-9 * std::abs(y[2][j])) << l << " " << j;
    }
  }
}

TEST(Solver, RefusesOptionsOutOfRange) {
  const auto model =
      Model::ReadDeck(std::string(YSPEC_EXAMPLES_DIR) + "/prem.200");
  auto options = SolverOptions{};
  options.integration.tolerance = 0;
  EXPECT_THROW(RadialSolver<>(model, 0.99, 0.995, Gravity::Full, true, options),
               std::invalid_argument);
}

TEST(Solver, StartsDeeperForLowerDegreesAndHigherFrequencies) {
  const auto model =
      Model::ReadDeck(std::string(YSPEC_EXAMPLES_DIR) + "/prem.200");
  const auto surface = model.SurfaceRadius();
  const auto solver =
      RadialSolver<>(model, surface - 0.005, surface - 0.001, Gravity::Full);
  const auto omega = 2 * std::numbers::pi * 0.02 * model.GetScales().Time();
  const auto top = solver.Source();
  const auto high = solver.StartPoint(600, omega, top, 0);
  const auto low = solver.StartPoint(100, omega, top, 0);
  EXPECT_GT(high.r, low.r);
  EXPECT_GT(high.r, 0.5);
  // Below about l = 10 at 20 mHz nothing is evanescent enough.
  EXPECT_EQ(solver.StartPoint(2, omega, top, 0).r, 0.0);
}

}  // namespace
