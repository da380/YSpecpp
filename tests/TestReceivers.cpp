#include <gtest/gtest.h>

#include <YSpec/Receivers.hpp>
#include <cmath>
#include <numbers>
#include <vector>

namespace {

using namespace YSpec;
constexpr auto pi = std::numbers::pi;

TEST(Geometry, PlacesAReceiverEastOfTheSourceOnTheEquator) {
  const auto g = Geometry({0, 30}, {0, 0});
  EXPECT_NEAR(g.distance, pi / 6, 1e-12);
  // The receiver lies due east (azimuth 90 degrees), so phi = 180 - 90.
  EXPECT_NEAR(g.phi, pi / 2, 1e-12);
  // The source lies due west of the receiver (azimuth 270 degrees).
  EXPECT_NEAR(g.rotation, pi / 2, 1e-12);
}

TEST(Geometry, PlacesAReceiverNorthOfTheSource) {
  const auto g = Geometry({40, 10}, {10, 10});
  EXPECT_NEAR(g.phi, pi, 1e-12);       // azimuth 0
  EXPECT_NEAR(g.rotation, pi, 1e-12);  // the source is due south
  // Geocentric colatitudes make the distance a little less than 30 degrees.
  EXPECT_NEAR(g.distance * 180 / pi, 29.9, 0.1);
}

TEST(Geometry, HandlesCoincidentAndAntipodalPoints) {
  const auto same = Geometry({12, 34}, {12, 34});
  EXPECT_EQ(same.distance, 0.0);
  // As for a receiver approaching the source from any direction.
  EXPECT_NEAR(same.rotation, same.phi, 1e-15);
  const auto near = Geometry({12.000001, 34.000001}, {12, 34});
  EXPECT_NEAR(near.rotation, near.phi, 1e-8);
  const auto antipode = Geometry({0, 180}, {0, 0});
  EXPECT_NEAR(antipode.distance, pi, 1e-12);
  EXPECT_TRUE(std::isfinite(antipode.phi));
  EXPECT_TRUE(std::isfinite(antipode.rotation));
}

TEST(AngularFunctions, AgreeWithTheStandardLibrary) {
  const auto lMax = 60;
  const auto theta = std::vector{0.3, 1.2, 2.0, 3.0};
  const auto angular = AngularFunctions(lMax, theta);
  const auto h = 1e-6;
  for (auto i = 0; i < std::ssize(theta); ++i) {
    const auto t = theta[i];
    for (auto l = 0; l <= lMax; ++l) {
      for (auto m = -std::min(l, 2); m <= std::min(l, 2); ++m) {
        const auto sign = (m < 0 && m % 2 != 0) ? -1.0 : 1.0;
        auto X = [&](double angle) {
          return sign * std::sph_legendre(l, std::abs(m), angle);
        };
        const auto& a = angular(i, l, m);
        const auto scale = std::sqrt((2.0 * l + 1) / (4 * pi));
        EXPECT_NEAR(a.x, X(t), 1e-12 * scale) << l << " " << m;
        const auto dx = (X(t + h) - X(t - h)) / (2 * h);
        EXPECT_NEAR(a.dx, dx, 1e-6 * scale * (l + 1)) << l << " " << m;
        EXPECT_NEAR(a.mx, m * X(t) / std::sin(t), 1e-11 * scale * (l + 1))
            << l << " " << m;
      }
    }
  }
}

TEST(AngularFunctions, ReachBeyondTheLimitOfADoubleTable) {
  // Above GSHTrans's limit for double the table is built in long double.
  const auto lMax = 2500;
  const auto theta = std::vector{0.4, 1.9};
  const auto angular = AngularFunctions(lMax, theta);
  for (auto i : {0, 1}) {
    for (auto l : {10, 1000, 2000, 2500}) {
      for (auto m : {-2, 0, 1, 2}) {
        const auto sign = (m < 0 && m % 2 != 0) ? -1.0 : 1.0;
        const auto expected =
            sign * std::sph_legendre(l, std::abs(m), theta[i]);
        EXPECT_NEAR(angular(i, l, m).x, expected,
                    1e-10 * std::sqrt((2.0 * l + 1) / (4 * pi)))
            << l << " " << m;
      }
    }
  }
}

TEST(AngularFunctions, AreRegularAtThePoles) {
  const auto angular = AngularFunctions(10, {0.0, pi});
  for (auto i : {0, 1}) {
    for (auto l = 1; l <= 10; ++l) {
      // Only m = 0 survives at a pole, and m = +-1 in the derivatives.
      EXPECT_NEAR(angular(i, l, 2).x, 0, 1e-13);
      EXPECT_NEAR(angular(i, l, 2).mx, 0, 1e-13);
      EXPECT_TRUE(std::isfinite(angular(i, l, 1).mx));
      // At the north pole m X / sin(theta) -> -dX/dtheta for m = 1 ... in
      // magnitude both equal sqrt((2l+1)/4pi) sqrt(l(l+1))/2.
      const auto expected =
          std::sqrt((2.0 * l + 1) / (4 * pi)) * std::sqrt(l * (l + 1.0)) / 2;
      EXPECT_NEAR(std::abs(angular(i, l, 1).mx), expected, 1e-12);
      EXPECT_NEAR(std::abs(angular(i, l, 1).dx), expected, 1e-12);
    }
  }
}

}  // namespace
