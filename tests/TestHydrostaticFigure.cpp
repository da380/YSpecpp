#include <gtest/gtest.h>

#include <YSpec/HydrostaticFigure.hpp>
#include <YSpec/Kernels.hpp>
#include <YSpec/Model.hpp>
#include <algorithm>
#include <limits>
#include <numbers>
#include <vector>

namespace {

using namespace YSpec;

const auto kPrem = std::string(YSPEC_EXAMPLES_DIR) + "/prem.200";

Model UniformDensity(double radius, double density, double inner = 0) {
  auto values = std::vector<Properties>(2);
  for (auto& p : values) {
    p[Property::Rho] = 1;
    p[Property::A] = p[Property::C] = p[Property::F] = 1;
    p[Property::L] = p[Property::N] = p[Property::Kappa] = 1;
    p[Property::Mu] = p[Property::QKappa] = p[Property::QMu] = 0;
    p[Property::QA] = p[Property::QF] = 0;
  }
  auto layer = Layer({inner, 1}, values, Phase::Solid, 0);
  return Model(Scales{radius, density}, 1, {std::move(layer)});
}

TEST(HydrostaticFigure, UniformDensityMatchesAnalyticProfile) {
  const auto model = UniformDensity(6.371e6, 5515);
  const auto figure = HydrostaticFigure::FromModel(model);
  const auto g = model(0, 1).g * model.GetScales().Acceleration();
  const auto radius = model.SurfaceRadius() * model.GetScales().Length();
  const auto m = EarthRotationRate * EarthRotationRate * radius / g;
  EXPECT_NEAR(figure.SurfaceEllipticity(), 1.25 * m, 2e-11 * m);
  EXPECT_NEAR(figure.Ellipticity(0), figure.SurfaceEllipticity(), 2e-11 * m);
  EXPECT_NEAR(figure.Ellipticity(radius / 2), figure.SurfaceEllipticity(),
              2e-11 * m);
}

TEST(HydrostaticFigure, PremSurfaceMatchesIndependentClairautReference) {
  const auto model = Model::ReadDeck(kPrem);
  const auto figure = HydrostaticFigure::FromModel(model);
  EXPECT_NEAR(figure.SurfaceEllipticity(), 0.003334826529126116, 2e-11);
  EXPECT_GT(figure.Nodes().size(), 200u);
  EXPECT_NEAR(figure.Ellipticity(0), figure.Nodes().front().ellipticity, 1e-15);
  for (auto i = std::size_t{1}; i < figure.Nodes().size(); ++i) {
    EXPECT_GT(figure.Nodes()[i].radius, figure.Nodes()[i - 1].radius);
  }
}

TEST(HydrostaticFigure, ProfileAndDerivativeAreContinuousAtLayerInterfaces) {
  const auto model = Model::ReadDeck(kPrem);
  const auto figure = HydrostaticFigure::FromModel(model);
  const auto fine =
      HydrostaticFigure::FromModel(model, EarthRotationRate, 1e-11);
  EXPECT_NEAR(fine.SurfaceEllipticity(), figure.SurfaceEllipticity(), 2e-11);
  const auto length = model.GetScales().Length();
  for (auto layer = Int{0}; layer + 1 < model.NumberOfLayers(); ++layer) {
    const auto radius = model.GetLayer(layer).Upper() * length;
    const auto node = std::ranges::find(figure.Nodes(), radius,
                                        &HydrostaticFigure::Node::radius);
    ASSERT_NE(node, figure.Nodes().end());
    constexpr auto h = 10.0;
    const auto left =
        (figure.Ellipticity(radius) - figure.Ellipticity(radius - h)) / h;
    const auto right =
        (figure.Ellipticity(radius + h) - figure.Ellipticity(radius)) / h;
    EXPECT_NEAR(left, right,
                1e-9 * std::max(std::abs(left), std::abs(right)) + 1e-14);
    EXPECT_NEAR(node->ellipticity, fine.Ellipticity(radius), 2e-10);
  }
}

TEST(HydrostaticFigure, RotationScalingAndSphericalIdentity) {
  const auto model = Model::ReadDeck(kPrem);
  const auto earth = HydrostaticFigure::FromModel(model);
  const auto half = HydrostaticFigure::FromModel(model, EarthRotationRate / 2);
  EXPECT_NEAR(half.SurfaceEllipticity() / earth.SurfaceEllipticity(), 0.25,
              1e-12);
  const auto sphere = HydrostaticFigure::Spherical();
  const auto zeroSpin = HydrostaticFigure::FromModel(model, 0);
  EXPECT_FALSE(sphere.Enabled());
  EXPECT_TRUE(zeroSpin.Enabled());
  EXPECT_DOUBLE_EQ(sphere.GeocentricLatitude(45), 45);
  EXPECT_DOUBLE_EQ(sphere.GeographicLatitude(-90), -90);
  for (const auto latitude : {-90.0, -45.0, 13.37, 45.0, 90.0}) {
    EXPECT_DOUBLE_EQ(zeroSpin.GeocentricLatitude(latitude), latitude);
    EXPECT_DOUBLE_EQ(zeroSpin.GeographicLatitude(latitude), latitude);
  }
  EXPECT_DOUBLE_EQ(sphere.Ellipticity(3e6), 0);
  EXPECT_THROW(sphere.Ellipticity(std::numeric_limits<double>::quiet_NaN()),
               std::invalid_argument);
}

TEST(HydrostaticFigure, SurfaceConversionsAreInverseAndPoleSafe) {
  const auto figure = HydrostaticFigure::FromModel(Model::ReadDeck(kPrem));
  for (auto latitude : {-90.0, -80.0, -13.82, 0.0, 45.0, 90.0}) {
    EXPECT_NEAR(figure.GeographicLatitude(figure.GeocentricLatitude(latitude)),
                latitude, 2e-14);
  }
  EXPECT_DOUBLE_EQ(figure.GeocentricLatitude(90), 90);
  EXPECT_DOUBLE_EQ(figure.GeocentricLatitude(-90), -90);
}

TEST(HydrostaticFigure, ModelAwareSamplingDefaultsToTheModelFigure) {
  const auto model = Model::ReadDeck(kPrem);
  const auto figure = HydrostaticFigure::FromModel(model);
  const SectionOptions options{
      .margin = 1, .halfWidth = 1, .spacing = 1, .radialSpacing = 1e6};
  const auto implicit = PathSection(model, {-13.82, -67.25}, {80, 0}, options);
  const auto explicitFigure =
      PathSection(model, {-13.82, -67.25}, {80, 0}, options, figure);
  ASSERT_EQ(implicit.places.size(), explicitFigure.places.size());
  for (auto i = std::size_t{0}; i < implicit.places.size(); ++i) {
    EXPECT_DOUBLE_EQ(implicit.places[i].colatitude,
                     explicitFigure.places[i].colatitude);
    EXPECT_DOUBLE_EQ(implicit.places[i].longitude,
                     explicitFigure.places[i].longitude);
  }
}

TEST(HydrostaticFigure, RejectsHollowModelsAndFastRotation) {
  const auto hollow = UniformDensity(6.371e6, 5515, 0.1);
  EXPECT_THROW(HydrostaticFigure::FromModel(hollow), std::invalid_argument);
  const auto model = Model::ReadDeck(kPrem);
  EXPECT_THROW(HydrostaticFigure::FromModel(model, 1), std::invalid_argument);
}

}  // namespace
