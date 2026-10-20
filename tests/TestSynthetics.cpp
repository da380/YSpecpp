#include <gtest/gtest.h>

#include <YSpec/Synthetics.hpp>
#include <algorithm>
#include <cmath>
#include <string>

namespace {

using namespace YSpec;

Settings Small() {
  auto s = Settings{};
  s.lMax = 40;
  s.fMin = 0.2;
  s.fMax = 5;
  s.length = 60;
  s.timeStep = 60;
  s.sourceDepth = 30;
  s.receiverDepth = 3;
  s.source = {20, 30};
  s.moment = {
      .rr = 1e18, .rt = 1e19, .rp = 5e18, .tt = -2e18, .tp = 3e18, .pp = 1e18};
  return s;
}

const auto kPrem = std::string(YSPEC_EXAMPLES_DIR) + "/prem.200";

TEST(Synthetics, AreContinuousAtTheEpicentre) {
  const auto model = Model::ReadDeck(kPrem);
  auto settings = Small();
  const auto d = 0.001;
  settings.receivers = {
      {20, 30}, {20 + d, 30}, {20 - d, 30}, {20, 30 + d}, {20, 30 - d}};
  const auto s = ComputeSeismograms(model, settings);
  for (auto c : {0, 1, 2}) {
    auto peak = 0.0;
    for (auto i = std::size_t{0}; i < s.Samples(); ++i) {
      peak = std::max(peak, std::abs(s(0, c, i)));
    }
    ASSERT_GT(peak, 0);
    for (auto k : {1, 2, 3, 4}) {
      for (auto i = std::size_t{0}; i < s.Samples(); ++i) {
        EXPECT_NEAR(s(k, c, i), s(0, c, i), 1e-2 * peak)
            << "component " << c << ", receiver " << k << ", sample " << i;
      }
    }
  }
}

TEST(Synthetics, ScaleWithTheMomentTensor) {
  const auto model = Model::ReadDeck(kPrem);
  auto settings = Small();
  settings.receivers = {{-10, 60}, {50, 100}};
  const auto a = ComputeSeismograms(model, settings);
  settings.moment = {.rr = 3 * settings.moment.rr,
                     .rt = 3 * settings.moment.rt,
                     .rp = 3 * settings.moment.rp,
                     .tt = 3 * settings.moment.tt,
                     .tp = 3 * settings.moment.tp,
                     .pp = 3 * settings.moment.pp};
  const auto b = ComputeSeismograms(model, settings);
  for (auto k : {0, 1}) {
    for (auto c : {0, 1, 2}) {
      for (auto i = std::size_t{0}; i < a.Samples(); ++i) {
        EXPECT_NEAR(b(k, c, i), 3 * a(k, c, i),
                    1e-12 * std::abs(b(k, c, i)) + 1e-30);
      }
    }
  }
}

}  // namespace
