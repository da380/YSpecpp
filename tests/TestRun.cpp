#include <gtest/gtest.h>

#include <YSpec/Run.hpp>
#include <array>
#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <string>

namespace {

using namespace YSpec;
namespace fs = std::filesystem;

TEST(Run, WritesAFilePerSourceAndReceiver) {
  const auto directory = fs::temp_directory_path() / "yspec-run-test";
  fs::remove_all(directory);
  fs::create_directories(directory);
  fs::copy_file(fs::path(YSPEC_EXAMPLES_DIR) / "prem.200",
                directory / "prem.200");
  std::ofstream(directory / "input.json") << R"({
  "model": "prem.200",
  "output": { "prefix": "out" },
  "degrees": { "max": 20 },
  "frequencies_mHz": { "min": 0.2, "max": 4 },
  "time_series": { "length_min": 30, "time_step_s": 60 },
  "sources": [
    { "name": "a", "depth_km": 30, "latitude": 0, "longitude": 0,
      "moment_tensor_Nm": { "rr": 1e18, "rt": 0, "rp": 0, "tt": 0, "tp": 0,
                            "pp": 0 } },
    { "name": "b", "depth_km": 100, "latitude": 10, "longitude": 10,
      "moment_tensor_Nm": { "rr": 0, "rt": 1e18, "rp": 0, "tt": 0, "tp": 0,
                            "pp": 0 } }
  ],
  "receivers": { "depth_km": 3, "locations": [[20, 30], [-20, 40]] }
})";
  const auto seismograms = YSpec::Run(directory / "input.json");
  ASSERT_EQ(seismograms.Sources(), 2u);
  ASSERT_EQ(seismograms.Receivers(), 2u);
  EXPECT_EQ(seismograms.Samples(), 31u);
  for (auto [s, name] : {std::pair{0, "a"}, {1, "b"}}) {
    for (auto k : {0, 1}) {
      const auto path = directory / std::format("out.{}.{}", name, k + 1);
      ASSERT_TRUE(fs::exists(path)) << path;
      auto file = std::ifstream(path);
      auto count = std::size_t{0};
      for (auto row = std::array<double, 4>{};
           file >> row[0] >> row[1] >> row[2] >> row[3]; ++count) {
        EXPECT_DOUBLE_EQ(row[0], 60.0 * count);
        for (auto c : {Vertical, North, East}) {
          const auto value = seismograms(s, k, c, count);
          EXPECT_NEAR(row[1 + c], value, 1e-6 * std::abs(value) + 1e-300);
        }
      }
      EXPECT_EQ(count, seismograms.Samples());
    }
  }
}

}  // namespace
