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

TEST(Run, WritesASeismogramFilePerSourceAndReceiver) {
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
  const auto seismograms = RunSeismograms(directory / "input.json");
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

TEST(Run, WritesASnapshotFilePerTime) {
  const auto directory = fs::temp_directory_path() / "yspec-run-snapshot-test";
  fs::remove_all(directory);
  fs::create_directories(directory);
  fs::copy_file(fs::path(YSPEC_EXAMPLES_DIR) / "prem.200",
                directory / "prem.200");
  std::ofstream(directory / "input.json") << R"({
  "model": "prem.200",
  "output": { "prefix": "snap" },
  "degrees": { "max": 20 },
  "frequencies_mHz": { "min": 0.2, "max": 4 },
  "time_series": { "length_min": 30 },
  "sources": [
    { "depth_km": 30, "latitude": 10, "longitude": 20,
      "moment_tensor_Nm": { "rr": 1e18, "rt": 0, "rp": 0, "tt": 0, "tp": 0,
                            "pp": 0 } }
  ],
  "snapshots": { "times_s": [600, 1200], "depths_km": [100, 500],
                 "distances_deg": { "min": 0, "max": 90, "count": 4 },
                 "azimuths_deg": [0, 180] }
})";
  const auto snapshots = RunSnapshots(directory / "input.json");
  for (auto time : {0, 1}) {
    const auto path = directory / std::format("snap.{}", time + 1);
    ASSERT_TRUE(fs::exists(path)) << path;
    auto file = std::ifstream(path);
    auto line = std::string{};
    std::getline(file, line);
    EXPECT_EQ(line.front(), '#');
    std::getline(file, line);
    EXPECT_EQ(line.front(), '#');
    auto count = std::size_t{0};
    for (auto row = std::array<double, 8>{};
         file >> row[0] >> row[1] >> row[2] >> row[3] >> row[4] >> row[5] >>
         row[6] >> row[7];
         ++count) {
      // Depth by depth, azimuth by azimuth, distance by distance.
      const auto depth = count / 8;
      const auto azimuth = (count / 4) % 2;
      const auto distance = count % 4;
      EXPECT_DOUBLE_EQ(row[0], depth == 0 ? 100 : 500);
      EXPECT_DOUBLE_EQ(row[1], 30.0 * static_cast<double>(distance));
      EXPECT_DOUBLE_EQ(row[2], azimuth == 0 ? 0 : 180);
      for (auto c : {Up, Away, Transverse}) {
        const auto value = snapshots(0, time, depth, distance, azimuth, c);
        EXPECT_NEAR(row[5 + c], value, 1e-6 * std::abs(value) + 1e-300);
      }
    }
    EXPECT_EQ(count, 16u);
  }
  // At the epicentre the location is the source's.
  auto file = std::ifstream(directory / "snap.1");
  auto line = std::string{};
  std::getline(file, line);
  std::getline(file, line);
  auto row = std::array<double, 8>{};
  file >> row[0] >> row[1] >> row[2] >> row[3] >> row[4];
  EXPECT_DOUBLE_EQ(row[3], 10);
  EXPECT_DOUBLE_EQ(row[4], 20);
}

}  // namespace
