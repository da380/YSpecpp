#include <gtest/gtest.h>

#include <YSpec/Input.hpp>
#include <filesystem>
#include <fstream>
#include <string>

namespace {

using namespace YSpec;
namespace fs = std::filesystem;

const auto kExample = fs::path(YSPEC_EXAMPLES_DIR) / "input.json";

// Writes a document to a file of its own and returns the path.
fs::path Write(const std::string& name, const std::string& text) {
  const auto directory = fs::temp_directory_path() / "yspec-input-tests";
  fs::create_directories(directory);
  const auto path = directory / name;
  std::ofstream(path) << text;
  return path;
}

// A minimal valid document, with an extra member spliced in at the top.
std::string Document(const std::string& extra = "") {
  return R"({
  "model": "prem.200",)" +
         extra + R"(
  "output": { "prefix": "out" },
  "degrees": { "max": 10 },
  "frequencies_mHz": { "min": 0.2, "max": 10 },
  "time_series": { "length_min": 30, "time_step_s": 2 },
  "source": {
    "depth_km": 10, "latitude": 0, "longitude": 0,
    "moment_tensor_Nm": { "rr": 1, "rt": 0, "rp": 0, "tt": 0, "tp": 0, "pp": 0 }
  },
  "receivers": { "depth_km": 0, "locations": [[10, 20]] }
})";
}

TEST(Input, ReadsTheExample) {
  const auto input = ReadInput(kExample);
  const auto& s = input.settings;
  EXPECT_EQ(input.model, kExample.parent_path() / "prem.200");
  EXPECT_EQ(input.prefix, kExample.parent_path() / "yspec.out");
  EXPECT_EQ(s.gravity, Gravity::Full);
  EXPECT_TRUE(s.attenuation);
  EXPECT_EQ(s.quantity, Quantity::Velocity);
  EXPECT_EQ(s.lMax, 600);
  EXPECT_DOUBLE_EQ(s.fMax, 40);
  EXPECT_DOUBLE_EQ(s.moment.tt, -2.976e18);
  EXPECT_EQ(s.receivers.size(), 60u);
  EXPECT_DOUBLE_EQ(s.receivers.back().longitude, 122);
}

TEST(Input, AppliesTheDefaults) {
  const auto input = ReadInput(Write("defaults.json", Document()));
  const auto& s = input.settings;
  EXPECT_EQ(s.gravity, Gravity::Full);
  EXPECT_TRUE(s.attenuation);
  EXPECT_EQ(s.quantity, Quantity::Displacement);
  EXPECT_FALSE(s.corrections);
  EXPECT_EQ(s.lMin, 0);
  EXPECT_EQ(s.filter, (std::array<double, 4>{}));
  EXPECT_EQ(s.solver.integration.tolerance,
            SolverOptions{}.integration.tolerance);
}

TEST(Input, RefusesAnUnknownMember) {
  const auto path =
      Write("unknown.json", Document(R"( "gravitation": "none",)"));
  try {
    (void)ReadInput(path);
    FAIL() << "no exception";
  } catch (const planetmodel::Error& e) {
    EXPECT_NE(std::string(e.what()).find("gravitation"), std::string::npos)
        << e.what();
  }
}

TEST(Input, RefusesAnUnknownChoice) {
  const auto path = Write("choice.json", Document(R"( "gravity": "self",)"));
  try {
    (void)ReadInput(path);
    FAIL() << "no exception";
  } catch (const planetmodel::Error& e) {
    const auto message = std::string(e.what());
    EXPECT_NE(message.find("cowling"), std::string::npos) << message;
  }
}

TEST(Input, RefusesValuesOutOfRange) {
  for (const auto& extra :
       {R"( "solver": { "tolerance": -1e-6 },)",
        R"( "solver": { "start_decay": 0 },)",
        R"( "filter_mHz": [5, 5, 5, 5],)", R"( "filter_mHz": [1, 2, 4, 3],)"}) {
    EXPECT_THROW((void)ReadInput(Write("range.json", Document(extra))),
                 planetmodel::Error)
        << extra;
  }
  auto text = Document();
  text.replace(text.find("[[10, 20]]"), 10, "[[100, 20]]");
  EXPECT_THROW((void)ReadInput(Write("latitude.json", text)),
               planetmodel::Error);
}

TEST(Input, RefusesAMissingMember) {
  auto text = Document();
  text.replace(text.find("\"degrees\""), std::string("\"degrees\"").size(),
               "\"degree\"");
  EXPECT_THROW((void)ReadInput(Write("missing.json", text)),
               planetmodel::Error);
}

}  // namespace
