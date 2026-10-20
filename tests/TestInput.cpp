#include <gtest/gtest.h>

#include <YSpec/Input.hpp>
#include <filesystem>
#include <fstream>
#include <stdexcept>
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

// A minimal valid document, with members spliced in at the top and among the
// first source's.
std::string Document(const std::string& extra = "",
                     const std::string& source = "") {
  return R"({
  "model": "prem.200",)" +
         extra + R"(
  "output": { "prefix": "out" },
  "degrees": { "max": 10 },
  "frequencies_mHz": { "min": 0.2, "max": 10 },
  "time_series": { "length_min": 30, "time_step_s": 2 },
  "sources": [{)" +
         source + R"(
    "depth_km": 10, "latitude": 0, "longitude": 0,
    "moment_tensor_Nm": { "rr": 1, "rt": 0, "rp": 0, "tt": 0, "tp": 0, "pp": 0 }
  }],
  "receivers": { "depth_km": 0, "locations": [[10, 20]] }
})";
}

TEST(Input, ReadsTheExample) {
  const auto input = ReadInput(kExample);
  const auto& c = input.calculation;
  EXPECT_EQ(input.model, kExample.parent_path() / "prem.200");
  EXPECT_EQ(input.prefix, kExample.parent_path() / "yspec.out");
  EXPECT_EQ(c.gravity, Gravity::Full);
  EXPECT_TRUE(c.attenuation);
  EXPECT_EQ(c.quantity, Quantity::Velocity);
  EXPECT_EQ(c.degrees.max, 600);
  // In SI units.
  EXPECT_DOUBLE_EQ(c.band.max, 40e-3);
  EXPECT_DOUBLE_EQ(c.sampling.duration, 3600);
  ASSERT_TRUE(c.taper.has_value());
  EXPECT_DOUBLE_EQ(c.taper->f4, 40e-3);
  ASSERT_EQ(input.sources.size(), 1u);
  EXPECT_DOUBLE_EQ(input.sources[0].depth, 29.5e3);
  EXPECT_DOUBLE_EQ(input.sources[0].moment.tt, -2.976e18);
  EXPECT_EQ(input.names, std::vector<std::string>{"1"});
  EXPECT_EQ(input.receivers.size(), 60u);
  EXPECT_DOUBLE_EQ(input.receivers.back().location.longitude, 122);
  EXPECT_DOUBLE_EQ(input.receivers.back().depth, 3e3);
}

TEST(Input, AppliesTheDefaults) {
  const auto input = ReadInput(Write("defaults.json", Document()));
  const auto& c = input.calculation;
  EXPECT_EQ(c.gravity, Gravity::Full);
  EXPECT_TRUE(c.attenuation);
  EXPECT_EQ(c.quantity, Quantity::Displacement);
  EXPECT_FALSE(c.corrections);
  EXPECT_EQ(c.degrees.min, 0);
  EXPECT_FALSE(c.taper.has_value());
  EXPECT_EQ(c.solver.tolerance, SolverOptions{}.tolerance);
}

TEST(Input, ReadsSeveralNamedSources) {
  auto text = Document("", R"( "name": "a",)");
  const auto second = std::string(R"(,{ "name": "b",
    "depth_km": 600, "latitude": -10, "longitude": 100,
    "moment_tensor_Nm": { "rr": 0, "rt": 1, "rp": 0, "tt": 0, "tp": 0, "pp": 0 }
  }])");
  text.replace(text.find("}],"), 2, "}" + second);
  const auto input = ReadInput(Write("sources.json", text));
  ASSERT_EQ(input.sources.size(), 2u);
  EXPECT_EQ(input.names, (std::vector<std::string>{"a", "b"}));
  EXPECT_DOUBLE_EQ(input.sources[1].depth, 600e3);
  EXPECT_DOUBLE_EQ(input.sources[1].moment.rt, 1);
}

TEST(Input, RefusesAnUnknownMember) {
  const auto path =
      Write("unknown.json", Document(R"( "gravitation": "none",)"));
  try {
    (void)ReadInput(path);
    FAIL() << "no exception";
  } catch (const std::runtime_error& e) {
    EXPECT_NE(std::string(e.what()).find("gravitation"), std::string::npos)
        << e.what();
  }
}

TEST(Input, RefusesAnUnknownChoice) {
  const auto path = Write("choice.json", Document(R"( "gravity": "self",)"));
  try {
    (void)ReadInput(path);
    FAIL() << "no exception";
  } catch (const std::runtime_error& e) {
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
                 std::runtime_error)
        << extra;
  }
  auto text = Document();
  text.replace(text.find("[[10, 20]]"), 10, "[[100, 20]]");
  EXPECT_THROW((void)ReadInput(Write("latitude.json", text)),
               std::runtime_error);
  EXPECT_THROW(
      (void)ReadInput(Write("name.json", Document("", R"( "name": "a/b",)"))),
      std::runtime_error);
}

TEST(Input, RefusesAMissingMember) {
  auto text = Document();
  text.replace(text.find("\"degrees\""), std::string("\"degrees\"").size(),
               "\"degree\"");
  EXPECT_THROW((void)ReadInput(Write("missing.json", text)),
               std::runtime_error);
}

}  // namespace
