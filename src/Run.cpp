#include <YSpec/Model.hpp>
#include <YSpec/Run.hpp>
#include <format>
#include <fstream>
#include <stdexcept>

namespace YSpec {

Seismograms Compute(const Input& input, const Progress& progress) {
  const auto model = Model::ReadDeck(input.model);
  return ComputeSeismograms(model, input.calculation, input.sources,
                            input.receivers, progress);
}

void WriteSeismograms(const Seismograms& seismograms, const Input& input) {
  for (auto source = std::size_t{0}; source < seismograms.Sources(); ++source) {
    for (auto receiver = std::size_t{0}; receiver < seismograms.Receivers();
         ++receiver) {
      auto path = input.prefix;
      path += seismograms.Sources() == 1
                  ? std::format(".{}", receiver + 1)
                  : std::format(".{}.{}", input.names.at(source), receiver + 1);
      auto file = std::ofstream(path);
      if (!file) throw std::runtime_error("cannot write " + path.string());
      const auto vertical = seismograms.Trace(source, receiver, Vertical);
      const auto north = seismograms.Trace(source, receiver, North);
      const auto east = seismograms.Trace(source, receiver, East);
      for (auto i = std::size_t{0}; i < seismograms.Samples(); ++i) {
        file << std::format("{:15.6E}{:15.6E}{:15.6E}{:15.6E}\n",
                            static_cast<double>(i) * seismograms.TimeStep(),
                            vertical[i], north[i], east[i]);
      }
      if (!file) throw std::runtime_error("cannot write " + path.string());
    }
  }
}

Seismograms Run(const std::filesystem::path& path, const Progress& progress) {
  const auto input = ReadInput(path);
  auto seismograms = Compute(input, progress);
  WriteSeismograms(seismograms, input);
  return seismograms;
}

}  // namespace YSpec
