// yspec-kernels: the sensitivity of a traveltime measured on a seismogram to
// spherically symmetric perturbations of the earth model, by the adjoint
// method.
//
//   yspec-kernels input.json
//
// Reads the calculation from a JSON file (see YSpec/Input.hpp for its
// members), computes the kernels and writes them (see YSpec/Run.hpp for the
// files), reporting progress on the standard error.

#include <YSpec/YSpec.hpp>
#include <filesystem>
#include <format>
#include <iostream>
#include <string>

#include "Program.hpp"

namespace {

const auto kUsage =
    "usage: yspec-kernels input.json\n"
    "\n"
    "Computes the radial kernels of a cross-correlation traveltime measured\n"
    "on the seismogram of one source at one receiver, as input.json describes\n"
    "(see README.md), by the adjoint method, and writes <prefix>.kernels,\n"
    "each line holding a radius, its depth, its weight in the integral over\n"
    "radius and the kernel of each parameter there, with their point parts at\n"
    "the source in the header, and <prefix>.measurement, the seismogram, the\n"
    "window and the adjoint source. The number of threads is set by\n"
    "OMP_NUM_THREADS.\n";

std::string Run(const std::filesystem::path& path) {
  using Program::Plural;
  const auto input = YSpec::ReadKernelInput(path);
  const auto& calculation = input.calculation;
  std::cerr << std::format(
      "yspec-kernels: model {}, degrees {} to {}, {} to {} mHz, traveltime "
      "on {} in {} to {} s, {}\n",
      input.model.string(), calculation.degrees.min, calculation.degrees.max,
      calculation.band.min * 1e3, calculation.band.max * 1e3, input.component,
      input.window.start, input.window.end,
      Plural(Program::Threads(), "thread"));
  const auto result =
      YSpec::ComputeKernels(input, Program::Reporter("yspec-kernels"));
  YSpec::WriteKernels(result, input);
  return std::format("kernels at {} radii with prefix {}",
                     result.kernels.Size(), input.prefix.string());
}

}  // namespace

int main(int argc, char** argv) {
  return Program::Main(argc, argv, "yspec-kernels", kUsage, Run);
}
