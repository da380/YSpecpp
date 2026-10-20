// yspec: synthetic seismograms for a spherically symmetric earth model.
//
//   yspec input.json
//
// Reads the calculation from a JSON file (see YSpec/Input.hpp for its
// members), computes it and writes the seismograms (see YSpec/Run.hpp for the
// files), reporting progress on the standard error.

#include <YSpec/YSpec.hpp>
#include <chrono>
#include <exception>
#include <format>
#include <iostream>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

const auto kUsage =
    "usage: yspec input.json\n"
    "\n"
    "Computes synthetic seismograms for the calculation described in\n"
    "input.json (see README.md) and writes one file per receiver,\n"
    "<prefix>.1, <prefix>.2, ... (<prefix>.<source>.1, ... for several\n"
    "sources), each line holding the time (s) and the vertical, north and\n"
    "east components in SI units. The number of threads is set by\n"
    "OMP_NUM_THREADS.\n";

std::string Plural(std::size_t n, const char* noun) {
  return std::format("{} {}{}", n, noun, n == 1 ? "" : "s");
}

}  // namespace

int main(int argc, char** argv) {
  const auto argument = argc == 2 ? std::string(argv[1]) : std::string{};
  if (argument == "--help" || argument == "-h") {
    std::cout << kUsage;
    return 0;
  }
  if (argument == "--version") {
    std::cout << "yspec " << YSPEC_VERSION << "\n";
    return 0;
  }
  if (argc != 2) {
    std::cerr << kUsage;
    return 2;
  }
  try {
    const auto start = std::chrono::steady_clock::now();
    const auto input = YSpec::ReadInput(argument);
    const auto& calculation = input.calculation;
    auto threads = std::size_t{1};
#ifdef _OPENMP
    threads = static_cast<std::size_t>(omp_get_max_threads());
#endif
    std::cerr << std::format(
        "yspec: model {}, degrees {} to {}, {} to {} mHz, {}, {}, {}\n",
        input.model.string(), calculation.degrees.min, calculation.degrees.max,
        calculation.band.min * 1e3, calculation.band.max * 1e3,
        Plural(input.sources.size(), "source"),
        Plural(input.receivers.size(), "receiver"), Plural(threads, "thread"));

    auto reported = 0;
    auto progress = [&reported](std::size_t done, std::size_t total) {
      const auto percent = static_cast<int>(100 * done / total);
      if (percent >= reported + 10 || done == total) {
        reported = percent - percent % 10;
        std::cerr << std::format("yspec: {:3}% of {} frequencies\n", percent,
                                 total);
      }
    };
    const auto seismograms = YSpec::Compute(input, progress);
    YSpec::WriteSeismograms(seismograms, input);

    const auto seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    std::cerr << std::format(
        "yspec: wrote {} with prefix {} in {:.2f} s\n",
        Plural(seismograms.Sources() * seismograms.Receivers(), "file"),
        input.prefix.string(), seconds);
  } catch (const std::exception& e) {
    std::cerr << "yspec: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
