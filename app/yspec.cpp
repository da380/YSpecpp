// yspec: synthetic seismograms for a spherically symmetric earth model.
//
//   yspec input.json
//
// Reads the calculation from a JSON file (see YSpec/Input.hpp for its
// members) and writes one file per receiver, <prefix>.1, <prefix>.2, ..., in
// the order the receivers are listed. Each line of a file holds the time
// after the origin (s) and the vertical (up), north and east components, in
// SI units of the quantity asked for.

#include <YSpec/Input.hpp>
#include <YSpec/Synthetics.hpp>
#include <chrono>
#include <cstdio>
#include <exception>
#include <format>
#include <fstream>
#include <iostream>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

void Write(const YSpec::Seismograms& seismograms,
           const std::filesystem::path& prefix) {
  for (auto k = std::size_t{0}; k < seismograms.Receivers(); ++k) {
    auto path = prefix;
    path += std::format(".{}", k + 1);
    auto file = std::ofstream(path);
    if (!file) throw std::runtime_error("cannot write " + path.string());
    for (auto i = std::size_t{0}; i < seismograms.Samples(); ++i) {
      file << std::format("{:15.6E}{:15.6E}{:15.6E}{:15.6E}\n",
                          static_cast<double>(i) * seismograms.TimeStep(),
                          seismograms(k, 0, i), seismograms(k, 1, i),
                          seismograms(k, 2, i));
    }
  }
}

}  // namespace

int main(int argc, char** argv) {
  const auto usage =
      "usage: yspec input.json\n"
      "\n"
      "Computes synthetic seismograms for the calculation described in\n"
      "input.json (see README.md) and writes one file per receiver,\n"
      "<prefix>.1, <prefix>.2, ..., each line holding the time (s) and the\n"
      "vertical, north and east components in SI units. The number of\n"
      "threads is set by OMP_NUM_THREADS.\n";
  if (argc == 2 &&
      (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
    std::cout << usage;
    return 0;
  }
  if (argc == 2 && std::string(argv[1]) == "--version") {
    std::cout << "yspec " << YSPEC_VERSION << "\n";
    return 0;
  }
  if (argc != 2) {
    std::cerr << usage;
    return 2;
  }
  try {
    const auto start = std::chrono::steady_clock::now();
    const auto input = YSpec::ReadInput(argv[1]);
    const auto model = YSpec::Model::ReadDeck(input.model);
    const auto& s = input.settings;
    auto threads = 1;
#ifdef _OPENMP
    threads = omp_get_max_threads();
#endif
    std::cerr << std::format(
        "yspec: model {}, degrees {} to {}, {} to {} mHz, {} receivers, {} "
        "thread{}\n",
        input.model.string(), s.lMin, s.lMax, s.fMin, s.fMax,
        s.receivers.size(), threads, threads == 1 ? "" : "s");

    auto reported = 0;
    auto progress = [&reported](std::size_t done, std::size_t total) {
      const auto percent = static_cast<int>(100 * done / total);
      if (percent >= reported + 10 || done == total) {
        reported = percent - percent % 10;
        std::cerr << std::format("yspec: {:3}% of {} frequencies\n", percent,
                                 total);
      }
    };
    const auto seismograms = YSpec::ComputeSeismograms(model, s, progress);
    Write(seismograms, input.prefix);

    const auto seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
            .count();
    std::cerr << std::format("yspec: wrote {}.1 to {}.{} in {:.2f} s\n",
                             input.prefix.string(), input.prefix.string(),
                             seismograms.Receivers(), seconds);
  } catch (const std::exception& e) {
    std::cerr << "yspec: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
