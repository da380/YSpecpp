#include "TimeSeries.hpp"

#include <FFTWpp/Ranges>
#include <algorithm>
#include <cmath>
#include <ranges>
#include <stdexcept>

namespace YSpec {

std::vector<double> InverseTransform(const FrequencyGrid& grid,
                                     const std::optional<Taper>& taper,
                                     const std::vector<Complex>& spectra,
                                     std::size_t count) {
  using namespace FFTWpp;
  const auto samples = static_cast<int>(grid.Samples());
  // The frequencies of a real series of that length, zero to Nyquist.
  const auto bins = samples / 2 + 1;
  const auto frequencies = grid.Count();
  if (spectra.size() != count * frequencies) {
    throw std::invalid_argument("InverseTransform: spectra of the wrong size");
  }
  const auto howMany = static_cast<int>(count);
  const auto complexLayout =
      Ranges::Layout(1, std::vector{bins}, howMany, std::vector{bins}, 1, bins);
  const auto realLayout = Ranges::Layout(1, std::vector{samples}, howMany,
                                         std::vector{samples}, 1, samples);
  auto spectrum = FFTWpp::vector<Complex>(complexLayout.size());
  auto result = FFTWpp::vector<double>(realLayout.size());
  auto plan = Ranges::Plan(Ranges::View(spectrum, complexLayout),
                           Ranges::View(result, realLayout), Estimate);

  // Frequencies above the Nyquist frequency are dropped.
  const auto last = std::min<std::size_t>(grid.Last(), bins - 1);
  for (auto series : std::views::iota(std::size_t{0}, count)) {
    for (auto i = grid.First(); i <= last; ++i) {
      const auto f = static_cast<double>(i) * grid.Spacing();
      const auto weight = taper ? TaperValue(*taper, f) : 1.0;
      spectrum[series * bins + i] =
          weight * spectra[series * frequencies + (i - grid.First())];
    }
  }
  plan.Execute();

  // The sum approximates the inverse Fourier integral with the frequency
  // spacing, and the damping is undone.
  const auto timeStep = grid.TimeStep();
  const auto spacing = grid.Spacing();
  for (auto series : std::views::iota(std::size_t{0}, count)) {
    for (auto i : std::views::iota(0, samples)) {
      result[series * samples + i] *=
          spacing * std::exp(grid.Damping() * i * timeStep);
    }
  }
  return {result.begin(), result.end()};
}

}  // namespace YSpec
