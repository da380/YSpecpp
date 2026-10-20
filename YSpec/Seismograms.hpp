#pragma once

/**
 * @file Seismograms.hpp
 * @brief Seismograms for every pair of a source and a receiver.
 */

#include <cstddef>
#include <span>
#include <vector>

namespace YSpec {

/** @brief The components of a seismogram. */
enum Component { Vertical = 0, North = 1, East = 2 };

/**
 * @brief Seismograms for a set of sources and receivers, in SI units of the
 * quantity recorded or, at a receiver with an instrument, in the units of its
 * output.
 *
 * @details Sample @f$i@f$ is at time @f$i\,\Delta t@f$ after the origin time.
 * Each trace is contiguous.
 */
class Seismograms {
 public:
  /** @brief Zero seismograms of a given size. */
  Seismograms(std::size_t sources, std::size_t receivers, std::size_t samples,
              double timeStep)
      : sources_{sources},
        receivers_{receivers},
        samples_{samples},
        timeStep_{timeStep},
        data_(sources * receivers * 3 * samples) {}

  /** @brief The number of sources. */
  std::size_t Sources() const { return sources_; }
  /** @brief The number of receivers. */
  std::size_t Receivers() const { return receivers_; }
  /** @brief The number of samples in each trace. */
  std::size_t Samples() const { return samples_; }
  /** @brief The interval between samples (s). */
  double TimeStep() const { return timeStep_; }

  /** @brief One component at one receiver from one source. */
  std::span<double> Trace(std::size_t source, std::size_t receiver,
                          Component component) {
    return {data_.data() + Offset(source, receiver, component), samples_};
  }
  /** @brief One component at one receiver from one source. */
  std::span<const double> Trace(std::size_t source, std::size_t receiver,
                                Component component) const {
    return {data_.data() + Offset(source, receiver, component), samples_};
  }

  /** @brief A sample. */
  double& operator()(std::size_t source, std::size_t receiver,
                     Component component, std::size_t i) {
    return data_[Offset(source, receiver, component) + i];
  }
  /** @brief A sample. */
  double operator()(std::size_t source, std::size_t receiver,
                    Component component, std::size_t i) const {
    return data_[Offset(source, receiver, component) + i];
  }

 private:
  std::size_t sources_;
  std::size_t receivers_;
  std::size_t samples_;
  double timeStep_;
  std::vector<double> data_;

  std::size_t Offset(std::size_t source, std::size_t receiver,
                     Component component) const {
    return ((source * receivers_ + receiver) * 3 + component) * samples_;
  }
};

}  // namespace YSpec
