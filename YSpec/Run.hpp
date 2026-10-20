#pragma once

/**
 * @file Run.hpp
 * @brief Calculations described by input files, from end to end: what the
 * programs do.
 */

#include <filesystem>

#include "Input.hpp"
#include "Kernels.hpp"
#include "Seismograms.hpp"
#include "Snapshots.hpp"
#include "Synthetics.hpp"

namespace YSpec {

/**
 * @brief Reads the model an input names and computes its seismograms.
 * @param input The calculation.
 * @param progress Called as frequencies are completed, if given.
 */
Seismograms ComputeSeismograms(const SeismogramInput& input,
                               const Progress& progress = {});

/**
 * @brief Writes seismograms as text, one file per source and receiver,
 * named after the input's output prefix: `prefix.1`, `prefix.2`, ... in the
 * order of the receivers when there is one source, and `prefix.name.1`, ...,
 * with the source's name, when there are several. Each line holds the time
 * after the origin (s) and the vertical, north and east components.
 * @param seismograms The seismograms.
 * @param input The calculation they are for, which gives the prefix and the
 * names of the sources.
 * @throws std::runtime_error if a file cannot be written.
 */
void WriteSeismograms(const Seismograms& seismograms,
                      const SeismogramInput& input);

/**
 * @brief Reads an input file, computes the seismograms it describes and
 * writes them.
 * @param path The input file.
 * @param progress Called as frequencies are completed, if given.
 * @return The seismograms.
 */
Seismograms RunSeismograms(const std::filesystem::path& path,
                           const Progress& progress = {});

/**
 * @brief Reads the model an input names and computes its snapshots.
 * @param input The calculation.
 * @param progress Called as frequencies are completed, if given.
 */
Snapshots ComputeSnapshots(const SnapshotInput& input,
                           const Progress& progress = {});

/**
 * @brief Writes snapshots as text, one file per source and time, named
 * `prefix.1`, `prefix.2`, ... in the order of the times when there is one
 * source, and `prefix.name.1`, ..., with the source's name, when there are
 * several. After two comment lines starting with '#', each line holds a point
 * of the grid, depth by depth, azimuth by azimuth and distance by distance:
 * its depth (km), distance (degrees) and azimuth (degrees) from the source,
 * its latitude and longitude, and the up, away and transverse components (see
 * SnapshotComponent).
 * @param snapshots The snapshots.
 * @param input The calculation they are for.
 * @throws std::runtime_error if a file cannot be written.
 */
void WriteSnapshots(const Snapshots& snapshots, const SnapshotInput& input);

/**
 * @brief Reads an input file, computes the snapshots it describes and writes
 * them.
 * @param path The input file.
 * @param progress Called as frequencies are completed, if given.
 * @return The snapshots.
 */
Snapshots RunSnapshots(const std::filesystem::path& path,
                       const Progress& progress = {});

/** @brief What the yspec-kernels program computes. */
struct KernelResult {
  Seismograms seismogram;   ///< The seismogram at the receiver.
  Measurement measurement;  ///< The traveltime, by its adjoint source.
  RadialKernels kernels;    ///< Its kernels in the model's parameters.
  NamedKernels named;       ///< And in those asked for.
};

/**
 * @brief Reads the model an input names and computes the kernels of its
 * traveltime.
 * @param input The calculation.
 * @param progress Called as frequencies of the kernels are completed.
 */
KernelResult ComputeKernels(const KernelInput& input,
                            const Progress& progress = {});

/**
 * @brief Writes the kernels and the measurement as text: `prefix.kernels`,
 * whose header gives the measurement, the units and the point parts at the
 * source (on a line beginning `# at_source`), and whose lines give the
 * radius and depth (km), the weight (m^3) and each kernel, in the order of
 * the column names on the last line of the header; and `prefix.measurement`,
 * whose lines give the time (s), the seismogram measured, the window and the
 * adjoint source.
 * @throws std::runtime_error if a file cannot be written.
 */
void WriteKernels(const KernelResult& result, const KernelInput& input);

/**
 * @brief Reads an input file, computes the kernels it describes and writes
 * them.
 */
KernelResult RunKernels(const std::filesystem::path& path,
                        const Progress& progress = {});

}  // namespace YSpec
