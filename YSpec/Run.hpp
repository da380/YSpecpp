#pragma once

/**
 * @file Run.hpp
 * @brief A calculation described by an input file, from end to end: what the
 * yspec program does.
 */

#include <filesystem>

#include "Input.hpp"
#include "Seismograms.hpp"
#include "Synthetics.hpp"

namespace YSpec {

/**
 * @brief Reads the model an input names and computes its seismograms.
 * @param input The calculation.
 * @param progress Called as frequencies are completed, if given.
 */
Seismograms Compute(const Input& input, const Progress& progress = {});

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
void WriteSeismograms(const Seismograms& seismograms, const Input& input);

/**
 * @brief Reads an input file, computes the seismograms it describes and
 * writes them.
 * @param path The input file.
 * @param progress Called as frequencies are completed, if given.
 * @return The seismograms.
 */
Seismograms Run(const std::filesystem::path& path,
                const Progress& progress = {});

}  // namespace YSpec
