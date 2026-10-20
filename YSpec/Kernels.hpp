#pragma once

/**
 * @file Kernels.hpp
 * @brief Sensitivity kernels of measurements made on seismograms, by the
 * adjoint method.
 *
 * @details A measurement @f$J@f$ on one component of the seismogram of a
 * source at a receiver is linearised as @f$\delta J = \int j(t)\,\delta
 * u(t)\,dt@f$, with @f$j@f$ its adjoint source and @f$u@f$ the seismogram as
 * ComputeSeismograms gives it. Its sensitivity to the model is formed from two
 * fields at each of the complex frequencies of the calculation: the forward
 * field of the source, and the adjoint field, which is the response at the
 * same complex frequency to a force at the receiver along the component
 * measured, with a spectrum formed from @f$j@f$. No time reversal is needed.
 * The derivation, from the weak form of the equations of motion with the
 * equilibrium gravitational potential as a variable constrained by Poisson's
 * equation, is in docs/theory/adjoint-kernels.tex.
 *
 * The radial kernels give the measurement's sensitivity to spherically
 * symmetric perturbations of the model: of its density, at fixed moduli and
 * with the equilibrium gravity that follows from it, of its five Love moduli,
 * at fixed density, and of their inverse quality factors. They are evaluated
 * on a mesh in radius, with Gauss-Legendre nodes in elements between the
 * knots of the model, the source and the receiver, and come with the weights
 * of the quadrature on that mesh.
 */

#include <array>
#include <cstddef>
#include <span>
#include <vector>

#include "Calculation.hpp"
#include "Model.hpp"
#include "Receiver.hpp"
#include "Seismograms.hpp"
#include "Source.hpp"

namespace YSpec {

/**
 * @brief A window in time, with cosine tapers inside it at each end.
 */
struct Window {
  double start;      ///< Its start (s), after the origin time.
  double end;        ///< Its end (s).
  double taper = 0;  ///< The length of each taper (s).

  /** @brief The value of the window at a time (s). */
  double operator()(double t) const;
};

/**
 * @brief A measurement on one component of a seismogram, through its adjoint
 * source.
 *
 * @details The measurement changes by @f$\delta J = \int j(t)\,\delta u(t)\,
 * dt@f$, approximated by the sum over the samples of the seismogram, with
 * @f$u@f$ the seismogram as ComputeSeismograms gives it for the same
 * calculation, recorded at the receiver by its instrument, if any.
 */
struct Measurement {
  Receiver receiver;    ///< Where it is made.
  Component component;  ///< On which component.
  /// @f$j@f$ at the samples of the seismogram, in units of the measurement
  /// per unit of the seismogram per second.
  std::vector<double> adjointSource;
};

/**
 * @brief The traveltime of a phase in a window, measured by cross-correlation
 * with the synthetic itself.
 *
 * @details A delay @f$\delta\tau@f$ changes the seismogram by
 * @f$-\dot u\,\delta\tau@f$, so to first order @f$\delta\tau = -\int W\dot
 * u\,\delta u\,dt / \int W\dot u^2\,dt@f$ (Dahlen et al. 2000), and the adjoint
 * source is @f$j = -W\dot u / \int W\dot u^2\,dt@f$. The rate @f$\dot u@f$
 * is formed from the spectrum of the synthetic, and is exact.
 * @param model The earth model.
 * @param calculation The calculation of the seismogram.
 * @param source The source.
 * @param receiver The receiver.
 * @param component The component measured.
 * @param window The window.
 * @throws std::invalid_argument if the synthetic is zero in the window, or
 * as ComputeSeismograms does.
 */
Measurement Traveltime(const Model& model, const Calculation& calculation,
                       const Source& source, const Receiver& receiver,
                       Component component, const Window& window);

/**
 * @brief The parameters of a spherically symmetric model, as the radial
 * kernels give them.
 *
 * @details The inverse quality factors are those of the model: @c QA
 * attenuates A and C, @c QF attenuates F, and @c QMu attenuates L and N.
 */
enum class RadialParameter : std::size_t {
  Rho,  ///< Density, at fixed moduli.
  A,    ///< Love modulus A, at fixed density.
  C,    ///< Love modulus C.
  F,    ///< Love modulus F.
  L,    ///< Love modulus L.
  N,    ///< Love modulus N.
  QA,   ///< Inverse quality factor of A and C.
  QF,   ///< Inverse quality factor of F.
  QMu,  ///< Inverse quality factor of L and N.
  Count
};

/** @brief How finely the radial kernels are evaluated. */
struct RadialMesh {
  int order = 5;             ///< Gauss-Legendre nodes in each element.
  double maxElement = 50e3;  ///< The largest element (m).
};

/**
 * @brief Kernels for spherically symmetric perturbations of the model.
 *
 * @details For each parameter @f$X@f$,
 *
 * @f[ \delta J = \int K_X(r)\,\delta X(r)\,r^2\,dr
 *     + P_X\,\delta X(r_s), @f]
 *
 * the integral over the model, which the weights approximate,
 * @f$\sum_i w_i K_X(r_i)\,\delta X(r_i)@f$, and the second term a point part
 * at the radius of the source. A moment tensor makes the displacement jump
 * across that radius, by @f$\mathbf{M}\cdot\hat{\mathbf{r}}@f$ divided by
 * C or L there, so its strain has a part concentrated on it, and the kernels
 * of C, F and L, and of their quality factors, a part @f$P_X@f$ there; the
 * others have none.
 *
 * In a fluid only perturbations that keep it isotropic have meaning:
 * @f$\delta A = \delta C = \delta F@f$ and @f$\delta q_A = \delta q_F@f$, to
 * which the sum of the kernels responds. Those of L, N and their quality
 * factor are zero there.
 *
 * In SI units: a kernel is in units of the measurement per unit of its
 * parameter per cubic metre, a weight in cubic metres, and a point part in
 * units of the measurement per unit of its parameter.
 */
class RadialKernels {
 public:
  /**
   * @brief Zero kernels at the given radii (m), with their weights (m^3),
   * for a source at radius @p sourceRadius (m).
   */
  RadialKernels(std::vector<double> radii, std::vector<double> weights,
                double sourceRadius);

  /** @brief The number of radii. */
  std::size_t Size() const { return radii_.size(); }
  /** @brief The radii (m), from the centre. */
  std::span<const double> Radii() const { return radii_; }
  /** @brief The weights of the quadrature over radius (m^3). */
  std::span<const double> Weights() const { return weights_; }

  /** @brief The kernel of a parameter, at the radii. */
  std::span<double> operator[](RadialParameter p) {
    return values_[static_cast<std::size_t>(p)];
  }
  /** @brief The kernel of a parameter, at the radii. */
  std::span<const double> operator[](RadialParameter p) const {
    return values_[static_cast<std::size_t>(p)];
  }

  /** @brief The radius of the source (m). */
  double SourceRadius() const { return sourceRadius_; }
  /** @brief The point part of the kernel of a parameter at the source. */
  double& AtSource(RadialParameter p) {
    return atSource_[static_cast<std::size_t>(p)];
  }
  /** @brief The point part of the kernel of a parameter at the source. */
  double AtSource(RadialParameter p) const {
    return atSource_[static_cast<std::size_t>(p)];
  }

 private:
  static constexpr auto Parameters =
      static_cast<std::size_t>(RadialParameter::Count);
  std::vector<double> radii_;
  std::vector<double> weights_;
  std::array<std::vector<double>, Parameters> values_;
  double sourceRadius_;
  std::array<double, Parameters> atSource_{};
};

/**
 * @brief Computes the radial kernels of a measurement.
 * @param model The earth model.
 * @param calculation The calculation of the seismogram measured.
 * @param source The source.
 * @param measurement The measurement.
 * @param mesh How finely to evaluate the kernels.
 * @param progress Called as frequencies are completed, if given.
 * @throws std::invalid_argument if the adjoint source does not have a value
 * at every sample of the seismogram, if the receiver is not in the solid
 * region of the source (the adjoint field's source is there), if the
 * seismogram has the potential and tilt corrections, or as
 * ComputeSeismograms does.
 */
RadialKernels ComputeRadialKernels(const Model& model,
                                   const Calculation& calculation,
                                   const Source& source,
                                   const Measurement& measurement,
                                   const RadialMesh& mesh = {},
                                   const Progress& progress = {});

}  // namespace YSpec
