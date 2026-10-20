#pragma once

/**
 * @file Source.hpp
 * @brief A moment-tensor point source: the jumps it imposes on the radial
 * solutions.
 *
 * @details In a frame with the source at the north pole only orders
 * @f$|m| \le 2@f$ are excited. For each degree and order the source imposes a
 * jump on the displacement and traction at the source radius; the radial
 * solvers give the response to unit jumps, and the response to the source is
 * their combination with the weights computed here. The time dependence is a
 * step in moment, hence the factor @f$1/(i\omega)@f$. The formulas are those
 * of the original YSpec.
 */

#include <array>
#include <cmath>
#include <complex>
#include <numbers>

#include "Equations.hpp"
#include "Units.hpp"

namespace YSpec {

/**
 * @brief The components of a moment tensor in spherical coordinates
 * @f$(r, \theta, \phi)@f$ at the source, @f$\theta@f$ being colatitude.
 */
struct MomentTensor {
  double rr = 0;  ///< @f$M_{rr}@f$.
  double rt = 0;  ///< @f$M_{r\theta}@f$.
  double rp = 0;  ///< @f$M_{r\phi}@f$.
  double tt = 0;  ///< @f$M_{\theta\theta}@f$.
  double tp = 0;  ///< @f$M_{\theta\phi}@f$.
  double pp = 0;  ///< @f$M_{\phi\phi}@f$.

  /** @brief The tensor divided by a scale. */
  MomentTensor operator/(double scale) const {
    return {rr / scale, rt / scale, rp / scale,
            tt / scale, tp / scale, pp / scale};
  }
};

/**
 * @brief For each order @f$m = -2, \dots, 2@f$ (index @f$m + 2@f$), the
 * weights of the responses to the two unit toroidal jumps.
 */
using ToroidalSource = std::array<std::array<Complex, 2>, 5>;

/**
 * @brief For each order @f$m = -2, \dots, 2@f$, the weights of the responses
 * to the four unit spheroidal jumps in @f$(U, V, P, S)@f$.
 */
using SpheroidalSource = std::array<std::array<Complex, 4>, 5>;

/** @brief The weights of the responses to the two unit radial jumps. */
using RadialSource = std::array<Complex, 2>;

namespace SourceDetails {

// sqrt((2l+1)/4pi) / (i omega r^2), the factor common to every jump.
inline Complex Factor(const Context& c, double r) {
  const auto nu = std::sqrt((2.0 * c.Degree() + 1) / (4 * std::numbers::pi));
  return nu / (Complex(0, 1) * c.Omega() * r * r);
}

}  // namespace SourceDetails

/**
 * @brief The toroidal weights for degree @f$l \ge 1@f$.
 * @param c The degree and frequency.
 * @param m The material at the source.
 * @param M The moment tensor, non-dimensional.
 */
inline ToroidalSource Toroidal(const Context& c, const ComplexMaterial& m,
                               const MomentTensor& M) {
  using namespace std::complex_literals;
  const auto rs = m.r;
  const auto s = std::sqrt(std::max(c.Zeta2() - 2, 0.0));
  const auto factor = SourceDetails::Factor(c, rs);
  auto w = ToroidalSource{};
  w[1][0] = 0.5 * rs * (-M.rp + 1i * M.rt) / m.L * factor;
  w[3][0] = 0.5 * rs * (M.rp + 1i * M.rt) / m.L * factor;
  w[0][1] = 0.25 * s * (-1i * M.tt + 1i * M.pp + 2.0 * M.tp) * factor;
  w[4][1] = 0.25 * s * (1i * M.tt - 1i * M.pp + 2.0 * M.tp) * factor;
  return w;
}

/**
 * @brief The spheroidal weights for degree @f$l \ge 1@f$.
 */
inline SpheroidalSource Spheroidal(const Context& c, const ComplexMaterial& m,
                                   const MomentTensor& M) {
  using namespace std::complex_literals;
  const auto rs = m.r;
  const auto zeta = c.Zeta();
  const auto s = std::sqrt(std::max(c.Zeta2() - 2, 0.0));
  const auto factor = SourceDetails::Factor(c, rs);
  const auto fc = m.F / m.C;
  auto w = SpheroidalSource{};
  // The jumps are in (U, V, P, S).
  w[0][3] = 0.25 * s * (-M.tt + M.pp - 2.0i * M.tp) * factor;
  w[1][1] = 0.5 * rs * (M.rt + 1i * M.rp) / m.L * factor;
  w[2][0] = rs * M.rr / m.C * factor;
  w[2][2] = (-M.tt - M.pp + 2.0 * M.rr * fc) * factor;
  w[2][3] = (0.5 * zeta * (M.tt + M.pp) - zeta * M.rr * fc) * factor;
  w[3][1] = 0.5 * rs * (-M.rt + 1i * M.rp) / m.L * factor;
  w[4][3] = 0.25 * s * (-M.tt + M.pp + 2.0i * M.tp) * factor;
  return w;
}

/**
 * @brief The radial (l = 0) weights.
 */
inline RadialSource Radial(const Context& c, const ComplexMaterial& m,
                           const MomentTensor& M) {
  const auto rs = m.r;
  const auto factor = SourceDetails::Factor(c, rs);
  return {rs * M.rr / m.C * factor,
          (-M.tt - M.pp + 2.0 * M.rr * m.F / m.C) * factor};
}

}  // namespace YSpec
