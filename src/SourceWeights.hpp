#pragma once

/**
 * @file SourceWeights.hpp
 * @brief The jumps a moment-tensor point source imposes on the radial
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

#include <YSpec/Source.hpp>
#include <array>
#include <cmath>
#include <complex>
#include <numbers>

#include "Equations.hpp"
#include "OrderArray.hpp"

namespace YSpec {

/** @brief For each order, the weights of the two unit toroidal jumps. */
using ToroidalWeights = OrderArray<std::array<Complex, 2>>;

/**
 * @brief For each order, the weights of the four unit spheroidal jumps, in
 * @f$(U, V, P, S)@f$.
 */
using SpheroidalWeights = OrderArray<std::array<Complex, 4>>;

/** @brief The weights of the two unit radial jumps. */
using RadialWeights = std::array<Complex, 2>;

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
 * @param moment The moment tensor, non-dimensional.
 */
inline ToroidalWeights Toroidal(const Context& c, const ComplexMaterial& m,
                                const MomentTensor& moment) {
  using namespace std::complex_literals;
  const auto rs = m.r;
  // sqrt((l-1)(l+2)), the factor of the terms of order +-2.
  const auto rootZeta2m2 = std::sqrt(std::max(c.Zeta2() - 2, 0.0));
  const auto factor = SourceDetails::Factor(c, rs);
  auto w = ToroidalWeights{};
  w[-1][0] = 0.5 * rs * (-moment.rp + 1i * moment.rt) / m.L * factor;
  w[1][0] = 0.5 * rs * (moment.rp + 1i * moment.rt) / m.L * factor;
  w[-2][1] = 0.25 * rootZeta2m2 *
             (-1i * moment.tt + 1i * moment.pp + 2.0 * moment.tp) * factor;
  w[2][1] = 0.25 * rootZeta2m2 *
            (1i * moment.tt - 1i * moment.pp + 2.0 * moment.tp) * factor;
  return w;
}

/**
 * @brief The spheroidal weights for degree @f$l \ge 1@f$.
 */
inline SpheroidalWeights Spheroidal(const Context& c, const ComplexMaterial& m,
                                    const MomentTensor& moment) {
  using namespace std::complex_literals;
  const auto rs = m.r;
  const auto zeta = c.Zeta();
  // sqrt((l-1)(l+2)), the factor of the terms of order +-2.
  const auto rootZeta2m2 = std::sqrt(std::max(c.Zeta2() - 2, 0.0));
  const auto factor = SourceDetails::Factor(c, rs);
  const auto fOverC = m.F / m.C;
  auto w = SpheroidalWeights{};
  // The jumps are in (U, V, P, S).
  w[-2][3] =
      0.25 * rootZeta2m2 * (-moment.tt + moment.pp - 2.0i * moment.tp) * factor;
  w[-1][1] = 0.5 * rs * (moment.rt + 1i * moment.rp) / m.L * factor;
  w[0][0] = rs * moment.rr / m.C * factor;
  w[0][2] = (-moment.tt - moment.pp + 2.0 * moment.rr * fOverC) * factor;
  w[0][3] = (0.5 * zeta * (moment.tt + moment.pp) - zeta * moment.rr * fOverC) *
            factor;
  w[1][1] = 0.5 * rs * (-moment.rt + 1i * moment.rp) / m.L * factor;
  w[2][3] =
      0.25 * rootZeta2m2 * (-moment.tt + moment.pp + 2.0i * moment.tp) * factor;
  return w;
}

/**
 * @brief The radial (l = 0) weights.
 */
inline RadialWeights Radial(const Context& c, const ComplexMaterial& m,
                            const MomentTensor& moment) {
  const auto rs = m.r;
  const auto factor = SourceDetails::Factor(c, rs);
  return {rs * moment.rr / m.C * factor,
          (-moment.tt - moment.pp + 2.0 * moment.rr * m.F / m.C) * factor};
}

}  // namespace YSpec
