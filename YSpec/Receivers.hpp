#pragma once

/**
 * @file Receivers.hpp
 * @brief Where the receivers are seen from the source, and the angular
 * functions that sum the response over degree and order.
 *
 * @details In a frame with the source at the north pole a receiver is at
 * colatitude @f$\Delta@f$, the epicentral distance, and longitude
 * @f$\phi = \pi - \zeta_s@f$ with @f$\zeta_s@f$ the azimuth of the receiver
 * from the source; the horizontal components there are rotated to north and
 * east with the azimuth of the source from the receiver. Geographic latitudes
 * are converted to geocentric with a flattening of about 1/297, and a point
 * exactly at a pole is moved 0.02 degrees off it, where the azimuths are
 * undefined.
 *
 * The angular functions are formed from the generalised Legendre functions
 * of GSHTrans: with @f$X^N_{lm}@f$ their orthonormalised values and
 * @f$\zeta = \sqrt{l(l+1)}@f$,
 *
 * @f[ \partial_\theta X_{lm} = \tfrac{\zeta}{2}(X^{-1}_{lm} - X^{+1}_{lm}),
 *     \qquad \frac{m}{\sin\theta} X_{lm} =
 *     -\tfrac{\zeta}{2}(X^{+1}_{lm} + X^{-1}_{lm}), @f]
 *
 * which are regular at both poles.
 */

#include <GSHTrans/Wigner.hpp>
#include <array>
#include <cmath>
#include <numbers>
#include <ranges>
#include <vector>

#include "Model.hpp"

namespace YSpec {

/** @brief A point on the surface, in geographic degrees. */
struct Location {
  double latitude;   ///< Geographic latitude, degrees north.
  double longitude;  ///< Longitude, degrees east.
};

/**
 * @brief A receiver as seen from the source.
 */
struct ReceiverGeometry {
  double distance;  ///< The epicentral distance (radians).
  double phi;       ///< The longitude in the frame with the source at the pole.
  double rotation;  ///< The angle that rotates the horizontal components.
};

/**
 * @brief The geometry of a receiver relative to a source.
 */
inline ReceiverGeometry Geometry(const Location& receiver,
                                 const Location& source) {
  constexpr auto pi = std::numbers::pi;
  constexpr auto degree = pi / 180;
  // The ratio of the tangents of geocentric and geographic latitude.
  constexpr auto geocentric = 0.993277;

  auto offPole = [](double latitude) {
    if (latitude == 90) return 89.98;
    if (latitude == -90) return -89.98;
    return latitude;
  };
  auto colatitude = [&](double latitude) {
    return pi / 2 -
           std::atan(geocentric * std::tan(offPole(latitude) * degree));
  };
  const auto thetaR = colatitude(receiver.latitude);
  const auto thetaS = colatitude(source.latitude);
  const auto phiR = receiver.longitude * degree;
  const auto phiS = source.longitude * degree;

  // The distance from the cross and dot products of the unit vectors, which
  // keeps its precision near 0 and pi where an arc cosine does not.
  auto unit = [](double theta, double phi) {
    return std::array{std::sin(theta) * std::cos(phi),
                      std::sin(theta) * std::sin(phi), std::cos(theta)};
  };
  const auto xR = unit(thetaR, phiR);
  const auto xS = unit(thetaS, phiS);
  const auto cross =
      std::hypot(xR[1] * xS[2] - xR[2] * xS[1], xR[2] * xS[0] - xR[0] * xS[2],
                 xR[0] * xS[1] - xR[1] * xS[0]);
  const auto dot = xR[0] * xS[0] + xR[1] * xS[1] + xR[2] * xS[2];
  const auto delta = std::atan2(cross, dot);

  // Where the two points coincide or are antipodal the azimuths are
  // undefined. Those chosen are the limits seen from a receiver approaching
  // along any great circle: at the source the receiver's azimuth plus pi is
  // the source's azimuth seen from the receiver, and at the antipode the two
  // azimuths are equal.
  if (cross == 0)
    return dot > 0 ? ReceiverGeometry{0, pi, pi}
                   : ReceiverGeometry{pi, pi, 2 * pi};

  // The azimuth at point a of point b, clockwise from north, in [0, 2 pi).
  auto azimuth = [](double thetaA, double phiA, double thetaB, double phiB) {
    const auto dPhi = phiB - phiA;
    const auto z =
        std::atan2(std::sin(dPhi) * std::sin(thetaB),
                   std::sin(thetaA) * std::cos(thetaB) -
                       std::cos(thetaA) * std::sin(thetaB) * std::cos(dPhi));
    return z < 0 ? z + 2 * pi : z;
  };
  const auto atReceiver = azimuth(thetaR, phiR, thetaS, phiS);
  const auto atSource = azimuth(thetaS, phiS, thetaR, phiR);
  return {delta, pi - atSource, 2 * pi - atReceiver};
}

/**
 * @brief The angular functions at one point, for one degree and order.
 */
struct Angular {
  double x;   ///< @f$X_{lm}(\theta)@f$.
  double dx;  ///< @f$\partial_\theta X_{lm}(\theta)@f$.
  double mx;  ///< @f$m X_{lm}(\theta) / \sin\theta@f$.
};

/**
 * @brief The angular functions at a set of colatitudes, for degrees 0 to
 * lMax and orders @f$|m| \le \min(l, 2)@f$ (zero for the others).
 */
class AngularFunctions {
 public:
  /**
   * @param lMax The largest degree.
   * @param theta The colatitudes, in @f$[0, \pi]@f$.
   */
  AngularFunctions(int lMax, std::vector<double> theta)
      : lMax_{lMax}, count_{std::ssize(theta)} {
    values_.assign(static_cast<std::size_t>(count_) * (lMax + 1) * 5, {});
    // GSHTrans refuses a double table above a degree at which a recursion for
    // high orders would underflow; long double reaches far higher.
    if (lMax <= GSHTrans::MaxSafeDegree<double>()) {
      Fill<double>(theta);
    } else {
      Fill<long double>(theta);
    }
  }

  /** @brief The largest degree. */
  int MaxDegree() const { return lMax_; }
  /** @brief The number of points. */
  Int Size() const { return count_; }

  /** @brief The functions at point @p i, degree @p l, order @p m. */
  const Angular& operator()(Int i, int l, int m) const {
    return values_[Index(i, l, m)];
  }

 private:
  int lMax_;
  Int count_;
  std::vector<Angular> values_;

  template <class Real>
  void Fill(const std::vector<double>& theta) {
    using Table = GSHTrans::Wigner<Real, GSHTrans::All, GSHTrans::All,
                                   GSHTrans::Multiple>;
    const auto nMax = lMax_ >= 1 ? 1 : 0;
    const auto angles = std::vector<Real>(theta.begin(), theta.end());
    const auto d = Table(lMax_, 2, nMax, angles);
    for (auto i : std::views::iota(Int{0}, count_)) {
      for (auto l : std::views::iota(0, lMax_ + 1)) {
        const auto zeta = std::sqrt(static_cast<double>(l) * (l + 1));
        const auto mMax = std::min(l, 2);
        for (auto m : std::views::iota(-mMax, mMax + 1)) {
          auto& a = Value(i, l, m);
          a.x = static_cast<double>(d[0, i][l][m]);
          if (l == 0) continue;
          const auto minus = static_cast<double>(d[-1, i][l][m]);
          const auto plus = static_cast<double>(d[1, i][l][m]);
          a.dx = zeta / 2 * (minus - plus);
          a.mx = -zeta / 2 * (plus + minus);
        }
      }
    }
  }

  // Degree outermost, so that a sum over the receivers at one degree reads
  // contiguous memory.
  std::size_t Index(Int i, int l, int m) const {
    return (static_cast<std::size_t>(l) * count_ + i) * 5 + (m + 2);
  }
  Angular& Value(Int i, int l, int m) { return values_[Index(i, l, m)]; }
};

}  // namespace YSpec
