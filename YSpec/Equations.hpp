#pragma once

/**
 * @file Equations.hpp
 * @brief The radial equations for one degree and one complex frequency.
 *
 * @details Each system has the form
 *
 * @f[ \frac{d}{dr}\begin{pmatrix} u \\ t \end{pmatrix} =
 *     \begin{pmatrix} T & K \\ S & -T^{\mathsf{T}} \end{pmatrix}
 *     \begin{pmatrix} u \\ t \end{pmatrix}, @f]
 *
 * with @f$K@f$ and @f$S@f$ symmetric, and is specified by the entries of the
 * blocks (Block2, Block4, Block6). The spheroidal systems are integrated as
 * minor vectors: second-order minors (five independent components) for the
 * 4x4 systems and third-order minors (fourteen, after the identities of a
 * solution regular at the centre) for the 6x6 system. The source is carried to
 * the receiver by the "b-vector" system. The orderings and the formulas are
 * those of the original YSpec (Al-Attar & Woodhouse 2008), transcribed
 * unchanged.
 *
 * The variables are scaled by the radius, as in the original code. With
 * @f$\zeta = \sqrt{l(l+1)}@f$ and @f$U, V, W@f$ the radial, consoidal and
 * toroidal displacement coefficients:
 *
 * - toroidal: @f$y_1 = \zeta r W@f$, and @f$y_2@f$ the matching traction;
 * - spheroidal, 6x6: @f$(y_1, \dots, y_6) = (rU, \zeta r V, r\phi, P, S,
 *   \psi)@f$, with @f$\phi@f$ the perturbation of the potential, @f$P@f$
 *   and @f$S@f$ the matching tractions, and the flux
 *   @f$\psi = r(\phi' + (l+1)\phi/r + 4\pi G\rho U)/(4\pi G)@f$, which
 *   is continuous and vanishes at a free surface;
 * - spheroidal, 4x4: @f$(rU, \zeta r V, P, S)@f$;
 * - radial: @f$y_1 = rU@f$, and @f$y_2@f$ the matching traction.
 *
 * Only the displacements, the potential and the flux at the receiver are
 * used.
 */

#include <array>
#include <cmath>
#include <complex>
#include <numbers>

#include "Bessel.hpp"
#include "Model.hpp"
#include "RungeKutta.hpp"
#include "Units.hpp"

namespace YSpec {

/** @brief How gravity enters the equations. */
enum class Gravity {
  None,     ///< No gravity at all.
  Cowling,  ///< The equilibrium gravity, but no perturbation of the potential.
  Full,     ///< Full self-gravitation.
};

/**
 * @brief The material at a radius at complex frequency: the moduli made
 * complex by the dispersion of a constant Q.
 */
struct ComplexMaterial {
  double r;       ///< The radius.
  double rho;     ///< Density.
  double g;       ///< Gravitational acceleration.
  Complex A;      ///< Love modulus A.
  Complex C;      ///< Love modulus C.
  Complex F;      ///< Love modulus F.
  Complex L;      ///< Love modulus L.
  Complex N;      ///< Love modulus N.
  Complex kappa;  ///< Voigt average bulk modulus.
  Complex mu;     ///< Voigt average shear modulus.
};

/**
 * @brief The degree and frequency of one solve, with the quantities derived
 * from them.
 *
 * @details The moduli at complex frequency @f$\omega@f$ are those of a constant
 * Q with logarithmic dispersion,
 *
 * @f[ M(\omega) = M_0\left[1 + q\,\frac{2}{\pi}
 *     \ln\frac{i\omega}{\omega_0}\right], @f]
 *
 * which is exact rather than first order in @f$q = 1/Q@f$.
 */
class Context {
 public:
  /**
   * @param model The earth model.
   * @param l The degree.
   * @param omega The complex angular frequency, non-dimensional, with a
   * negative imaginary part.
   * @param attenuation Whether to include attenuation.
   */
  Context(const Model& model, int l, Complex omega, bool attenuation = true)
      : model_{&model},
        l_{l},
        zeta2_{static_cast<double>(l) * (l + 1)},
        zeta_{std::sqrt(zeta2_)},
        omega_{omega},
        omega2_{omega * omega},
        omegaInverse2_{1.0 / omega2_},
        dispersion_{attenuation ? 2 / std::numbers::pi *
                                      std::log(Complex(0, 1) * omega /
                                               model.ReferenceFrequency())
                                : Complex{}} {}

  /** @brief The model. */
  const Model& GetModel() const { return *model_; }
  /** @brief The degree @f$l@f$. */
  int Degree() const { return l_; }
  /** @brief @f$l(l+1)@f$. */
  double Zeta2() const { return zeta2_; }
  /** @brief @f$\sqrt{l(l+1)}@f$. */
  double Zeta() const { return zeta_; }
  /** @brief The complex frequency. */
  Complex Omega() const { return omega_; }
  /** @brief The square of the frequency. */
  Complex Omega2() const { return omega2_; }
  /** @brief The inverse square of the frequency. */
  Complex OmegaInverse2() const { return omegaInverse2_; }
  /** @brief @f$(2/\pi)\ln(i\omega/\omega_0)@f$, or zero without attenuation. */
  Complex Dispersion() const { return dispersion_; }

  /** @brief The material at a radius in a layer, at this frequency. */
  ComplexMaterial operator()(Int layer, double r) const {
    return Complexify((*model_)(layer, r), r);
  }

  /** @brief A real material made complex at this frequency. */
  ComplexMaterial Complexify(const Material& m, double r) const {
    const auto d = dispersion_;
    return {.r = r,
            .rho = m.rho,
            .g = m.g,
            .A = m.A * (1.0 + m.qA * d),
            .C = m.C * (1.0 + m.qA * d),
            .F = m.F * (1.0 + m.qF * d),
            .L = m.L * (1.0 + m.qMu * d),
            .N = m.N * (1.0 + m.qMu * d),
            .kappa = m.kappa * (1.0 + m.qKappa * d),
            .mu = m.mu * (1.0 + m.qMu * d)};
  }

 private:
  const Model* model_;
  int l_;
  double zeta2_;
  double zeta_;
  Complex omega_;
  Complex omega2_;
  Complex omegaInverse2_;
  Complex dispersion_;
};

//----------------------------------------------------------------------------
// The blocks of the coefficient matrices.
//----------------------------------------------------------------------------

/** @brief A 2x2 system: @f$T = t@f$, @f$K = k@f$, @f$S = s@f$. */
struct Block2 {
  Complex t;  ///< T.
  Complex k;  ///< K.
  Complex s;  ///< S.
};

/**
 * @brief A 4x4 system, with the entries of T, K and S that are not zero.
 * Those that are real in both the solid and the fluid forms are held as
 * real, which halves the cost of the products they enter.
 */
struct Block4 {
  Complex t11;  ///< T(1,1).
  Complex t12;  ///< T(1,2).
  double t21;   ///< T(2,1).
  double t22;   ///< T(2,2).
  Complex k11;  ///< K(1,1).
  Complex k22;  ///< K(2,2).
  Complex s11;  ///< S(1,1).
  Complex s12;  ///< S(1,2) = S(2,1).
  Complex s22;  ///< S(2,2).
};

/** @brief The 6x6 self-gravitating system in a solid. */
struct Block6 {
  Complex t11;  ///< T(1,1).
  Complex t12;  ///< T(1,2).
  double t21;   ///< T(2,1).
  double t22;   ///< T(2,2).
  double t31;   ///< T(3,1).
  double t33;   ///< T(3,3).
  Complex k11;  ///< K(1,1).
  Complex k22;  ///< K(2,2).
  double k33;   ///< K(3,3).
  Complex s11;  ///< S(1,1).
  Complex s12;  ///< S(1,2) = S(2,1).
  double s13;   ///< S(1,3) = S(3,1).
  Complex s22;  ///< S(2,2).
  double s23;   ///< S(2,3) = S(3,2).
};

/** @brief The toroidal system. */
inline Block2 ToroidalBlock(const Context& c, const ComplexMaterial& m) {
  const auto ri = 1 / m.r;
  return {.t = 2 * ri,
          .k = 1.0 / m.L,
          .s = -c.Omega2() * m.rho + ri * ri * (c.Zeta2() - 2) * m.N};
}

/**
 * @brief The radial (l = 0) system. With gravity, the perturbation of the
 * potential has been eliminated: for l = 0 its gradient cancels the density
 * term of the gradient of gravity.
 */
inline Block2 RadialBlock(const Context& c, const ComplexMaterial& m,
                          bool gravity) {
  const auto ri = 1 / m.r;
  const auto cInverse = 1.0 / m.C;
  const auto gamma = m.A - m.N - m.F * m.F * cInverse;
  const auto g = gravity ? m.g : 0.0;
  return {.t = (1.0 - 2.0 * m.F * cInverse) * ri,
          .k = cInverse,
          .s = -c.Omega2() * m.rho + 4.0 * (gamma - m.rho * g * m.r) * ri * ri};
}

/** @brief The spheroidal system in a solid, with full self-gravitation. */
inline Block6 SolidBlock6(const Context& c, const ComplexMaterial& m) {
  const auto ri = 1 / m.r;
  const auto zeta = c.Zeta();
  const auto zeta2 = c.Zeta2();
  const auto l = static_cast<double>(c.Degree());
  const auto cInverse = 1.0 / m.C;
  const auto gamma = m.A - m.N - m.F * m.F * cInverse;
  const auto rho = m.rho;
  const auto gg = m.g;
  return {
      .t11 = (1.0 - 2.0 * m.F * cInverse) * ri,
      .t12 = zeta * ri * m.F * cInverse,
      .t21 = -zeta * ri,
      .t22 = 2 * ri,
      .t31 = -4 * PiG * rho,
      .t33 = -l * ri,
      .k11 = cInverse,
      .k22 = 1.0 / m.L,
      .k33 = 4 * PiG,
      .s11 = -c.Omega2() * rho + 4.0 * (gamma - rho * gg * m.r) * ri * ri,
      .s12 = zeta * (rho * gg * m.r - 2.0 * gamma) * ri * ri,
      .s13 = -rho * (l + 1) * ri,
      .s22 = -c.Omega2() * rho + (zeta2 * (gamma + m.N) - 2.0 * m.N) * ri * ri,
      .s23 = rho * zeta * ri};
}

/**
 * @brief The spheroidal system in a solid, without gravity or in the Cowling
 * approximation.
 */
inline Block4 SolidBlock4(const Context& c, const ComplexMaterial& m,
                          bool cowling) {
  const auto ri = 1 / m.r;
  const auto zeta = c.Zeta();
  const auto zeta2 = c.Zeta2();
  const auto cInverse = 1.0 / m.C;
  const auto gamma = m.A - m.N - m.F * m.F * cInverse;
  const auto rho = m.rho;
  const auto rg = cowling ? rho * m.g * m.r : 0.0;
  return {
      .t11 = (1.0 - 2.0 * m.F * cInverse) * ri,
      .t12 = zeta * ri * m.F * cInverse,
      .t21 = -zeta * ri,
      .t22 = 2 * ri,
      .k11 = cInverse,
      .k22 = 1.0 / m.L,
      .s11 = -c.Omega2() * rho + 4.0 * (gamma - rg) * ri * ri,
      .s12 = zeta * (rg - 2.0 * gamma) * ri * ri,
      .s22 = -c.Omega2() * rho + (zeta2 * (gamma + m.N) - 2.0 * m.N) * ri * ri};
}

/** @brief The spheroidal system in a fluid, with full self-gravitation. */
inline Block4 FluidBlock4(const Context& c, const ComplexMaterial& m) {
  const auto ri = 1 / m.r;
  const auto zeta2 = c.Zeta2();
  const auto l = static_cast<double>(c.Degree());
  const auto wi2 = c.OmegaInverse2();
  const auto rho = m.rho;
  const auto gg = m.g;
  return {
      .t11 = (gg * zeta2 * wi2 * ri - 1.0) * ri,
      .t12 = zeta2 * wi2 * ri * ri,
      .t21 = -4 * PiG * rho,
      .t22 = -l * ri,
      .k11 = 1.0 / m.C - zeta2 * wi2 * ri * ri / rho,
      .k22 = 4 * PiG,
      .s11 = -c.Omega2() * rho + rho * gg * (gg * zeta2 * ri * wi2 - 4.0) * ri,
      .s12 = rho * (gg * zeta2 * ri * wi2 - (l + 1)) * ri,
      .s22 = rho * zeta2 * ri * ri * wi2};
}

/**
 * @brief The spheroidal system in a fluid, without gravity or in the Cowling
 * approximation.
 */
inline Block2 FluidBlock2(const Context& c, const ComplexMaterial& m,
                          bool cowling) {
  const auto ri = 1 / m.r;
  const auto zeta2 = c.Zeta2();
  const auto wi2 = c.OmegaInverse2();
  const auto rho = m.rho;
  const auto gg = cowling ? m.g : 0.0;
  return {
      .t = (gg * zeta2 * wi2 * ri - 1.0) * ri,
      .k = 1.0 / m.C - zeta2 * wi2 * ri * ri / rho,
      .s = -c.Omega2() * rho + rho * gg * (gg * zeta2 * ri * wi2 - 4.0) * ri};
}

//----------------------------------------------------------------------------
// Derivatives.
//----------------------------------------------------------------------------

/** @brief The derivative of the solution of a 2x2 system. */
inline Vector<2> Derivative(const Block2& c, const Vector<2>& y) {
  return {c.t * y[0] + c.k * y[1], c.s * y[0] - c.t * y[1]};
}

/** @brief The derivative of the second-order minors of a 4x4 system. */
inline Vector<5> MinorDerivative(const Block4& c, const Vector<5>& m) {
  const auto& [t11, t12, t21, t22, k11, k22, s11, s12, s22] = c;
  const auto a11 = t11 + t22;
  const auto a22 = t11 - t22;
  auto dmdr = Vector<5>{};
  dmdr[0] = a11 * m[0] + k22 * m[2] - k11 * m[3];
  dmdr[1] = s12 * m[0] - t21 * m[2] + t12 * m[3];
  dmdr[2] = s22 * m[0] - 2.0 * t12 * m[1] + a22 * m[2] + k11 * m[4];
  dmdr[3] = -s11 * m[0] + 2.0 * t21 * m[1] - a22 * m[3] - k22 * m[4];
  dmdr[4] = -2.0 * s12 * m[1] + s11 * m[2] - s22 * m[3] - a11 * m[4];
  return dmdr;
}

/** @brief The derivative of the reduced third-order minors of the 6x6 system.
 */
inline Vector<14> MinorDerivative(const Block6& c, const Vector<14>& m) {
  const auto& [t11, t12, t21, t22, t31, t33, k11, k22, k33, s11, s12, s13, s22,
               s23] = c;
  const auto a123 = t11 + t22 + t33;
  const auto a126 = t11 + t22 - t33;
  const auto a135 = t11 - t22 + t33;
  const auto a234 = -t11 + t22 + t33;
  auto dmdr = Vector<14>{};
  dmdr[0] = k33 * m[3] - k22 * m[5] + k11 * m[9] + m[0] * a123;
  dmdr[1] = -k22 * m[6] + m[0] * s13 - m[2] * t21 + m[1] * t22 - m[3] * t31;
  dmdr[2] = -k11 * m[10] + m[0] * s23 + m[2] * t11 - m[1] * t12;
  dmdr[3] = k22 * m[8] - k11 * m[11] + m[3] * a126;
  dmdr[4] = -k33 * m[7] - m[0] * s12 + m[9] * t12 - m[5] * t21 + m[2] * t31 +
            m[4] * t33;
  dmdr[5] =
      -k33 * m[8] - k11 * m[12] - m[0] * s22 - 2.0 * m[4] * t12 + m[5] * a135;
  dmdr[6] = m[2] * s12 + m[5] * s13 - m[1] * s22 - m[4] * s23 + m[10] * t12 -
            m[6] * t22 + m[8] * t31;
  dmdr[7] = m[3] * s12 - m[2] * s13 - m[1] * s23 + m[11] * t12 - m[8] * t21 -
            m[7] * t33;
  dmdr[8] = k11 * m[13] + m[3] * s22 - 2.0 * m[2] * s23 - 2.0 * m[7] * t12 -
            m[8] * a234;
  dmdr[9] = -k33 * m[11] + k22 * m[12] + m[0] * s11 + 2.0 * m[4] * t21 -
            2.0 * m[1] * t31 + m[9] * a234;
  dmdr[10] = -m[2] * s11 + m[1] * s12 - m[4] * s13 - m[9] * s23 - m[10] * t11 +
             m[6] * t21 - m[7] * t31;
  dmdr[11] = -k22 * m[13] - m[3] * s11 + 2.0 * m[1] * s13 + 2.0 * m[7] * t21 -
             m[11] * a135;
  dmdr[12] = k33 * m[13] - m[5] * s11 + 2.0 * m[4] * s12 + m[9] * s22 +
             2.0 * m[6] * t31 - m[12] * a126;
  dmdr[13] = m[8] * s11 - 2.0 * m[7] * s12 + 2.0 * m[6] * s13 - m[11] * s22 +
             2.0 * m[10] * s23 - m[13] * a123;
  return dmdr;
}

/** @brief The derivative of the b-vector of a 4x4 system. */
inline Vector<4> BDerivative(const Block4& c, const Vector<4>& b) {
  const auto& [t11, t12, t21, t22, k11, k22, s11, s12, s22] = c;
  auto dbdr = Vector<4>{};
  dbdr[0] = -s11 * b[2] - s12 * b[3] - t11 * b[0] - t21 * b[1];
  dbdr[1] = -s12 * b[2] - s22 * b[3] - t12 * b[0] - t22 * b[1];
  dbdr[2] = -k11 * b[0] + t11 * b[2] + t12 * b[3];
  dbdr[3] = -k22 * b[1] + t21 * b[2] + t22 * b[3];
  return dbdr;
}

/** @brief The derivative of the b-vector of the 6x6 system. */
inline Vector<14> BDerivative(const Block6& c, const Vector<14>& b) {
  const auto& [t11, t12, t21, t22, t31, t33, k11, k22, k33, s11, s12, s13, s22,
               s23] = c;
  auto dbdr = Vector<14>{};
  dbdr[0] = s11 * b[6] + s12 * (b[7] - b[2]) + s13 * b[8] - s22 * b[3] -
            s23 * b[4] - (t11 + t22) * b[0] + t31 * b[5];
  dbdr[1] = s11 * b[9] + s12 * b[10] - s13 * (2.0 * b[2] + b[7]) - s23 * b[3] -
            t21 * b[5] - (t11 + t33) * b[1];
  dbdr[2] = s12 * b[11] + s13 * b[12] + t12 * b[3] - t21 * b[6] - t31 * b[9];
  dbdr[3] = -k22 * b[0] - s11 * b[11] + s13 * b[13] + t21 * (b[2] - b[7]) +
            (t22 - t11) * b[3] - t31 * b[10];
  dbdr[4] = -k33 * b[1] - s11 * b[12] - s12 * b[13] - t21 * b[8] +
            t31 * (2.0 * b[2] + b[7]) + (t33 - t11) * b[4];
  dbdr[5] = s12 * b[9] - s13 * b[6] + s22 * b[10] - s23 * (b[2] + 2.0 * b[7]) -
            t12 * b[1] - (t22 + t33) * b[5];
  dbdr[6] = k11 * b[0] + s22 * b[11] + s23 * b[12] - t12 * (b[2] - b[7]) +
            (t11 - t22) * b[6];
  dbdr[7] = -s12 * b[11] + s23 * b[13] - t12 * b[3] + t21 * b[6];
  dbdr[8] = -k33 * b[5] - s12 * b[12] - s22 * b[13] - t12 * b[4] - t22 * b[8] +
            t31 * b[6] + t33 * b[8];
  dbdr[9] = k11 * b[1] + s23 * b[11] + t12 * b[10] + (t11 - t33) * b[9];
  dbdr[10] = k22 * b[5] - s13 * b[11] + t21 * b[9] + (t22 - t33) * b[10];
  dbdr[11] = -k11 * b[3] + k22 * b[6] + (t11 + t22) * b[11];
  dbdr[12] = -k11 * b[4] + k33 * b[9] + t12 * b[13] + (t11 + t33) * b[12];
  dbdr[13] = -k22 * b[8] + k33 * b[10] + t21 * b[12] - t31 * b[11] +
             (t22 + t33) * b[13];
  return dbdr;
}

//----------------------------------------------------------------------------
// Minors of starting solutions.
//----------------------------------------------------------------------------

namespace EquationDetails {

inline Complex Det2(Complex a11, Complex a12, Complex a21, Complex a22) {
  return a11 * a22 - a12 * a21;
}

// The minor of rows i and j of the 4x2 matrix whose columns are y1 and y2.
inline Complex Minor2(const Vector<4>& y1, const Vector<4>& y2, int i, int j) {
  return Det2(y1[i], y2[i], y1[j], y2[j]);
}

// The minor of rows i, j and k of the 6x3 matrix with columns y1, y2, y3.
inline Complex Minor3(const Vector<6>& y1, const Vector<6>& y2,
                      const Vector<6>& y3, int i, int j, int k) {
  return y1[i] * Det2(y2[j], y3[j], y2[k], y3[k]) -
         y2[i] * Det2(y1[j], y3[j], y1[k], y3[k]) +
         y3[i] * Det2(y1[j], y2[j], y1[k], y2[k]);
}

}  // namespace EquationDetails

/** @brief The five independent second-order minors of two 4-vectors. */
inline Vector<5> Minors(const Vector<4>& y1, const Vector<4>& y2) {
  using EquationDetails::Minor2;
  return {Minor2(y1, y2, 0, 1), Minor2(y1, y2, 0, 2), Minor2(y1, y2, 0, 3),
          Minor2(y1, y2, 1, 2), Minor2(y1, y2, 2, 3)};
}

/** @brief The fourteen independent third-order minors of three 6-vectors. */
inline Vector<14> Minors(const Vector<6>& y1, const Vector<6>& y2,
                         const Vector<6>& y3) {
  auto m = [&](int i, int j, int k) {
    return EquationDetails::Minor3(y1, y2, y3, i - 1, j - 1, k - 1);
  };
  return {m(1, 2, 3), m(1, 2, 4), m(1, 2, 5), m(1, 2, 6), m(1, 3, 4),
          m(1, 3, 5), m(1, 4, 5), m(1, 4, 6), m(1, 5, 6), m(2, 3, 4),
          m(2, 4, 5), m(2, 4, 6), m(3, 4, 5), m(4, 5, 6)};
}

//----------------------------------------------------------------------------
// Starting solutions: those regular at the centre of a homogeneous sphere.
//----------------------------------------------------------------------------

/** @brief The toroidal starting solution. */
inline Vector<2> ToroidalStart(const Context& c, const ComplexMaterial& m) {
  const auto l = static_cast<double>(c.Degree());
  const auto r = m.r;
  const auto q = BesselRatio(c.Degree(), c.Omega2() * m.rho * r * r / m.mu);
  return {c.Zeta() * r, c.Zeta() * m.mu * ((l - 1) - q)};
}

/** @brief The radial starting solution, with or without gravity. */
inline Vector<2> RadialStart(const Context& c, const ComplexMaterial& m,
                             bool gravity) {
  const auto r = m.r;
  const auto modulus = m.kappa + 4.0 * m.mu / 3.0;
  const auto shift = gravity ? 16 * PiG * m.rho / 3 : 0.0;
  const auto gamma = m.rho * (c.Omega2() + shift) / modulus;
  const auto q = BesselRatio(0, gamma * r * r);
  return {q, r * modulus * gamma - 4.0 * m.mu * q / r};
}

/** @brief The minors of the spheroidal starting solutions in a solid, with
 * full self-gravitation. */
inline Vector<14> SolidStart6(const Context& c, const ComplexMaterial& m) {
  const auto l = static_cast<double>(c.Degree());
  const auto zeta = c.Zeta();
  const auto zeta2 = c.Zeta2();
  const auto omega2 = c.Omega2();
  const auto r = m.r;
  const auto rho = m.rho;
  const auto alpha2 = (m.kappa + 4.0 * m.mu / 3.0) / rho;
  const auto beta2 = m.mu / rho;
  const auto tmp1 = omega2 / beta2;
  const auto tmp2 = (omega2 + 16 * PiG * rho / 3) / alpha2;
  const auto tmp3 = std::pow(8 * PiG * zeta * rho / 3, 2) / (alpha2 * beta2);

  auto y = std::array<Vector<6>, 3>{};
  for (auto i : {0, 1}) {
    const auto sign = i == 0 ? 1.0 : -1.0;
    const auto gamma =
        0.5 * (tmp1 + tmp2) +
        0.5 * sign * std::sqrt((tmp1 - tmp2) * (tmp1 - tmp2) + tmp3);
    const auto zt = 0.75 * beta2 * (gamma - tmp1) / (PiG * rho);
    const auto xi = zt - (l + 1);
    const auto q = BesselRatio(c.Degree(), gamma * r * r);
    y[i][0] = l * xi - zt * q;
    y[i][1] = zeta * (xi + q);
    y[i][2] = -4 * PiG * rho * r * zt;
    y[i][3] = -rho * alpha2 * zt * r * gamma +
              2.0 * m.mu * l * (l - 1) * xi / r +
              2.0 * m.mu * (2.0 * zt + zeta2) * q / r;
    y[i][4] = m.mu * zeta *
              (r * gamma + 2.0 * (l - 1) * xi / r - 2.0 * (zt + 1.0) * q / r);
    y[i][5] = -rho * (zeta2 + (l + 1) * zt);
  }
  y[2][0] = l;
  y[2][1] = zeta;
  y[2][2] = (omega2 - 4 * PiG * rho * l / 3) * r;
  y[2][3] = 2.0 * l * (l - 1) * m.mu / r;
  y[2][4] = 2.0 * zeta * (l - 1) * m.mu / r;
  y[2][5] =
      0.25 * ((2 * l + 1) * omega2 - 8 * PiG * l * (l - 1) * rho / 3) / PiG;
  return Minors(y[0], y[1], y[2]);
}

/** @brief The minors of the spheroidal starting solutions in a solid, without
 * gravity (used in the Cowling approximation too). */
inline Vector<5> SolidStart4(const Context& c, const ComplexMaterial& m) {
  const auto l = static_cast<double>(c.Degree());
  const auto zeta = c.Zeta();
  const auto omega2 = c.Omega2();
  const auto r = m.r;
  const auto rho = m.rho;
  const auto alpha2 = (m.kappa + 4.0 * m.mu / 3.0) / rho;
  const auto beta2 = m.mu / rho;
  auto y1 = Vector<4>{};
  auto y2 = Vector<4>{};
  {
    // With z = (kr)^2 for the P wavenumber.
    const auto z = omega2 * r * r / alpha2;
    const auto q = BesselRatio(c.Degree(), z);
    y1[0] = l - q;
    y1[1] = zeta;
    y1[2] = (-rho * alpha2 * z + 2.0 * m.mu * (l * (l - 1) + 2.0 * q)) / r;
    y1[3] = 2.0 * zeta * m.mu * ((l - 1) - q) / r;
  }
  {
    // With z = (kr)^2 for the S wavenumber.
    const auto z = omega2 * r * r / beta2;
    const auto q = BesselRatio(c.Degree(), z);
    y2[0] = -l * (l + 1);
    y2[1] = -zeta * (l + 1) + q;
    y2[2] = -2.0 * m.mu * (l * (l * l - 1) - l * (l + 1) * q) / r;
    y2[3] = zeta * m.mu * (z - 2.0 * (l * l - 1) - 2.0 * q) / r;
  }
  return Minors(y1, y2);
}

/** @brief The minors of the spheroidal starting solutions in a fluid, with
 * full self-gravitation. */
inline Vector<5> FluidStart4(const Context& c, const ComplexMaterial& m) {
  const auto l = static_cast<double>(c.Degree());
  const auto zeta2 = c.Zeta2();
  const auto omega2 = c.Omega2();
  const auto r = m.r;
  const auto rho = m.rho;
  const auto alpha2 = m.kappa / rho;
  const auto gamg = 4 * PiG * rho / 3;
  const auto gamma =
      (omega2 + 4 * gamg - zeta2 * gamg * gamg * c.OmegaInverse2()) / alpha2;
  const auto zt = -omega2 / gamg;
  const auto xi = zt - (l + 1);
  const auto q = BesselRatio(c.Degree(), gamma * r * r);
  const auto y1 =
      Vector<4>{l * xi - zt * q, -4 * PiG * rho * r * zt,
                -rho * alpha2 * zt * r * gamma, -rho * (zeta2 + (l + 1) * zt)};
  const auto y2 = Vector<4>{
      l, (omega2 - 4 * PiG * rho * l / 3) * r, 0.0,
      0.25 * ((2 * l + 1) * omega2 - 8 * PiG * l * (l - 1) * rho / 3) / PiG};
  return Minors(y1, y2);
}

/** @brief The spheroidal starting solution in a fluid, without gravity (used
 * in the Cowling approximation too). */
inline Vector<2> FluidStart2(const Context& c, const ComplexMaterial& m) {
  const auto l = static_cast<double>(c.Degree());
  const auto r = m.r;
  const auto alpha2 = m.kappa / m.rho;
  const auto z = c.Omega2() * r * r / alpha2;
  const auto q = BesselRatio(c.Degree(), z);
  return {l - q, -m.rho * alpha2 * z / r};
}

//----------------------------------------------------------------------------
// Solid-fluid boundaries, the secular determinant, the b-vector and the
// solution.
//----------------------------------------------------------------------------

/** @brief Third-order minors in a solid to second-order minors in a fluid. */
inline Vector<5> SolidToFluid(const Vector<14>& m1) {
  return {-m1[5], -m1[6], m1[8], -m1[12], m1[13]};
}

/** @brief Second-order minors in a solid to the fluid solution (no gravity or
 * Cowling). */
inline Vector<2> SolidToFluid(const Vector<5>& m1) { return {-m1[2], -m1[4]}; }

/** @brief Second-order minors in a fluid to third-order minors in a solid. */
inline Vector<14> FluidToSolid(const Vector<5>& m1) {
  auto m2 = Vector<14>{};
  m2[0] = m1[0];
  m2[1] = m1[1];
  m2[3] = m1[2];
  m2[9] = -m1[3];
  m2[11] = -m1[4];
  return m2;
}

/** @brief The fluid solution to second-order minors in a solid (no gravity or
 * Cowling). */
inline Vector<5> FluidToSolid(const Vector<2>& m1) {
  auto m2 = Vector<5>{};
  m2[0] = m1[0];
  m2[3] = -m1[1];
  return m2;
}

/** @brief The secular determinant from minors integrated from below and from
 * above (4x4 systems). */
inline Complex SecularDeterminant(const Vector<5>& m1, const Vector<5>& m2) {
  return m1[4] * m2[0] + 2.0 * m1[1] * m2[1] + m1[3] * m2[2] + m1[2] * m2[3] +
         m1[0] * m2[4];
}

/** @brief The secular determinant from minors integrated from below and from
 * above (6x6 system). */
inline Complex SecularDeterminant(const Vector<14>& m1, const Vector<14>& m2) {
  return -m1[13] * m2[0] + 2.0 * m1[6] * m2[1] + 2.0 * m1[10] * m2[2] +
         m1[12] * m2[3] + 2.0 * m1[7] * m2[4] + m1[11] * m2[5] -
         2.0 * m1[1] * m2[6] - 2.0 * m1[4] * m2[7] - m1[9] * m2[8] +
         m1[8] * m2[9] - 2.0 * m1[2] * m2[10] - m1[5] * m2[11] -
         m1[3] * m2[12] + m1[0] * m2[13];
}

/** @brief The starting b-vector at the source from the minors and a jump (4x4
 * systems). */
inline Vector<4> BStart(const Vector<5>& m, const Vector<4>& s) {
  return {-m[4] * s[1] - m[1] * s[2] - m[3] * s[3],
          m[4] * s[0] - m[2] * s[2] + m[1] * s[3],
          m[1] * s[0] + m[2] * s[1] - m[0] * s[3],
          m[3] * s[0] - m[1] * s[1] + m[0] * s[2]};
}

/** @brief The starting b-vector at the source from the minors and a jump (6x6
 * system). */
inline Vector<14> BStart(const Vector<14>& m, const Vector<6>& s) {
  auto b = Vector<14>{};
  b[0] = m[13] * s[2] - m[6] * s[3] - m[10] * s[4] - m[12] * s[5];
  b[1] = -m[13] * s[1] - m[7] * s[3] - m[11] * s[4] + m[10] * s[5];
  b[2] = m[6] * s[1] + m[7] * s[2] + m[1] * s[4] + m[4] * s[5];
  b[3] = m[10] * s[1] + m[11] * s[2] - m[1] * s[3] + m[9] * s[5];
  b[4] = m[12] * s[1] - m[10] * s[2] - m[4] * s[3] - m[9] * s[4];
  b[5] = m[13] * s[0] - m[8] * s[3] + m[7] * s[4] - m[6] * s[5];
  b[6] = -m[6] * s[0] + m[8] * s[2] + m[2] * s[4] + m[5] * s[5];
  b[7] = -m[10] * s[0] - m[7] * s[2] - m[2] * s[3] - m[4] * s[5];
  b[8] = -m[12] * s[0] + m[6] * s[2] - m[5] * s[3] + m[4] * s[4];
  b[9] = -m[7] * s[0] - m[8] * s[1] + m[3] * s[4] - m[2] * s[5];
  b[10] = -m[11] * s[0] + m[7] * s[1] - m[3] * s[3] + m[1] * s[5];
  b[11] = m[1] * s[0] + m[2] * s[1] + m[3] * s[2] - m[0] * s[5];
  b[12] = m[4] * s[0] + m[5] * s[1] - m[2] * s[2] + m[0] * s[4];
  b[13] = m[9] * s[0] - m[4] * s[1] + m[1] * s[2] - m[0] * s[3];
  return b;
}

/** @brief The solution from the minors and the b-vector (4x4 systems). */
inline Vector<4> Solution(const Vector<5>& m, const Vector<4>& b) {
  return {m[0] * b[1] + m[1] * b[2] + m[2] * b[3],
          -m[0] * b[0] + m[3] * b[2] - m[1] * b[3],
          -m[1] * b[0] - m[3] * b[1] + m[4] * b[3],
          -m[2] * b[0] + m[1] * b[1] - m[4] * b[2]};
}

/** @brief The solution from the minors and the b-vector (6x6 system). */
inline Vector<6> Solution(const Vector<14>& m, const Vector<14>& b) {
  auto y = Vector<6>{};
  y[0] = m[2] * b[2] + m[0] * b[5] + m[1] * b[6] + 2.0 * m[2] * b[7] +
         m[3] * b[8] + m[4] * b[9] + m[5] * b[10] + m[6] * b[11] +
         m[7] * b[12] + m[8] * b[13];
  y[1] = -m[0] * b[1] - m[1] * (2.0 * b[2] + b[7]) - m[2] * b[3] - m[3] * b[4] -
         m[4] * b[10] - m[7] * b[13] + m[9] * b[9] + m[10] * b[11] +
         m[11] * b[12];
  y[2] = m[0] * b[0] - m[1] * b[8] + m[2] * b[4] + m[4] * (b[7] - b[2]) -
         m[5] * b[3] + m[6] * b[13] - m[9] * b[6] - m[10] * b[12] +
         m[12] * b[11];
  y[3] = m[1] * b[0] + m[4] * b[1] - m[6] * b[3] - m[7] * b[4] + m[9] * b[5] -
         m[10] * (b[2] + 2.0 * b[7]) - m[11] * b[8] - m[12] * b[10] +
         m[13] * b[13];
  y[4] = m[2] * b[0] - m[4] * b[5] + m[5] * b[1] + m[6] * (2.0 * b[2] + b[7]) +
         m[7] * b[8] - m[8] * b[4] + m[10] * b[6] + m[12] * b[9] -
         m[13] * b[12];
  y[5] = m[1] * b[5] - m[2] * b[1] + m[3] * b[0] + m[6] * b[10] +
         m[7] * (b[2] - b[7]) + m[8] * b[3] - m[10] * b[9] + m[11] * b[6] +
         m[13] * b[11];
  return y;
}

}  // namespace YSpec
