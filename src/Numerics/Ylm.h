//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_YLM_H
#define QMCPLUSPLUS_YLM_H

#include <cmath>
#include <complex>
#include "OhmmsPETE/TinyVector.h"
#include "config/stdlib/Constants.h"

namespace qmcplusplus
{
template<typename T>
inline T LegendrePll(int l, T x)
{
  if (l == 0)
    return 1.0;
  else
  {
    T sqt     = std::sqrt(1.0 - x) * std::sqrt(1.0 + x);
    T val     = 1.0;
    T dblfact = 1.0;
    for (int i = 1; i <= l; i++)
    {
      val *= -sqt;
      val *= dblfact;
      dblfact += 2.0;
    }
    return val;
  }
}

template<typename T>
inline T LegendrePlm(int l, int m, T x)
{
  if (m < 0)
  {
    m        = std::abs(m);
    T posval = LegendrePlm(l, m, x);
    T sign   = (m % 2 == 0) ? 1.0 : -1.0;
    T mfact  = 1.0;
    for (int i = 2; i <= (l - m); i++)
      mfact *= static_cast<T>(i);
    T pfact = 1.0;
    for (int i = 2; i <= (l + m); i++)
      pfact *= static_cast<T>(i);
    return posval * sign * mfact / pfact;
  }
  // Now we can assume that m is not negative
  T pmm   = LegendrePll(m, x);
  T pmp1m = x * (2 * m + 1) * pmm;
  if (m == l)
    return pmm;
  else if (l == (m + 1))
    return pmp1m;
  else
  // Use recursive formula
  {
    T Plm2m = pmm;
    T Plm1m = pmp1m;
    T Pl    = 0;
    for (int i = m + 2; i <= l; i++)
    {
      Pl    = (1.0 / static_cast<T>(i - m)) * (x * (2 * i - 1) * Plm1m - (i + m - 1) * Plm2m);
      Plm2m = Plm1m;
      Plm1m = Pl;
    }
    return Pl;
  }
}

/** calculates Ylm
 * param[in] l angular momentum
 * param[in] m magnetic quantum number
 * param[in] r position vector. Note: This must be a unit vector and in the order [z,x,y] to be correct
 */
template<typename T>
inline std::complex<T> Ylm(int l, int m, const TinyVector<T, 3>& r)
{
  T costheta = r[0];
  T phi      = std::atan2(r[2], r[1]);
  int lmm    = l - m;
  int lpm    = l + m;
  T mfact    = 1.0;
  T pfact    = 1.0;
  for (int i = lmm; i > 0; i--)
    mfact *= static_cast<T>(i);
  for (int i = lpm; i > 0; i--)
    pfact *= static_cast<T>(i);
  T prefactor = std::sqrt(static_cast<T>(2 * l + 1) * mfact / (4.0 * M_PI * pfact));
  prefactor *= LegendrePlm(l, m, costheta);
  return std::complex<T>(prefactor * std::cos(m * phi), prefactor * std::sin(m * phi));
}

/** calculate the derivative of a Ylm with respect to theta and phi
 * param[in] l: angular momentum
 * param[in] m: magnetic quantum number
 * param[in] r: cartesian position align with [z,x,y]. note: MUST BE UNIT VECTOR to be consistent with Ylm above
 * param[out] theta_deriv: derivative of Ylm with respect to theta
 * param[out] phi_deriv: derivative of Ylm with respect to phi
 * param[in] conj: true if we are taking derivatives of conj(Ylm)
 */
template<typename T>
inline void derivYlmSpherical(const int l,
                              const int m,
                              const TinyVector<T, 3>& r,
                              std::complex<T>& theta_deriv,
                              std::complex<T>& phi_deriv,
                              const bool conj)
{
  T theta                = std::acos(r[0]);
  T phi                  = std::atan2(r[2], r[1]);
  std::complex<T> emiphi = std::complex<T>(std::cos(phi), -std::sin(phi));
  theta_deriv            = static_cast<T>(m) / std::tan(theta) * Ylm(l, m, r);
  theta_deriv += std::sqrt(static_cast<T>((l - m) * (l + m + 1))) * emiphi * Ylm(l, m + 1, r);
  phi_deriv = std::complex<T>(0.0, static_cast<T>(m)) * Ylm(l, m, r);
  if (conj)
  {
    theta_deriv = std::conj(theta_deriv);
    phi_deriv   = std::conj(phi_deriv);
  }
}

/** wrapper for Ylm, which can take any normal position vector as an argument
 * param[in] l angular momentum
 * param[in] m magnetic quantum number
 * param[in] r is a position vector. This does not have to be normalized and is in the standard ordering [x,y,z]
 */
template<typename T>
inline std::complex<T> sphericalHarmonic(const int l, const int m, const TinyVector<T, 3>& r)
{
  TinyVector<T, 3> unit;
  T norm  = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
  unit[0] = r[2] / norm;
  unit[1] = r[0] / norm;
  unit[2] = r[1] / norm;
  return Ylm(l, m, unit);
}

/** get cartesian derivatives of spherical Harmonics. This takes a arbitrary position vector (x,y,z) and returns (d/dx, d/dy, d/dz)Ylm
 * param[in] l angular momentum 
 * param[in] m magnetic quantum number
 * param[in] r position vector. This does not have to be normalized and is in the standard ordering [x,y,z]
 * param[out] grad (d/dx, d/dy, d/dz) of Ylm
 */
template<typename T>
inline void sphericalHarmonicGrad(const int l,
                                  const int m,
                                  const TinyVector<T, 3>& r,
                                  TinyVector<std::complex<T>, 3>& grad)
{
  TinyVector<T, 3> unit;
  T norm  = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
  unit[0] = r[2] / norm;
  unit[1] = r[0] / norm;
  unit[2] = r[1] / norm;
  std::complex<T> dth, dph;
  derivYlmSpherical(l, m, unit, dth, dph, false);

  T dth_dx = r[0] * r[2] / (norm * norm * std::sqrt(r[0] * r[0] + r[1] * r[1]));
  T dph_dx = -r[1] / (r[0] * r[0] + r[1] * r[1]);
  T dth_dy = r[1] * r[2] / (norm * norm * std::sqrt(r[0] * r[0] + r[1] * r[1]));
  T dph_dy = r[0] / (r[0] * r[0] + r[1] * r[1]);
  T dth_dz = -std::sqrt(r[0] * r[0] + r[1] * r[1]) / (norm * norm);
  T dph_dz = 0;

  grad[0] = dth * dth_dx + dph * dph_dx;
  grad[1] = dth * dth_dy + dph * dph_dy;
  grad[2] = dth * dth_dz + dph * dph_dz;
}

} // namespace qmcplusplus
#endif
