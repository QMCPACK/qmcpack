//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////

/*
 DO NOT MAKE PERMANENT EDITS IN THIS FILE
 This file is generated from src/Numerics/codegen/gen_spherical_tensor.py and
 SoaSphericalTensor.cpp.in.

 Edit the template or generator and rerun gen_spherical_tensor.py.
*/

#include "SoaSphericalTensor.h"

namespace qmcplusplus
{
PRAGMA_OFFLOAD("omp begin declare target")
template<typename T>
void SoaSphericalTensor<T>::gradient_recurrence(const int l,
                                                       const int m,
                                                       const T fac,
                                                       const T* restrict lower,
                                                       T& gx,
                                                       T& gy,
                                                       T& gz)
{
  constexpr T czero(0);
  constexpr T ahalf(0.5);
  const int ma = std::abs(m);
  const int lm = index(l - 1, 0);
  const T cp = std::sqrt(fac * (l - ma - 1) * (l - ma));
  const T cm = std::sqrt(fac * (l + ma - 1) * (l + ma));
  const T c0 = std::sqrt(fac * (l - ma) * (l + ma));

  T dpr, dpi, dmr, dmi;
  gz = (l > ma) ? c0 * lower[lm + m] : czero;
  if (l > ma + 1)
  {
    dpr = cp * lower[lm + ma + 1];
    dpi = cp * lower[lm - ma - 1];
  }
  else
  {
    dpr = czero;
    dpi = czero;
  }
  if (l > 1)
  {
    switch (ma)
    {
    case 0:
      dmr = -cm * lower[lm + 1];
      dmi = cm * lower[lm - 1];
      break;
    case 1:
      dmr = cm * lower[lm];
      dmi = czero;
      break;
    default:
      dmr = cm * lower[lm + ma - 1];
      dmi = cm * lower[lm - ma + 1];
    }
  }
  else
  {
    dmr = cm * lower[lm];
    dmi = czero;
  }
  if (m < 0)
  {
    gx = ahalf * (dpi - dmi);
    gy = -ahalf * (dpr + dmr);
  }
  else
  {
    gx = ahalf * (dpr - dmr);
    gy = ahalf * (dpi + dmi);
  }
}
PRAGMA_OFFLOAD("omp end declare target")
template<typename T>
void SoaSphericalTensor<T>::evaluateVGH(T x, T y, T z)
{
  // The Hessian is the same recurrence applied a second time to the gradients, and
  // the recurrence consumes unnormalized quantities. Stopping short of normalizing
  // hands them over in exactly that form, so nothing has to be undone here;
  // norm_factor_ goes on at the end, once the Hessian loop is done reading them.
  const int Nlm = cYlm.size();
  evaluateVGL_bare(x, y, z, cYlm.data(), Lmax, factorL_.data(), factorLM_.data(), factor2L_.data(), cYlm.capacity());

  constexpr T czero(0);
  constexpr T ahalf(0.5);
  // not restrict-qualified: normalize_vg writes these same rows below
  const T* gYlmX = cYlm.data(1);
  const T* gYlmY = cYlm.data(2);
  const T* gYlmZ = cYlm.data(3);
  T* restrict hYlmXX      = cYlm.data(4);
  T* restrict hYlmXY      = cYlm.data(5);
  T* restrict hYlmXZ      = cYlm.data(6);
  T* restrict hYlmYY      = cYlm.data(7);
  T* restrict hYlmYZ      = cYlm.data(8);
  T* restrict hYlmZZ      = cYlm.data(9);

  hYlmXX[0] = czero;
  hYlmXY[0] = czero;
  hYlmXZ[0] = czero;
  hYlmYY[0] = czero;
  hYlmYZ[0] = czero;
  hYlmZZ[0] = czero;

  for (int l = 1; l <= Lmax; ++l)
  {
    const T fac = factor2L_[l];
    for (int m = -l; m <= l; ++m)
    {
      const int lm = index(l, m);

      // The recurrence coefficients do not depend on position, so differentiating
      // a second time is the same recurrence applied to the bare l-1 gradients that
      // evaluateVGL_bare just stored, with the target harmonic's norm_factor_
      // applied on the way out.
      const auto differentiate_gradient = [&](const T* lower_derivative) -> TinyVector<T, 3> {
        T dx, dy, dz;
        gradient_recurrence(l, m, fac, lower_derivative, dx, dy, dz);
        return {dx, dy, dz};
      };

      // d_gx[i] is the i-th derivative of the stored x gradient, so the mixed
      // second derivatives each come out of two separate recurrence passes.
      const TinyVector<T, 3> d_gx = differentiate_gradient(gYlmX);
      const TinyVector<T, 3> d_gy = differentiate_gradient(gYlmY);
      const TinyVector<T, 3> d_gz = differentiate_gradient(gYlmZ);

      // The Hessian is symmetric analytically, so d_gx[1] and d_gy[0] differ
      // only by roundoff. Both are already in hand, so averaging them is free
      // and keeps the stored Hessian exactly symmetric.
      const T norm = norm_factor_[lm];
      hYlmXX[lm]   = norm * d_gx[0];
      hYlmXY[lm]   = norm * ahalf * (d_gx[1] + d_gy[0]);
      hYlmXZ[lm]   = norm * ahalf * (d_gx[2] + d_gz[0]);
      hYlmYY[lm]   = norm * d_gy[1];
      hYlmYZ[lm]   = norm * ahalf * (d_gy[2] + d_gz[1]);
      hYlmZZ[lm]   = norm * d_gz[2];
    }
  }

  normalize_vg(cYlm.data(), Nlm, norm_factor_.data(), cYlm.capacity());
}

template<typename T>
void SoaSphericalTensor<T>::evaluateVGHGH(T x, T y, T z)
{
  throw std::runtime_error("SoaSphericalTensor<T>::evaluateVGHGH(x,y,z):  Not implemented\n");
}

template class SoaSphericalTensor<float>;
template class SoaSphericalTensor<double>;

} // namespace qmcplusplus
