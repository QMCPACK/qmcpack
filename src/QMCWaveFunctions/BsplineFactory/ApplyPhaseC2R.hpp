//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "OMPTarget/OMPTargetMath.hpp"

namespace qmcplusplus
{
namespace C2R
{
template<typename ST, typename TT>
inline void assign_v(ST x,
                     ST y,
                     ST z,
                     TT* restrict results_scratch_ptr,
                     const ST* restrict offload_scratch_ptr,
                     const ST* restrict myKcart_ptr,
                     size_t myKcart_padded_size,
                     size_t first_spo,
                     int nComplexBands,
                     int index)
{
  const ST* restrict kx = myKcart_ptr;
  const ST* restrict ky = myKcart_ptr + myKcart_padded_size;
  const ST* restrict kz = myKcart_ptr + myKcart_padded_size * 2;

  const ST* restrict val = offload_scratch_ptr;
  TT* restrict psi_s     = results_scratch_ptr;

  const size_t jr = index << 1;
  const size_t ji = jr + 1;
  //phase
  ST s, c, p = -(x * kx[index] + y * ky[index] + z * kz[index]);
  omptarget::sincos(p, &s, &c);

  const ST val_r        = val[jr];
  const ST val_i        = val[ji];
  const size_t psiIndex = first_spo + index + omptarget::min(index, nComplexBands);
  psi_s[psiIndex]       = val_r * c - val_i * s;
  if (index < nComplexBands)
    psi_s[psiIndex + 1] = val_i * c + val_r * s;
}

/** assign_vgl
   */
template<typename ST, typename TT>
inline void assign_vgl(ST x,
                       ST y,
                       ST z,
                       TT* restrict results_scratch_ptr,
                       size_t orb_padded_size,
                       const ST* mKK_ptr,
                       const ST* restrict offload_scratch_ptr,
                       size_t spline_padded_size,
                       const ST G[9],
                       const ST* myKcart_ptr,
                       size_t myKcart_padded_size,
                       size_t first_spo,
                       int nComplexBands,
                       int index)
{
  constexpr ST two(2);
  const ST &g00 = G[0], &g01 = G[1], &g02 = G[2], &g10 = G[3], &g11 = G[4], &g12 = G[5], &g20 = G[6], &g21 = G[7],
           &g22 = G[8];

  const ST* restrict k0 = myKcart_ptr;
  const ST* restrict k1 = myKcart_ptr + myKcart_padded_size;
  const ST* restrict k2 = myKcart_ptr + myKcart_padded_size * 2;

  const ST* restrict val   = offload_scratch_ptr + spline_padded_size * SoAFields3D::VAL;
  const ST* restrict g0    = offload_scratch_ptr + spline_padded_size * SoAFields3D::GRAD0;
  const ST* restrict g1    = offload_scratch_ptr + spline_padded_size * SoAFields3D::GRAD1;
  const ST* restrict g2    = offload_scratch_ptr + spline_padded_size * SoAFields3D::GRAD2;
  const ST* restrict lcart = offload_scratch_ptr + spline_padded_size * SoAFields3D::LAPL;

  const size_t jr = index << 1;
  const size_t ji = jr + 1;

  const ST kX    = k0[index];
  const ST kY    = k1[index];
  const ST kZ    = k2[index];
  const ST val_r = val[jr];
  const ST val_i = val[ji];

  //phase
  ST s, c, p = -(x * kX + y * kY + z * kZ);
  omptarget::sincos(p, &s, &c);

  //dot(PrimLattice.G,myG[j])
  const ST dX_r = g00 * g0[jr] + g01 * g1[jr] + g02 * g2[jr];
  const ST dY_r = g10 * g0[jr] + g11 * g1[jr] + g12 * g2[jr];
  const ST dZ_r = g20 * g0[jr] + g21 * g1[jr] + g22 * g2[jr];

  const ST dX_i = g00 * g0[ji] + g01 * g1[ji] + g02 * g2[ji];
  const ST dY_i = g10 * g0[ji] + g11 * g1[ji] + g12 * g2[ji];
  const ST dZ_i = g20 * g0[ji] + g21 * g1[ji] + g22 * g2[ji];

  // \f$\nabla \psi_r + {\bf k}\psi_i\f$
  const ST gX_r = dX_r + val_i * kX;
  const ST gY_r = dY_r + val_i * kY;
  const ST gZ_r = dZ_r + val_i * kZ;
  const ST gX_i = dX_i - val_r * kX;
  const ST gY_i = dY_i - val_r * kY;
  const ST gZ_i = dZ_i - val_r * kZ;

  const ST lap_r = lcart[jr] + mKK_ptr[index] * val_r + two * (kX * dX_i + kY * dY_i + kZ * dZ_i);
  const ST lap_i = lcart[ji] + mKK_ptr[index] * val_i - two * (kX * dX_r + kY * dY_r + kZ * dZ_r);

  TT* restrict psi    = results_scratch_ptr;
  TT* restrict dpsi_x = results_scratch_ptr + orb_padded_size;
  TT* restrict dpsi_y = results_scratch_ptr + orb_padded_size * 2;
  TT* restrict dpsi_z = results_scratch_ptr + orb_padded_size * 3;
  TT* restrict d2psi  = results_scratch_ptr + orb_padded_size * 4;

  const size_t psiIndex = first_spo + index + omptarget::min(index, nComplexBands);

  psi[psiIndex]    = c * val_r - s * val_i;
  d2psi[psiIndex]  = c * lap_r - s * lap_i;
  dpsi_x[psiIndex] = c * gX_r - s * gX_i;
  dpsi_y[psiIndex] = c * gY_r - s * gY_i;
  dpsi_z[psiIndex] = c * gZ_r - s * gZ_i;

  if (index < nComplexBands)
  {
    psi[psiIndex + 1]    = c * val_i + s * val_r;
    d2psi[psiIndex + 1]  = c * lap_i + s * lap_r;
    dpsi_x[psiIndex + 1] = c * gX_i + s * gX_r;
    dpsi_y[psiIndex + 1] = c * gY_i + s * gY_r;
    dpsi_z[psiIndex + 1] = c * gZ_i + s * gZ_r;
  }
}
} // namespace C2R
} // namespace qmcplusplus
