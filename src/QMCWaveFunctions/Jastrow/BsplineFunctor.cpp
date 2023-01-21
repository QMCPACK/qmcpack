//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: John R. Gergely,  University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#include "BsplineFunctor.h"

namespace qmcplusplus
{
template<typename REAL>
void BsplineFunctor<REAL>::mw_evaluateVGL(const int iat,
                                          const int num_groups,
                                          const BsplineFunctor* const functors[],
                                          const int n_src,
                                          const int* grp_ids,
                                          const int nw,
                                          REAL* mw_vgl,        // [nw][DIM+2]
                                          const int n_padded,
                                          const REAL* mw_dist, // [nw][DIM+1][n_padded]
                                          REAL* mw_cur_allu,   // [nw][3][n_padded]
                                          Vector<char, OffloadPinnedAllocator<char>>& transfer_buffer)
{
  constexpr unsigned DIM = OHMMS_DIM;
  static_assert(DIM == 3, "only support 3D due to explicit x,y,z coded.");
  const size_t dist_stride = n_padded * (DIM + 1);

  /* transfer buffer used for
     * Bspline coefs device pointer sizeof(T*), DeltaRInv sizeof(T) and cutoff_radius sizeof(T)
     * these contents change based on the group of the target particle, so it is prepared per call.
     */
  transfer_buffer.resize((sizeof(REAL*) + sizeof(REAL) * 2) * num_groups);
  REAL** mw_coefs_ptr        = reinterpret_cast<REAL**>(transfer_buffer.data());
  REAL* mw_DeltaRInv_ptr     = reinterpret_cast<REAL*>(transfer_buffer.data() + sizeof(REAL*) * num_groups);
  REAL* mw_cutoff_radius_ptr = mw_DeltaRInv_ptr + num_groups;
  for (int ig = 0; ig < num_groups; ig++)
    if (functors[ig])
    {
      mw_coefs_ptr[ig]         = functors[ig]->spline_coefs_->device_data();
      mw_DeltaRInv_ptr[ig]     = functors[ig]->DeltaRInv;
      mw_cutoff_radius_ptr[ig] = functors[ig]->cutoff_radius;
    }
    else
    {
      mw_coefs_ptr[ig]         = nullptr;
      mw_DeltaRInv_ptr[ig]     = 0.0;
      mw_cutoff_radius_ptr[ig] = 0.0; // important! Prevent spline evaluation to access nullptr.
    }

  auto* transfer_buffer_ptr = transfer_buffer.data();

  PRAGMA_OFFLOAD("omp target teams distribute map(always, to: transfer_buffer_ptr[:transfer_buffer.size()]) \
                    map(to: grp_ids[:n_src]) \
                    map(to: mw_dist[:dist_stride*nw]) \
                    map(from: mw_cur_allu[:n_padded*3*nw]) \
                    map(always, from: mw_vgl[:(DIM+2)*nw])")
  for (int ip = 0; ip < nw; ip++)
  {
    REAL val_sum(0);
    REAL grad_x(0);
    REAL grad_y(0);
    REAL grad_z(0);
    REAL lapl(0);

    const REAL* dist   = mw_dist + ip * dist_stride;
    const REAL* dipl_x = dist + n_padded;
    const REAL* dipl_y = dist + n_padded * 2;
    const REAL* dipl_z = dist + n_padded * 3;

    REAL** mw_coefs        = reinterpret_cast<REAL**>(transfer_buffer_ptr);
    REAL* mw_DeltaRInv     = reinterpret_cast<REAL*>(transfer_buffer_ptr + sizeof(REAL*) * num_groups);
    REAL* mw_cutoff_radius = mw_DeltaRInv + num_groups;

    REAL* cur_allu = mw_cur_allu + ip * n_padded * 3;

    PRAGMA_OFFLOAD("omp parallel for reduction(+: val_sum, grad_x, grad_y, grad_z, lapl)")
    for (int j = 0; j < n_src; j++)
    {
      if (j == iat)
        continue;
      const int ig       = grp_ids[j];
      const REAL* coefs  = mw_coefs[ig];
      REAL DeltaRInv     = mw_DeltaRInv[ig];
      REAL cutoff_radius = mw_cutoff_radius[ig];

      REAL r = dist[j];
      REAL u(0);
      REAL dudr(0);
      REAL d2udr2(0);
      if (r < cutoff_radius)
      {
        u = evaluate_impl(dist[j], coefs, DeltaRInv, dudr, d2udr2);
        dudr *= REAL(1) / r;
      }
      // save u, dudr/r and d2udr2 to cur_allu
      cur_allu[j]                = u;
      cur_allu[j + n_padded]     = dudr;
      cur_allu[j + n_padded * 2] = d2udr2;
      val_sum += u;
      lapl += d2udr2 + (DIM - 1) * dudr;
      grad_x += dudr * dipl_x[j];
      grad_y += dudr * dipl_y[j];
      grad_z += dudr * dipl_z[j];
    }

    REAL* vgl = mw_vgl + ip * (DIM + 2);
    vgl[0]    = val_sum;
    vgl[1]    = grad_x;
    vgl[2]    = grad_y;
    vgl[3]    = grad_z;
    vgl[4]    = -lapl;
  }
}

template<typename REAL>
void BsplineFunctor<REAL>::mw_evaluateV(const int num_groups,
                                        const BsplineFunctor<REAL>* const functors[],
                                        const int n_src,
                                        const int* grp_ids,
                                        const int num_pairs,
                                        const int* ref_at,
                                        const REAL* mw_dist,
                                        const int dist_stride,
                                        REAL* mw_vals,
                                        Vector<char, OffloadPinnedAllocator<char>>& transfer_buffer)
{
  /* transfer buffer used for
     * Bspline coefs device pointer sizeof(REAL*), DeltaRInv sizeof(REAL), cutoff_radius sizeof(REAL)
     * these contents change based on the group of the target particle, so it is prepared per call.
     */
  transfer_buffer.resize((sizeof(REAL*) + sizeof(REAL) * 2) * num_groups);
  REAL** mw_coefs_ptr        = reinterpret_cast<REAL**>(transfer_buffer.data());
  REAL* mw_DeltaRInv_ptr     = reinterpret_cast<REAL*>(transfer_buffer.data() + sizeof(REAL*) * num_groups);
  REAL* mw_cutoff_radius_ptr = mw_DeltaRInv_ptr + num_groups;
  for (int ig = 0; ig < num_groups; ig++)
    if (functors[ig])
    {
      mw_coefs_ptr[ig]         = functors[ig]->spline_coefs_->device_data();
      mw_DeltaRInv_ptr[ig]     = functors[ig]->DeltaRInv;
      mw_cutoff_radius_ptr[ig] = functors[ig]->cutoff_radius;
    }
    else
    {
      mw_coefs_ptr[ig]         = nullptr;
      mw_DeltaRInv_ptr[ig]     = 0.0;
      mw_cutoff_radius_ptr[ig] = 0.0; // important! Prevent spline evaluation to access nullptr.
    }

  auto* transfer_buffer_ptr = transfer_buffer.data();

  PRAGMA_OFFLOAD("omp target teams distribute map(always, to:transfer_buffer_ptr[:transfer_buffer.size()]) \
                    map(to: grp_ids[:n_src]) \
                    map(to:ref_at[:num_pairs], mw_dist[:dist_stride*num_pairs]) \
                    map(always, from:mw_vals[:num_pairs])")
  for (int ip = 0; ip < num_pairs; ip++)
  {
    REAL sum               = 0;
    const REAL* dist       = mw_dist + ip * dist_stride;
    REAL** mw_coefs        = reinterpret_cast<REAL**>(transfer_buffer_ptr);
    REAL* mw_DeltaRInv     = reinterpret_cast<REAL*>(transfer_buffer_ptr + sizeof(REAL*) * num_groups);
    REAL* mw_cutoff_radius = mw_DeltaRInv + num_groups;
    PRAGMA_OFFLOAD("omp parallel for reduction(+: sum)")
    for (int j = 0; j < n_src; j++)
    {
      const int ig       = grp_ids[j];
      const REAL* coefs  = mw_coefs[ig];
      REAL DeltaRInv     = mw_DeltaRInv[ig];
      REAL cutoff_radius = mw_cutoff_radius[ig];

      REAL r = dist[j];
      if (j != ref_at[ip] && r < cutoff_radius)
      {
        r *= DeltaRInv;
        REAL ipart;
        const REAL t = std::modf(r, &ipart);
        const int i  = (int)ipart;
        sum += coefs[i + 0] * (((A0 * t + A1) * t + A2) * t + A3) + coefs[i + 1] * (((A4 * t + A5) * t + A6) * t + A7) +
            coefs[i + 2] * (((A8 * t + A9) * t + A10) * t + A11) +
            coefs[i + 3] * (((A12 * t + A13) * t + A14) * t + A15);
      }
    }
    mw_vals[ip] = sum;
  }
}

template<typename REAL>
void BsplineFunctor<REAL>::mw_updateVGL(const int iat,
                                        const std::vector<bool>& isAccepted,
                                        const int num_groups,
                                        const BsplineFunctor* const functors[],
                                        const int n_src,
                                        const int* grp_ids,
                                        const int nw,
                                        REAL* mw_vgl,        // [nw][DIM+2]
                                        const int n_padded,
                                        const REAL* mw_dist, // [nw][DIM+1][n_padded]
                                        REAL* mw_allUat,     // [nw][DIM+2][n_padded]
                                        REAL* mw_cur_allu,   // [nw][3][n_padded]
                                        Vector<char, OffloadPinnedAllocator<char>>& transfer_buffer)
{
  constexpr unsigned DIM = OHMMS_DIM;
  static_assert(DIM == 3, "only support 3D due to explicit x,y,z coded.");
  const size_t dist_stride = n_padded * (DIM + 1);

  /* transfer buffer used for
     * Bspline coefs device pointer sizeof(REAL*), DeltaRInv sizeof(REAL), cutoff_radius sizeof(REAL)
     * and packed accept list at most nw * sizeof(int)
     * these contents change based on the group of the target particle, so it is prepared per call.
     */
  transfer_buffer.resize((sizeof(REAL*) + sizeof(REAL) * 2) * num_groups + nw * sizeof(int));
  REAL** mw_coefs_ptr        = reinterpret_cast<REAL**>(transfer_buffer.data());
  REAL* mw_DeltaRInv_ptr     = reinterpret_cast<REAL*>(transfer_buffer.data() + sizeof(REAL*) * num_groups);
  REAL* mw_cutoff_radius_ptr = mw_DeltaRInv_ptr + num_groups;
  int* accepted_indices =
      reinterpret_cast<int*>(transfer_buffer.data() + (sizeof(REAL*) + sizeof(REAL) * 2) * num_groups);

  for (int ig = 0; ig < num_groups; ig++)
    if (functors[ig])
    {
      mw_coefs_ptr[ig]         = functors[ig]->spline_coefs_->device_data();
      mw_DeltaRInv_ptr[ig]     = functors[ig]->DeltaRInv;
      mw_cutoff_radius_ptr[ig] = functors[ig]->cutoff_radius;
    }
    else
    {
      mw_coefs_ptr[ig]         = nullptr;
      mw_DeltaRInv_ptr[ig]     = 0.0;
      mw_cutoff_radius_ptr[ig] = 0.0; // important! Prevent spline evaluation to access nullptr.
    }

  int nw_accepted = 0;
  for (int iw = 0; iw < nw; iw++)
    if (isAccepted[iw])
      accepted_indices[nw_accepted++] = iw;

  auto* transfer_buffer_ptr = transfer_buffer.data();

  PRAGMA_OFFLOAD("omp target teams distribute map(always, to: transfer_buffer_ptr[:transfer_buffer.size()]) \
                    map(to: grp_ids[:n_src]) \
                    map(to: mw_dist[:dist_stride*nw]) \
                    map(to: mw_vgl[:(DIM+2)*nw]) \
                    map(always, from: mw_allUat[:nw * n_padded * (DIM + 2)])")
  for (int iw = 0; iw < nw_accepted; iw++)
  {
    REAL** mw_coefs        = reinterpret_cast<REAL**>(transfer_buffer_ptr);
    REAL* mw_DeltaRInv     = reinterpret_cast<REAL*>(transfer_buffer_ptr + sizeof(REAL*) * num_groups);
    REAL* mw_cutoff_radius = mw_DeltaRInv + num_groups;
    int* accepted_indices =
        reinterpret_cast<int*>(transfer_buffer_ptr + (sizeof(REAL*) + sizeof(REAL) * 2) * num_groups);
    const int ip = accepted_indices[iw];

    const REAL* dist_new   = mw_dist + ip * dist_stride;
    const REAL* dipl_x_new = dist_new + n_padded;
    const REAL* dipl_y_new = dist_new + n_padded * 2;
    const REAL* dipl_z_new = dist_new + n_padded * 3;

    const REAL* dist_old   = mw_dist + ip * dist_stride + dist_stride * nw;
    const REAL* dipl_x_old = dist_old + n_padded;
    const REAL* dipl_y_old = dist_old + n_padded * 2;
    const REAL* dipl_z_old = dist_old + n_padded * 3;

    REAL* Uat    = mw_allUat + ip * n_padded;
    REAL* dUat_x = mw_allUat + n_padded * nw + ip * n_padded * DIM;
    REAL* dUat_y = dUat_x + n_padded;
    REAL* dUat_z = dUat_y + n_padded;
    REAL* d2Uat  = mw_allUat + n_padded * (DIM + 1) * nw + ip * n_padded;

    REAL* cur_allu = mw_cur_allu + ip * n_padded * 3;

    PRAGMA_OFFLOAD("omp parallel for")
    for (int j = 0; j < n_src; j++)
    {
      if (j == iat)
        continue;
      const int ig       = grp_ids[j];
      const REAL* coefs  = mw_coefs[ig];
      REAL DeltaRInv     = mw_DeltaRInv[ig];
      REAL cutoff_radius = mw_cutoff_radius[ig];

      REAL r = dist_old[j];
      REAL u(0);
      REAL dudr(0);
      REAL d2udr2(0);
      if (r < cutoff_radius)
      {
        u = evaluate_impl(dist_old[j], coefs, DeltaRInv, dudr, d2udr2);
        dudr *= REAL(1) / r;
      }
      // update Uat, dUat, d2Uat
      REAL cur_u      = cur_allu[j];
      REAL cur_dudr   = cur_allu[j + n_padded];
      REAL cur_d2udr2 = cur_allu[j + n_padded * 2];
      Uat[j] += cur_u - u;
      dUat_x[j] -= dipl_x_new[j] * cur_dudr - dipl_x_old[j] * dudr;
      dUat_y[j] -= dipl_y_new[j] * cur_dudr - dipl_y_old[j] * dudr;
      dUat_z[j] -= dipl_z_new[j] * cur_dudr - dipl_z_old[j] * dudr;
      constexpr REAL lapfac(DIM - 1);
      d2Uat[j] -= cur_d2udr2 + lapfac * cur_dudr - (d2udr2 + lapfac * dudr);
    }
    REAL* vgl   = mw_vgl + ip * (DIM + 2);
    Uat[iat]    = vgl[0];
    dUat_x[iat] = vgl[1];
    dUat_y[iat] = vgl[2];
    dUat_z[iat] = vgl[3];
    d2Uat[iat]  = vgl[4];
  }
}

template struct BsplineFunctor<QMCTraits::RealType>;

} // namespace qmcplusplus
