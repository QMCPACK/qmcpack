//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "SplineC2ROMPTarget.h"
#include "spline2/MultiBsplineEval.hpp"
#include "spline2/MultiBsplineEval_OMPoffload.hpp"
#include "QMCWaveFunctions/BsplineFactory/contraction_helper.hpp"
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
                     size_t orb_size,
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
  const size_t psiIndex = first_spo + index + (index < nComplexBands ? index : nComplexBands);
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
                       const ST* mKK_ptr,
                       size_t orb_size,
                       const ST* restrict offload_scratch_ptr,
                       size_t spline_padded_size,
                       const ST symGGt[6],
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

  const ST* restrict val = offload_scratch_ptr;
  const ST* restrict g0  = offload_scratch_ptr + spline_padded_size;
  const ST* restrict g1  = offload_scratch_ptr + spline_padded_size * 2;
  const ST* restrict g2  = offload_scratch_ptr + spline_padded_size * 3;
  const ST* restrict h00 = offload_scratch_ptr + spline_padded_size * 4;
  const ST* restrict h01 = offload_scratch_ptr + spline_padded_size * 5;
  const ST* restrict h02 = offload_scratch_ptr + spline_padded_size * 6;
  const ST* restrict h11 = offload_scratch_ptr + spline_padded_size * 7;
  const ST* restrict h12 = offload_scratch_ptr + spline_padded_size * 8;
  const ST* restrict h22 = offload_scratch_ptr + spline_padded_size * 9;

  TT* restrict psi   = results_scratch_ptr;
  TT* restrict dpsi  = results_scratch_ptr + orb_size;
  TT* restrict d2psi = results_scratch_ptr + orb_size * 4;

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

  const ST lcart_r = SymTrace(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], symGGt);
  const ST lcart_i = SymTrace(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], symGGt);
  const ST lap_r   = lcart_r + mKK_ptr[index] * val_r + two * (kX * dX_i + kY * dY_i + kZ * dZ_i);
  const ST lap_i   = lcart_i + mKK_ptr[index] * val_i - two * (kX * dX_r + kY * dY_r + kZ * dZ_r);

  const size_t psiIndex  = first_spo + index + (index < nComplexBands ? index : nComplexBands);
  psi[psiIndex]          = c * val_r - s * val_i;
  d2psi[psiIndex]        = c * lap_r - s * lap_i;
  dpsi[psiIndex * 3]     = c * gX_r - s * gX_i;
  dpsi[psiIndex * 3 + 1] = c * gY_r - s * gY_i;
  dpsi[psiIndex * 3 + 2] = c * gZ_r - s * gZ_i;

  if (index < nComplexBands)
  {
    psi[psiIndex + 1]      = c * val_i + s * val_r;
    d2psi[psiIndex + 1]    = c * lap_i + s * lap_r;
    dpsi[psiIndex * 3 + 3] = c * gX_i + s * gX_r;
    dpsi[psiIndex * 3 + 4] = c * gY_i + s * gY_r;
    dpsi[psiIndex * 3 + 5] = c * gZ_i + s * gZ_r;
  }
}
} // namespace C2R

template<typename ST>
inline void SplineC2ROMPTarget<ST>::set_spline(SingleSplineType* spline_r,
                                               SingleSplineType* spline_i,
                                               int twist,
                                               int ispline,
                                               int level)
{
  SplineInst->copy_spline(spline_r, 2 * ispline);
  SplineInst->copy_spline(spline_i, 2 * ispline + 1);
}

template<typename ST>
bool SplineC2ROMPTarget<ST>::read_splines(hdf_archive& h5f)
{
  std::ostringstream o;
  o << "spline_" << MyIndex;
  einspline_engine<SplineType> bigtable(SplineInst->getSplinePtr());
  return h5f.readEntry(bigtable, o.str().c_str()); //"spline_0");
}

template<typename ST>
bool SplineC2ROMPTarget<ST>::write_splines(hdf_archive& h5f)
{
  std::ostringstream o;
  o << "spline_" << MyIndex;
  einspline_engine<SplineType> bigtable(SplineInst->getSplinePtr());
  return h5f.writeEntry(bigtable, o.str().c_str()); //"spline_0");
}

template<typename ST>
inline void SplineC2ROMPTarget<ST>::assign_v(const PointType& r,
                                             const vContainer_type& myV,
                                             ValueVector& psi,
                                             int first,
                                             int last) const
{
  // protect last
  last = last > kPoints.size() ? kPoints.size() : last;

  const ST x = r[0], y = r[1], z = r[2];
  const ST* restrict kx = myKcart->data(0);
  const ST* restrict ky = myKcart->data(1);
  const ST* restrict kz = myKcart->data(2);

  TT* restrict psi_s = psi.data() + first_spo;
#pragma omp simd
  for (size_t j = first; j < std::min(nComplexBands, last); j++)
  {
    ST s, c;
    const size_t jr = j << 1;
    const size_t ji = jr + 1;
    const ST val_r  = myV[jr];
    const ST val_i  = myV[ji];
    omptarget::sincos(-(x * kx[j] + y * ky[j] + z * kz[j]), &s, &c);
    psi_s[jr] = val_r * c - val_i * s;
    psi_s[ji] = val_i * c + val_r * s;
  }

  psi_s += nComplexBands;
#pragma omp simd
  for (size_t j = std::max(nComplexBands, first); j < last; j++)
  {
    ST s, c;
    const ST val_r = myV[2 * j];
    const ST val_i = myV[2 * j + 1];
    omptarget::sincos(-(x * kx[j] + y * ky[j] + z * kz[j]), &s, &c);
    psi_s[j] = val_r * c - val_i * s;
  }
}

template<typename ST>
void SplineC2ROMPTarget<ST>::evaluateValue(const ParticleSet& P, const int iat, ValueVector& psi)
{
  const PointType& r = P.activeR(iat);
  PointType ru(PrimLattice.toUnit_floor(r));

  if (true)
  {
#pragma omp parallel
    {
      int first, last;
      FairDivideAligned(myV.size(), getAlignment<ST>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

      spline2::evaluate3d(SplineInst->getSplinePtr(), ru, myV, first, last);
      assign_v(r, myV, psi, first / 2, last / 2);
    }
  }
  else
  {
    const int ChunkSizePerTeam = 128;
    const int NumTeams         = (myV.size() + ChunkSizePerTeam - 1) / ChunkSizePerTeam;

    const auto padded_size = myV.size();
    offload_scratch.resize(padded_size);

    // Ye: need to extract sizes and pointers before entering target region
    const auto orb_size       = psi.size();
    const auto* spline_ptr    = SplineInst->getSplinePtr();
    auto* offload_scratch_ptr = offload_scratch.data();
    auto* psi_ptr             = psi.data();
    const auto x = r[0], y = r[1], z = r[2];
    const auto rux = ru[0], ruy = ru[1], ruz = ru[2];
    const auto myKcart_padded_size = myKcart->capacity();
    auto* myKcart_ptr              = myKcart->data();
    const size_t first_spo_local   = first_spo;
    const int nComplexBands_local  = nComplexBands;

    {
      ScopedTimer offload(offload_timer_);
      PRAGMA_OFFLOAD("omp target teams distribute num_teams(NumTeams) \
                  map(always, from: psi_ptr[0:orb_size])")
      for (int team_id = 0; team_id < NumTeams; team_id++)
      {
        const int first = ChunkSizePerTeam * team_id;
        const int last  = (first + ChunkSizePerTeam) > padded_size ? padded_size : first + ChunkSizePerTeam;

        int ix, iy, iz;
        ST a[4], b[4], c[4];
        spline2::computeLocationAndFractional(spline_ptr, rux, ruy, ruz, ix, iy, iz, a, b, c);

        PRAGMA_OFFLOAD("omp parallel for")
        for (int index = 0; index < last - first; index++)
          spline2offload::evaluate_v_impl_v2(spline_ptr, ix, iy, iz, a, b, c, offload_scratch_ptr + first, first,
                                             index);
        const int first_cplx = first / 2;
        const int last_cplx  = last / 2 < orb_size ? last / 2 : orb_size;
        PRAGMA_OFFLOAD("omp parallel for")
        for (int index = first_cplx; index < last_cplx; index++)
          C2R::assign_v(x, y, z, psi_ptr, orb_size, offload_scratch_ptr, myKcart_ptr, myKcart_padded_size,
                        first_spo_local, nComplexBands_local, index);
      }
    }
  }
}

template<typename ST>
void SplineC2ROMPTarget<ST>::evaluateDetRatios(const VirtualParticleSet& VP,
                                               ValueVector& psi,
                                               const ValueVector& psiinv,
                                               std::vector<ValueType>& ratios)
{
  const int nVP = VP.getTotalNum();
  psiinv_pos_copy.resize(psiinv.size() + nVP * 6);

  // stage psiinv to psiinv_pos_copy
  std::copy_n(psiinv.data(), psiinv.size(), psiinv_pos_copy.data());

  // pack particle positions
  auto* restrict pos_scratch = psiinv_pos_copy.data() + psiinv.size();
  for (int iat = 0; iat < nVP; ++iat)
  {
    const PointType& r = VP.activeR(iat);
    PointType ru(PrimLattice.toUnit_floor(r));
    pos_scratch[iat * 6]     = r[0];
    pos_scratch[iat * 6 + 1] = r[1];
    pos_scratch[iat * 6 + 2] = r[2];
    pos_scratch[iat * 6 + 3] = ru[0];
    pos_scratch[iat * 6 + 4] = ru[1];
    pos_scratch[iat * 6 + 5] = ru[2];
  }

  const int ChunkSizePerTeam = 128;
  const int NumTeams         = (myV.size() + ChunkSizePerTeam - 1) / ChunkSizePerTeam;
  ratios_private.resize(nVP, NumTeams);
  const auto padded_size = myV.size();
  offload_scratch.resize(padded_size * nVP);
  const auto orb_size = psiinv.size();
  results_scratch.resize(orb_size * nVP);

  // Ye: need to extract sizes and pointers before entering target region
  const auto* spline_ptr         = SplineInst->getSplinePtr();
  auto* offload_scratch_ptr      = offload_scratch.data();
  auto* results_scratch_ptr      = results_scratch.data();
  const auto myKcart_padded_size = myKcart->capacity();
  auto* myKcart_ptr              = myKcart->data();
  auto* psiinv_ptr               = psiinv_pos_copy.data();
  auto* ratios_private_ptr       = ratios_private.data();
  const size_t first_spo_local   = first_spo;
  const int nComplexBands_local  = nComplexBands;

  {
    ScopedTimer offload(offload_timer_);
    PRAGMA_OFFLOAD("omp target teams distribute collapse(2) num_teams(NumTeams*nVP) \
                map(always, to: psiinv_ptr[0:psiinv_pos_copy.size()]) \
                map(always, from: ratios_private_ptr[0:NumTeams*nVP])")
    for (int iat = 0; iat < nVP; iat++)
      for (int team_id = 0; team_id < NumTeams; team_id++)
      {
        const int first = ChunkSizePerTeam * team_id;
        const int last  = (first + ChunkSizePerTeam) > padded_size ? padded_size : first + ChunkSizePerTeam;
        auto* restrict offload_scratch_iat_ptr = offload_scratch_ptr + padded_size * iat;
        auto* restrict psi_iat_ptr             = results_scratch_ptr + orb_size * iat;
        auto* restrict pos_scratch             = psiinv_ptr + orb_size;

        int ix, iy, iz;
        ST a[4], b[4], c[4];
        spline2::computeLocationAndFractional(spline_ptr, ST(pos_scratch[iat * 6 + 3]), ST(pos_scratch[iat * 6 + 4]),
                                              ST(pos_scratch[iat * 6 + 5]), ix, iy, iz, a, b, c);

        PRAGMA_OFFLOAD("omp parallel for")
        for (int index = 0; index < last - first; index++)
          spline2offload::evaluate_v_impl_v2(spline_ptr, ix, iy, iz, a, b, c, offload_scratch_iat_ptr + first, first,
                                             index);
        const int first_cplx = first / 2;
        const int last_cplx  = last / 2 < orb_size ? last / 2 : orb_size;
        PRAGMA_OFFLOAD("omp parallel for")
        for (int index = first_cplx; index < last_cplx; index++)
          C2R::assign_v(ST(pos_scratch[iat * 6]), ST(pos_scratch[iat * 6 + 1]), ST(pos_scratch[iat * 6 + 2]),
                        psi_iat_ptr, orb_size, offload_scratch_iat_ptr, myKcart_ptr, myKcart_padded_size,
                        first_spo_local, nComplexBands_local, index);

        const int first_real = first_cplx + std::min(nComplexBands_local, first_cplx);
        const int last_real  = last_cplx + std::min(nComplexBands_local, last_cplx);
        TT sum(0);
        PRAGMA_OFFLOAD("omp parallel for simd reduction(+:sum)")
        for (int i = first_real; i < last_real; i++)
          sum += psi_iat_ptr[i] * psiinv_ptr[i];
        ratios_private_ptr[iat * NumTeams + team_id] = sum;
      }
  }

  // do the reduction manually
  for (int iat = 0; iat < nVP; ++iat)
  {
    ratios[iat] = TT(0);
    for (int tid = 0; tid < NumTeams; tid++)
      ratios[iat] += ratios_private[iat][tid];
  }
}

template<typename ST>
void SplineC2ROMPTarget<ST>::mw_evaluateDetRatios(const RefVectorWithLeader<SPOSet>& spo_list,
                                                  const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                                                  const RefVector<ValueVector>& psi_list,
                                                  const std::vector<const ValueType*>& invRow_ptr_list,
                                                  std::vector<std::vector<ValueType>>& ratios_list) const
{
  assert(this == &spo_list.getLeader());
  auto& phi_leader            = spo_list.getCastedLeader<SplineC2ROMPTarget<ST>>();
  auto& mw_mem                = *phi_leader.mw_mem_;
  auto& det_ratios_buffer_H2D = mw_mem.det_ratios_buffer_H2D;
  auto& mw_ratios_private     = mw_mem.mw_ratios_private;
  auto& mw_offload_scratch    = mw_mem.mw_offload_scratch;
  auto& mw_results_scratch    = mw_mem.mw_results_scratch;
  const size_t nw             = spo_list.size();
  const size_t orb_size       = phi_leader.size();

  size_t mw_nVP = 0;
  for (const VirtualParticleSet& VP : vp_list)
    mw_nVP += VP.getTotalNum();

  const size_t packed_size = nw * sizeof(ValueType*) + mw_nVP * (6 * sizeof(TT) + sizeof(int));
  det_ratios_buffer_H2D.resize(packed_size);

  // pack invRow_ptr_list to det_ratios_buffer_H2D
  Vector<const ValueType*> ptr_buffer(reinterpret_cast<const ValueType**>(det_ratios_buffer_H2D.data()), nw);
  for (size_t iw = 0; iw < nw; iw++)
    ptr_buffer[iw] = invRow_ptr_list[iw];

  // pack particle positions
  auto* pos_ptr = reinterpret_cast<TT*>(det_ratios_buffer_H2D.data() + nw * sizeof(ValueType*));
  auto* ref_id_ptr =
      reinterpret_cast<int*>(det_ratios_buffer_H2D.data() + nw * sizeof(ValueType*) + mw_nVP * 6 * sizeof(TT));
  size_t iVP = 0;
  for (size_t iw = 0; iw < nw; iw++)
  {
    const VirtualParticleSet& VP = vp_list[iw];
    assert(ratios_list[iw].size() == VP.getTotalNum());
    for (size_t iat = 0; iat < VP.getTotalNum(); ++iat, ++iVP)
    {
      ref_id_ptr[iVP]    = iw;
      const PointType& r = VP.activeR(iat);
      PointType ru(PrimLattice.toUnit_floor(r));
      pos_ptr[0] = r[0];
      pos_ptr[1] = r[1];
      pos_ptr[2] = r[2];
      pos_ptr[3] = ru[0];
      pos_ptr[4] = ru[1];
      pos_ptr[5] = ru[2];
      pos_ptr += 6;
    }
  }

  const int ChunkSizePerTeam = 128;
  const int NumTeams         = (myV.size() + ChunkSizePerTeam - 1) / ChunkSizePerTeam;
  mw_ratios_private.resize(mw_nVP, NumTeams);
  const auto padded_size = myV.size();
  mw_offload_scratch.resize(padded_size * mw_nVP);
  mw_results_scratch.resize(orb_size * mw_nVP);

  // Ye: need to extract sizes and pointers before entering target region
  const auto* spline_ptr         = SplineInst->getSplinePtr();
  auto* offload_scratch_ptr      = mw_offload_scratch.data();
  auto* results_scratch_ptr      = mw_results_scratch.data();
  const auto myKcart_padded_size = myKcart->capacity();
  auto* myKcart_ptr              = myKcart->data();
  auto* buffer_H2D_ptr           = det_ratios_buffer_H2D.data();
  auto* ratios_private_ptr       = mw_ratios_private.data();
  const size_t first_spo_local   = first_spo;
  const int nComplexBands_local  = nComplexBands;

  {
    ScopedTimer offload(offload_timer_);
    PRAGMA_OFFLOAD("omp target teams distribute collapse(2) num_teams(NumTeams*mw_nVP) \
                map(always, to: buffer_H2D_ptr[0:det_ratios_buffer_H2D.size()]) \
                map(always, from: ratios_private_ptr[0:NumTeams*mw_nVP])")
    for (int iat = 0; iat < mw_nVP; iat++)
      for (int team_id = 0; team_id < NumTeams; team_id++)
      {
        const int first = ChunkSizePerTeam * team_id;
        const int last  = (first + ChunkSizePerTeam) > padded_size ? padded_size : first + ChunkSizePerTeam;
        auto* restrict offload_scratch_iat_ptr = offload_scratch_ptr + padded_size * iat;
        auto* restrict psi_iat_ptr             = results_scratch_ptr + orb_size * iat;
        auto* ref_id_ptr = reinterpret_cast<int*>(buffer_H2D_ptr + nw * sizeof(ValueType*) + mw_nVP * 6 * sizeof(TT));
        auto* restrict psiinv_ptr  = reinterpret_cast<const ValueType**>(buffer_H2D_ptr)[ref_id_ptr[iat]];
        auto* restrict pos_scratch = reinterpret_cast<TT*>(buffer_H2D_ptr + nw * sizeof(ValueType*));

        int ix, iy, iz;
        ST a[4], b[4], c[4];
        spline2::computeLocationAndFractional(spline_ptr, ST(pos_scratch[iat * 6 + 3]), ST(pos_scratch[iat * 6 + 4]),
                                              ST(pos_scratch[iat * 6 + 5]), ix, iy, iz, a, b, c);

        PRAGMA_OFFLOAD("omp parallel for")
        for (int index = 0; index < last - first; index++)
          spline2offload::evaluate_v_impl_v2(spline_ptr, ix, iy, iz, a, b, c, offload_scratch_iat_ptr + first, first,
                                             index);
        const int first_cplx = first / 2;
        const int last_cplx  = last / 2 < orb_size ? last / 2 : orb_size;
        PRAGMA_OFFLOAD("omp parallel for")
        for (int index = first_cplx; index < last_cplx; index++)
          C2R::assign_v(ST(pos_scratch[iat * 6]), ST(pos_scratch[iat * 6 + 1]), ST(pos_scratch[iat * 6 + 2]),
                        psi_iat_ptr, orb_size, offload_scratch_iat_ptr, myKcart_ptr, myKcart_padded_size,
                        first_spo_local, nComplexBands_local, index);

        const int first_real = first_cplx + std::min(nComplexBands_local, first_cplx);
        const int last_real  = last_cplx + std::min(nComplexBands_local, last_cplx);
        TT sum(0);
        PRAGMA_OFFLOAD("omp parallel for simd reduction(+:sum)")
        for (int i = first_real; i < last_real; i++)
          sum += psi_iat_ptr[i] * psiinv_ptr[i];
        ratios_private_ptr[iat * NumTeams + team_id] = sum;
      }
  }

  // do the reduction manually
  iVP = 0;
  for (size_t iw = 0; iw < nw; iw++)
  {
    auto& ratios = ratios_list[iw];
    for (size_t iat = 0; iat < ratios.size(); iat++, iVP++)
    {
      ratios[iat] = TT(0);
      for (int tid = 0; tid < NumTeams; ++tid)
        ratios[iat] += mw_ratios_private[iVP][tid];
    }
  }
}

/** assign_vgl_from_l can be used when myL is precomputed and myV,myG,myL in cartesian
   */
template<typename ST>
inline void SplineC2ROMPTarget<ST>::assign_vgl_from_l(const PointType& r,
                                                      ValueVector& psi,
                                                      GradVector& dpsi,
                                                      ValueVector& d2psi)
{
  constexpr ST two(2);
  const ST x = r[0], y = r[1], z = r[2];

  const ST* restrict k0 = myKcart->data(0);
  ASSUME_ALIGNED(k0);
  const ST* restrict k1 = myKcart->data(1);
  ASSUME_ALIGNED(k1);
  const ST* restrict k2 = myKcart->data(2);
  ASSUME_ALIGNED(k2);

  const ST* restrict g0 = myG.data(0);
  ASSUME_ALIGNED(g0);
  const ST* restrict g1 = myG.data(1);
  ASSUME_ALIGNED(g1);
  const ST* restrict g2 = myG.data(2);
  ASSUME_ALIGNED(g2);

  const size_t N = kPoints.size();

#pragma omp simd
  for (size_t j = 0; j < nComplexBands; j++)
  {
    const size_t jr = j << 1;
    const size_t ji = jr + 1;

    const ST kX    = k0[j];
    const ST kY    = k1[j];
    const ST kZ    = k2[j];
    const ST val_r = myV[jr];
    const ST val_i = myV[ji];

    //phase
    ST s, c;
    omptarget::sincos(-(x * kX + y * kY + z * kZ), &s, &c);

    //dot(PrimLattice.G,myG[j])
    const ST dX_r = g0[jr];
    const ST dY_r = g1[jr];
    const ST dZ_r = g2[jr];

    const ST dX_i = g0[ji];
    const ST dY_i = g1[ji];
    const ST dZ_i = g2[ji];

    // \f$\nabla \psi_r + {\bf k}\psi_i\f$
    const ST gX_r = dX_r + val_i * kX;
    const ST gY_r = dY_r + val_i * kY;
    const ST gZ_r = dZ_r + val_i * kZ;
    const ST gX_i = dX_i - val_r * kX;
    const ST gY_i = dY_i - val_r * kY;
    const ST gZ_i = dZ_i - val_r * kZ;

    const ST lap_r = myL[jr] + (*mKK)[j] * val_r + two * (kX * dX_i + kY * dY_i + kZ * dZ_i);
    const ST lap_i = myL[ji] + (*mKK)[j] * val_i - two * (kX * dX_r + kY * dY_r + kZ * dZ_r);

    const size_t psiIndex = first_spo + jr;
    psi[psiIndex]         = c * val_r - s * val_i;
    psi[psiIndex + 1]     = c * val_i + s * val_r;
    d2psi[psiIndex]       = c * lap_r - s * lap_i;
    d2psi[psiIndex + 1]   = c * lap_i + s * lap_r;
    dpsi[psiIndex][0]     = c * gX_r - s * gX_i;
    dpsi[psiIndex][1]     = c * gY_r - s * gY_i;
    dpsi[psiIndex][2]     = c * gZ_r - s * gZ_i;
    dpsi[psiIndex + 1][0] = c * gX_i + s * gX_r;
    dpsi[psiIndex + 1][1] = c * gY_i + s * gY_r;
    dpsi[psiIndex + 1][2] = c * gZ_i + s * gZ_r;
  }

#pragma omp simd
  for (size_t j = nComplexBands; j < N; j++)
  {
    const size_t jr = j << 1;
    const size_t ji = jr + 1;

    const ST kX    = k0[j];
    const ST kY    = k1[j];
    const ST kZ    = k2[j];
    const ST val_r = myV[jr];
    const ST val_i = myV[ji];

    //phase
    ST s, c;
    omptarget::sincos(-(x * kX + y * kY + z * kZ), &s, &c);

    //dot(PrimLattice.G,myG[j])
    const ST dX_r = g0[jr];
    const ST dY_r = g1[jr];
    const ST dZ_r = g2[jr];

    const ST dX_i = g0[ji];
    const ST dY_i = g1[ji];
    const ST dZ_i = g2[ji];

    // \f$\nabla \psi_r + {\bf k}\psi_i\f$
    const ST gX_r         = dX_r + val_i * kX;
    const ST gY_r         = dY_r + val_i * kY;
    const ST gZ_r         = dZ_r + val_i * kZ;
    const ST gX_i         = dX_i - val_r * kX;
    const ST gY_i         = dY_i - val_r * kY;
    const ST gZ_i         = dZ_i - val_r * kZ;
    const size_t psiIndex = first_spo + nComplexBands + j;
    psi[psiIndex]         = c * val_r - s * val_i;
    dpsi[psiIndex][0]     = c * gX_r - s * gX_i;
    dpsi[psiIndex][1]     = c * gY_r - s * gY_i;
    dpsi[psiIndex][2]     = c * gZ_r - s * gZ_i;

    const ST lap_r  = myL[jr] + (*mKK)[j] * val_r + two * (kX * dX_i + kY * dY_i + kZ * dZ_i);
    const ST lap_i  = myL[ji] + (*mKK)[j] * val_i - two * (kX * dX_r + kY * dY_r + kZ * dZ_r);
    d2psi[psiIndex] = c * lap_r - s * lap_i;
  }
}

template<typename ST>
void SplineC2ROMPTarget<ST>::evaluateVGL(const ParticleSet& P,
                                         const int iat,
                                         ValueVector& psi,
                                         GradVector& dpsi,
                                         ValueVector& d2psi)
{
  const PointType& r = P.activeR(iat);
  PointType ru(PrimLattice.toUnit_floor(r));

  const int ChunkSizePerTeam = 128;
  const int NumTeams         = (myV.size() + ChunkSizePerTeam - 1) / ChunkSizePerTeam;

  const auto padded_size = myV.size();
  // for V(1)G(3)H(6) intermediate result
  offload_scratch.resize(padded_size * 10);
  const auto orb_size = psi.size();
  // for V(1)G(3)L(1) final result
  results_scratch.resize(orb_size * 5);

  // Ye: need to extract sizes and pointers before entering target region
  const auto* spline_ptr    = SplineInst->getSplinePtr();
  auto* offload_scratch_ptr = offload_scratch.data();
  auto* results_scratch_ptr = results_scratch.data();
  const auto x = r[0], y = r[1], z = r[2];
  const auto rux = ru[0], ruy = ru[1], ruz = ru[2];
  const auto myKcart_padded_size = myKcart->capacity();
  auto* mKK_ptr                  = mKK->data();
  auto* GGt_ptr                  = GGt_offload->data();
  auto* PrimLattice_G_ptr        = PrimLattice_G_offload->data();
  auto* myKcart_ptr              = myKcart->data();
  const size_t first_spo_local   = first_spo;
  const int nComplexBands_local  = nComplexBands;

  {
    ScopedTimer offload(offload_timer_);
    PRAGMA_OFFLOAD("omp target teams distribute num_teams(NumTeams) \
                map(always, from: results_scratch_ptr[0:orb_size*5])")
    for (int team_id = 0; team_id < NumTeams; team_id++)
    {
      const int first = ChunkSizePerTeam * team_id;
      const int last  = (first + ChunkSizePerTeam) > padded_size ? padded_size : first + ChunkSizePerTeam;

      int ix, iy, iz;
      ST a[4], b[4], c[4], da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
      spline2::computeLocationAndFractional(spline_ptr, rux, ruy, ruz, ix, iy, iz, a, b, c, da, db, dc, d2a, d2b, d2c);

      const ST G[9]      = {PrimLattice_G_ptr[0], PrimLattice_G_ptr[1], PrimLattice_G_ptr[2],
                            PrimLattice_G_ptr[3], PrimLattice_G_ptr[4], PrimLattice_G_ptr[5],
                            PrimLattice_G_ptr[6], PrimLattice_G_ptr[7], PrimLattice_G_ptr[8]};
      const ST symGGt[6] = {GGt_ptr[0], GGt_ptr[1] + GGt_ptr[3], GGt_ptr[2] + GGt_ptr[6],
                            GGt_ptr[4], GGt_ptr[5] + GGt_ptr[7], GGt_ptr[8]};

      PRAGMA_OFFLOAD("omp parallel for")
      for (int index = 0; index < last - first; index++)
        spline2offload::evaluate_vgh_impl_v2(spline_ptr, ix, iy, iz, a, b, c, da, db, dc, d2a, d2b, d2c,
                                             offload_scratch_ptr + first, offload_scratch_ptr + padded_size + first,
                                             offload_scratch_ptr + padded_size * 4 + first, padded_size, first, index);
      const int first_cplx = first / 2;
      const int last_cplx  = last / 2 < orb_size ? last / 2 : orb_size;
      PRAGMA_OFFLOAD("omp parallel for")
      for (int index = first_cplx; index < last_cplx; index++)
        C2R::assign_vgl(x, y, z, results_scratch_ptr, mKK_ptr, orb_size, offload_scratch_ptr, padded_size, symGGt, G,
                        myKcart_ptr, myKcart_padded_size, first_spo_local, nComplexBands_local, index);
    }
  }

  for (size_t i = 0; i < orb_size; i++)
  {
    psi[i]     = results_scratch[i];
    dpsi[i][0] = results_scratch[orb_size + i * 3];
    dpsi[i][1] = results_scratch[orb_size + i * 3 + 1];
    dpsi[i][2] = results_scratch[orb_size + i * 3 + 2];
    d2psi[i]   = results_scratch[orb_size * 4 + i];
  }
}

template<typename ST>
void SplineC2ROMPTarget<ST>::evaluateVGLMultiPos(const Vector<ST, OffloadPinnedAllocator<ST>>& multi_pos,
                                                 Vector<ST, OffloadPinnedAllocator<ST>>& offload_scratch,
                                                 Vector<TT, OffloadPinnedAllocator<TT>>& results_scratch,
                                                 const RefVector<ValueVector>& psi_v_list,
                                                 const RefVector<GradVector>& dpsi_v_list,
                                                 const RefVector<ValueVector>& d2psi_v_list) const
{
  const size_t num_pos       = psi_v_list.size();
  const int ChunkSizePerTeam = 128;
  const int NumTeams         = (myV.size() + ChunkSizePerTeam - 1) / ChunkSizePerTeam;
  const auto padded_size     = myV.size();
  // for V(1)G(3)H(6) intermediate result
  offload_scratch.resize(padded_size * num_pos * 10);
  const auto orb_size = psi_v_list[0].get().size();
  // for V(1)G(3)L(1) final result
  results_scratch.resize(orb_size * num_pos * 5);

  // Ye: need to extract sizes and pointers before entering target region
  const auto* spline_ptr         = SplineInst->getSplinePtr();
  auto* pos_copy_ptr             = multi_pos.data();
  auto* offload_scratch_ptr      = offload_scratch.data();
  auto* results_scratch_ptr      = results_scratch.data();
  const auto myKcart_padded_size = myKcart->capacity();
  auto* mKK_ptr                  = mKK->data();
  auto* GGt_ptr                  = GGt_offload->data();
  auto* PrimLattice_G_ptr        = PrimLattice_G_offload->data();
  auto* myKcart_ptr              = myKcart->data();
  const size_t first_spo_local   = first_spo;
  const int nComplexBands_local  = nComplexBands;

  {
    ScopedTimer offload(offload_timer_);
    PRAGMA_OFFLOAD("omp target teams distribute collapse(2) num_teams(NumTeams*num_pos) \
                    map(always, to: pos_copy_ptr[0:num_pos*6]) \
                    map(always, from: results_scratch_ptr[0:orb_size*num_pos*5])")
    for (int iw = 0; iw < num_pos; iw++)
      for (int team_id = 0; team_id < NumTeams; team_id++)
      {
        const int first = ChunkSizePerTeam * team_id;
        const int last  = (first + ChunkSizePerTeam) > padded_size ? padded_size : first + ChunkSizePerTeam;
        auto* restrict offload_scratch_iw_ptr = offload_scratch_ptr + padded_size * iw * 10;
        auto* restrict psi_iw_ptr             = results_scratch_ptr + orb_size * iw * 5;

        int ix, iy, iz;
        ST a[4], b[4], c[4], da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
        spline2::computeLocationAndFractional(spline_ptr, pos_copy_ptr[iw * 6 + 3], pos_copy_ptr[iw * 6 + 4],
                                              pos_copy_ptr[iw * 6 + 5], ix, iy, iz, a, b, c, da, db, dc, d2a, d2b, d2c);

        const ST G[9]      = {PrimLattice_G_ptr[0], PrimLattice_G_ptr[1], PrimLattice_G_ptr[2],
                              PrimLattice_G_ptr[3], PrimLattice_G_ptr[4], PrimLattice_G_ptr[5],
                              PrimLattice_G_ptr[6], PrimLattice_G_ptr[7], PrimLattice_G_ptr[8]};
        const ST symGGt[6] = {GGt_ptr[0], GGt_ptr[1] + GGt_ptr[3], GGt_ptr[2] + GGt_ptr[6],
                              GGt_ptr[4], GGt_ptr[5] + GGt_ptr[7], GGt_ptr[8]};

        PRAGMA_OFFLOAD("omp parallel for")
        for (int index = 0; index < last - first; index++)
          spline2offload::evaluate_vgh_impl_v2(spline_ptr, ix, iy, iz, a, b, c, da, db, dc, d2a, d2b, d2c,
                                               offload_scratch_iw_ptr + first,
                                               offload_scratch_iw_ptr + padded_size + first,
                                               offload_scratch_iw_ptr + padded_size * 4 + first, padded_size, first,
                                               index);
        const int first_cplx = first / 2;
        const int last_cplx  = last / 2 < orb_size ? last / 2 : orb_size;
        PRAGMA_OFFLOAD("omp parallel for")
        for (int index = first_cplx; index < last_cplx; index++)
          C2R::assign_vgl(pos_copy_ptr[iw * 6], pos_copy_ptr[iw * 6 + 1], pos_copy_ptr[iw * 6 + 2], psi_iw_ptr, mKK_ptr,
                          orb_size, offload_scratch_iw_ptr, padded_size, symGGt, G, myKcart_ptr, myKcart_padded_size,
                          first_spo_local, nComplexBands_local, index);
      }
  }

  for (int iw = 0; iw < num_pos; ++iw)
  {
    auto* restrict results_iw_ptr = results_scratch_ptr + orb_size * iw * 5;
    ValueVector& psi_v(psi_v_list[iw]);
    GradVector& dpsi_v(dpsi_v_list[iw]);
    ValueVector& d2psi_v(d2psi_v_list[iw]);
    for (size_t i = 0; i < orb_size; i++)
    {
      psi_v[i]     = results_iw_ptr[i];
      dpsi_v[i][0] = results_iw_ptr[orb_size + i * 3];
      dpsi_v[i][1] = results_iw_ptr[orb_size + i * 3 + 1];
      dpsi_v[i][2] = results_iw_ptr[orb_size + i * 3 + 2];
      d2psi_v[i]   = results_iw_ptr[orb_size * 4 + i];
    }
  }
}

template<typename ST>
void SplineC2ROMPTarget<ST>::mw_evaluateVGL(const RefVectorWithLeader<SPOSet>& sa_list,
                                            const RefVectorWithLeader<ParticleSet>& P_list,
                                            int iat,
                                            const RefVector<ValueVector>& psi_v_list,
                                            const RefVector<GradVector>& dpsi_v_list,
                                            const RefVector<ValueVector>& d2psi_v_list) const
{
  assert(this == &sa_list.getLeader());
  auto& phi_leader = sa_list.getCastedLeader<SplineC2ROMPTarget<ST>>();
  // make this class unit tests friendly without the need of setup resources.
  if (!phi_leader.mw_mem_)
  {
    app_warning() << "SplineC2ROMPTarget : This message should not be seen in production (performance bug) runs but "
                     "only unit tests (expected)."
                  << std::endl;
    phi_leader.mw_mem_ = std::make_unique<SplineOMPTargetMultiWalkerMem<ST, TT>>();
  }
  auto& mw_mem             = *phi_leader.mw_mem_;
  auto& mw_pos_copy        = mw_mem.mw_pos_copy;
  auto& mw_offload_scratch = mw_mem.mw_offload_scratch;
  auto& mw_results_scratch = mw_mem.mw_results_scratch;
  const int nwalkers       = sa_list.size();
  mw_pos_copy.resize(nwalkers * 6);

  // pack particle positions
  for (int iw = 0; iw < nwalkers; ++iw)
  {
    const PointType& r = P_list[iw].activeR(iat);
    PointType ru(PrimLattice.toUnit_floor(r));
    mw_pos_copy[iw * 6]     = r[0];
    mw_pos_copy[iw * 6 + 1] = r[1];
    mw_pos_copy[iw * 6 + 2] = r[2];
    mw_pos_copy[iw * 6 + 3] = ru[0];
    mw_pos_copy[iw * 6 + 4] = ru[1];
    mw_pos_copy[iw * 6 + 5] = ru[2];
  }

  phi_leader.evaluateVGLMultiPos(mw_pos_copy, mw_offload_scratch, mw_results_scratch, psi_v_list, dpsi_v_list,
                                 d2psi_v_list);
}

template<typename ST>
void SplineC2ROMPTarget<ST>::mw_evaluateVGLandDetRatioGrads(const RefVectorWithLeader<SPOSet>& spo_list,
                                                            const RefVectorWithLeader<ParticleSet>& P_list,
                                                            int iat,
                                                            const std::vector<const ValueType*>& invRow_ptr_list,
                                                            VGLVector& phi_vgl_v,
                                                            std::vector<ValueType>& ratios,
                                                            std::vector<GradType>& grads) const
{
  assert(this == &spo_list.getLeader());
  auto& phi_leader         = spo_list.getCastedLeader<SplineC2ROMPTarget<ST>>();
  auto& mw_mem             = *phi_leader.mw_mem_;
  auto& buffer_H2D         = mw_mem.buffer_H2D;
  auto& rg_private         = mw_mem.rg_private;
  auto& mw_offload_scratch = mw_mem.mw_offload_scratch;
  auto& mw_results_scratch = mw_mem.mw_results_scratch;
  const int nwalkers       = spo_list.size();
  buffer_H2D.resize(nwalkers, sizeof(ST) * 6 + sizeof(ValueType*));

  // pack particle positions and invRow pointers.
  for (int iw = 0; iw < nwalkers; ++iw)
  {
    const PointType& r = P_list[iw].activeR(iat);
    PointType ru(PrimLattice.toUnit_floor(r));
    Vector<ST> pos_copy(reinterpret_cast<ST*>(buffer_H2D[iw]), 6);

    pos_copy[0] = r[0];
    pos_copy[1] = r[1];
    pos_copy[2] = r[2];
    pos_copy[3] = ru[0];
    pos_copy[4] = ru[1];
    pos_copy[5] = ru[2];

    auto& invRow_ptr = *reinterpret_cast<const ValueType**>(buffer_H2D[iw] + sizeof(ST) * 6);
    invRow_ptr       = invRow_ptr_list[iw];
  }

  const size_t num_pos       = nwalkers;
  const int ChunkSizePerTeam = 128;
  const int NumTeams         = (myV.size() + ChunkSizePerTeam - 1) / ChunkSizePerTeam;
  const auto padded_size     = myV.size();
  // for V(1)G(3)H(6) intermediate result
  mw_offload_scratch.resize(padded_size * num_pos * 10);
  const auto orb_size = phi_vgl_v.size() / num_pos;
  // for V(1)G(3)L(1) final result
  mw_results_scratch.resize(orb_size * num_pos * 5);
  // per team ratio and grads
  rg_private.resize(num_pos, NumTeams * 4);

  // Ye: need to extract sizes and pointers before entering target region
  const auto* spline_ptr         = SplineInst->getSplinePtr();
  auto* buffer_H2D_ptr           = buffer_H2D.data();
  auto* offload_scratch_ptr      = mw_offload_scratch.data();
  auto* results_scratch_ptr      = mw_results_scratch.data();
  const auto myKcart_padded_size = myKcart->capacity();
  auto* mKK_ptr                  = mKK->data();
  auto* GGt_ptr                  = GGt_offload->data();
  auto* PrimLattice_G_ptr        = PrimLattice_G_offload->data();
  auto* myKcart_ptr              = myKcart->data();
  auto* phi_vgl_ptr              = phi_vgl_v.data();
  auto* rg_private_ptr           = rg_private.data();
  const size_t buffer_H2D_stride = buffer_H2D.cols();
  const size_t first_spo_local   = first_spo;
  const size_t phi_vgl_stride    = phi_vgl_v.capacity();
  const int nComplexBands_local  = nComplexBands;

  {
    ScopedTimer offload(offload_timer_);
    PRAGMA_OFFLOAD("omp target teams distribute collapse(2) num_teams(NumTeams*num_pos) \
                    map(always, to: buffer_H2D_ptr[:buffer_H2D.size()]) \
                    map(always, from: rg_private_ptr[0:rg_private.size()])")
    for (int iw = 0; iw < num_pos; iw++)
      for (int team_id = 0; team_id < NumTeams; team_id++)
      {
        const int first = ChunkSizePerTeam * team_id;
        const int last  = (first + ChunkSizePerTeam) > padded_size ? padded_size : first + ChunkSizePerTeam;
        auto* restrict offload_scratch_iw_ptr = offload_scratch_ptr + padded_size * iw * 10;
        auto* restrict psi_iw_ptr             = results_scratch_ptr + orb_size * iw * 5;
        const auto* restrict pos_iw_ptr       = reinterpret_cast<ST*>(buffer_H2D_ptr + buffer_H2D_stride * iw);
        const auto* restrict invRow_iw_ptr =
            *reinterpret_cast<ValueType**>(buffer_H2D_ptr + buffer_H2D_stride * iw + sizeof(ST) * 6);

        int ix, iy, iz;
        ST a[4], b[4], c[4], da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
        spline2::computeLocationAndFractional(spline_ptr, pos_iw_ptr[3], pos_iw_ptr[4], pos_iw_ptr[5], ix, iy, iz, a, b,
                                              c, da, db, dc, d2a, d2b, d2c);

        const ST G[9]      = {PrimLattice_G_ptr[0], PrimLattice_G_ptr[1], PrimLattice_G_ptr[2],
                              PrimLattice_G_ptr[3], PrimLattice_G_ptr[4], PrimLattice_G_ptr[5],
                              PrimLattice_G_ptr[6], PrimLattice_G_ptr[7], PrimLattice_G_ptr[8]};
        const ST symGGt[6] = {GGt_ptr[0], GGt_ptr[1] + GGt_ptr[3], GGt_ptr[2] + GGt_ptr[6],
                              GGt_ptr[4], GGt_ptr[5] + GGt_ptr[7], GGt_ptr[8]};

        PRAGMA_OFFLOAD("omp parallel for")
        for (int index = 0; index < last - first; index++)
          spline2offload::evaluate_vgh_impl_v2(spline_ptr, ix, iy, iz, a, b, c, da, db, dc, d2a, d2b, d2c,
                                               offload_scratch_iw_ptr + first,
                                               offload_scratch_iw_ptr + padded_size + first,
                                               offload_scratch_iw_ptr + padded_size * 4 + first, padded_size, first,
                                               index);
        const int first_cplx = first / 2;
        const int last_cplx  = last / 2 < orb_size ? last / 2 : orb_size;
        PRAGMA_OFFLOAD("omp parallel for")
        for (int index = first_cplx; index < last_cplx; index++)
          C2R::assign_vgl(pos_iw_ptr[0], pos_iw_ptr[1], pos_iw_ptr[2], psi_iw_ptr, mKK_ptr, orb_size,
                          offload_scratch_iw_ptr, padded_size, symGGt, G, myKcart_ptr, myKcart_padded_size,
                          first_spo_local, nComplexBands_local, index);

        ValueType* restrict psi   = psi_iw_ptr;
        ValueType* restrict dpsi  = psi_iw_ptr + orb_size;
        ValueType* restrict d2psi = psi_iw_ptr + orb_size * 4;

        ValueType* restrict out_phi_v = phi_vgl_ptr + iw * orb_size;
        ValueType* restrict out_phi_g = phi_vgl_ptr + phi_vgl_stride + iw * orb_size * 3;
        ValueType* restrict out_phi_l = phi_vgl_ptr + phi_vgl_stride * 4 + iw * orb_size;

        const int first_real = first_cplx + std::min(nComplexBands_local, first_cplx);
        const int last_real  = last_cplx + std::min(nComplexBands_local, last_cplx);
        ValueType ratio(0), grad_x(0), grad_y(0), grad_z(0);
        PRAGMA_OFFLOAD("omp parallel for reduction(+: ratio, grad_x, grad_y, grad_z)")
        for (size_t j = first_real; j < last_real; j++)
        {
          out_phi_v[j]         = psi[j];
          out_phi_l[j]         = d2psi[j];
          out_phi_g[j * 3]     = dpsi[j * 3];
          out_phi_g[j * 3 + 1] = dpsi[j * 3 + 1];
          out_phi_g[j * 3 + 2] = dpsi[j * 3 + 2];

          ratio += psi[j] * invRow_iw_ptr[j];
          grad_x += dpsi[j * 3] * invRow_iw_ptr[j];
          grad_y += dpsi[j * 3 + 1] * invRow_iw_ptr[j];
          grad_z += dpsi[j * 3 + 2] * invRow_iw_ptr[j];
        }

        rg_private_ptr[(iw * NumTeams + team_id) * 4]     = ratio;
        rg_private_ptr[(iw * NumTeams + team_id) * 4 + 1] = grad_x;
        rg_private_ptr[(iw * NumTeams + team_id) * 4 + 2] = grad_y;
        rg_private_ptr[(iw * NumTeams + team_id) * 4 + 3] = grad_z;
      }
  }

  for (int iw = 0; iw < num_pos; iw++)
  {
    ValueType ratio(0);
    for (int team_id = 0; team_id < NumTeams; team_id++)
      ratio += rg_private[iw][team_id * 4];
    ratios[iw] = ratio;

    ValueType grad_x(0), grad_y(0), grad_z(0);
    for (int team_id = 0; team_id < NumTeams; team_id++)
    {
      grad_x += rg_private[iw][team_id * 4 + 1];
      grad_y += rg_private[iw][team_id * 4 + 2];
      grad_z += rg_private[iw][team_id * 4 + 3];
    }
    grads[iw] = GradType{grad_x / ratio, grad_y / ratio, grad_z / ratio};
  }
}

template<typename ST>
void SplineC2ROMPTarget<ST>::assign_vgh(const PointType& r,
                                        ValueVector& psi,
                                        GradVector& dpsi,
                                        HessVector& grad_grad_psi,
                                        int first,
                                        int last) const
{
  // protect last
  last = last > kPoints.size() ? kPoints.size() : last;

  const ST g00 = PrimLattice.G(0), g01 = PrimLattice.G(1), g02 = PrimLattice.G(2), g10 = PrimLattice.G(3),
           g11 = PrimLattice.G(4), g12 = PrimLattice.G(5), g20 = PrimLattice.G(6), g21 = PrimLattice.G(7),
           g22 = PrimLattice.G(8);
  const ST x = r[0], y = r[1], z = r[2];

  const ST* restrict k0 = myKcart->data(0);
  const ST* restrict k1 = myKcart->data(1);
  const ST* restrict k2 = myKcart->data(2);

  const ST* restrict g0  = myG.data(0);
  const ST* restrict g1  = myG.data(1);
  const ST* restrict g2  = myG.data(2);
  const ST* restrict h00 = myH.data(0);
  const ST* restrict h01 = myH.data(1);
  const ST* restrict h02 = myH.data(2);
  const ST* restrict h11 = myH.data(3);
  const ST* restrict h12 = myH.data(4);
  const ST* restrict h22 = myH.data(5);

#pragma omp simd
  for (size_t j = first; j < std::min(nComplexBands, last); j++)
  {
    int jr = j << 1;
    int ji = jr + 1;

    const ST kX    = k0[j];
    const ST kY    = k1[j];
    const ST kZ    = k2[j];
    const ST val_r = myV[jr];
    const ST val_i = myV[ji];

    //phase
    ST s, c;
    omptarget::sincos(-(x * kX + y * kY + z * kZ), &s, &c);

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

    const size_t psiIndex = first_spo + jr;

    psi[psiIndex]     = c * val_r - s * val_i;
    dpsi[psiIndex][0] = c * gX_r - s * gX_i;
    dpsi[psiIndex][1] = c * gY_r - s * gY_i;
    dpsi[psiIndex][2] = c * gZ_r - s * gZ_i;

    psi[psiIndex + 1]     = c * val_i + s * val_r;
    dpsi[psiIndex + 1][0] = c * gX_i + s * gX_r;
    dpsi[psiIndex + 1][1] = c * gY_i + s * gY_r;
    dpsi[psiIndex + 1][2] = c * gZ_i + s * gZ_r;

    const ST h_xx_r =
        v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g00, g01, g02) + kX * (gX_i + dX_i);
    const ST h_xy_r =
        v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g10, g11, g12) + kX * (gY_i + dY_i);
    const ST h_xz_r =
        v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g20, g21, g22) + kX * (gZ_i + dZ_i);
    const ST h_yx_r =
        v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g10, g11, g12, g00, g01, g02) + kY * (gX_i + dX_i);
    const ST h_yy_r =
        v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g10, g11, g12, g10, g11, g12) + kY * (gY_i + dY_i);
    const ST h_yz_r =
        v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g10, g11, g12, g20, g21, g22) + kY * (gZ_i + dZ_i);
    const ST h_zx_r =
        v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g20, g21, g22, g00, g01, g02) + kZ * (gX_i + dX_i);
    const ST h_zy_r =
        v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g20, g21, g22, g10, g11, g12) + kZ * (gY_i + dY_i);
    const ST h_zz_r =
        v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g20, g21, g22, g20, g21, g22) + kZ * (gZ_i + dZ_i);

    const ST h_xx_i =
        v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g00, g01, g02) - kX * (gX_r + dX_r);
    const ST h_xy_i =
        v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g10, g11, g12) - kX * (gY_r + dY_r);
    const ST h_xz_i =
        v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g20, g21, g22) - kX * (gZ_r + dZ_r);
    const ST h_yx_i =
        v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g10, g11, g12, g00, g01, g02) - kY * (gX_r + dX_r);
    const ST h_yy_i =
        v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g10, g11, g12, g10, g11, g12) - kY * (gY_r + dY_r);
    const ST h_yz_i =
        v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g10, g11, g12, g20, g21, g22) - kY * (gZ_r + dZ_r);
    const ST h_zx_i =
        v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g20, g21, g22, g00, g01, g02) - kZ * (gX_r + dX_r);
    const ST h_zy_i =
        v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g20, g21, g22, g10, g11, g12) - kZ * (gY_r + dY_r);
    const ST h_zz_i =
        v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g20, g21, g22, g20, g21, g22) - kZ * (gZ_r + dZ_r);

    grad_grad_psi[psiIndex][0] = c * h_xx_r - s * h_xx_i;
    grad_grad_psi[psiIndex][1] = c * h_xy_r - s * h_xy_i;
    grad_grad_psi[psiIndex][2] = c * h_xz_r - s * h_xz_i;
    grad_grad_psi[psiIndex][3] = c * h_yx_r - s * h_yx_i;
    grad_grad_psi[psiIndex][4] = c * h_yy_r - s * h_yy_i;
    grad_grad_psi[psiIndex][5] = c * h_yz_r - s * h_yz_i;
    grad_grad_psi[psiIndex][6] = c * h_zx_r - s * h_zx_i;
    grad_grad_psi[psiIndex][7] = c * h_zy_r - s * h_zy_i;
    grad_grad_psi[psiIndex][8] = c * h_zz_r - s * h_zz_i;

    grad_grad_psi[psiIndex + 1][0] = c * h_xx_i + s * h_xx_r;
    grad_grad_psi[psiIndex + 1][1] = c * h_xy_i + s * h_xy_r;
    grad_grad_psi[psiIndex + 1][2] = c * h_xz_i + s * h_xz_r;
    grad_grad_psi[psiIndex + 1][3] = c * h_yx_i + s * h_yx_r;
    grad_grad_psi[psiIndex + 1][4] = c * h_yy_i + s * h_yy_r;
    grad_grad_psi[psiIndex + 1][5] = c * h_yz_i + s * h_yz_r;
    grad_grad_psi[psiIndex + 1][6] = c * h_zx_i + s * h_zx_r;
    grad_grad_psi[psiIndex + 1][7] = c * h_zy_i + s * h_zy_r;
    grad_grad_psi[psiIndex + 1][8] = c * h_zz_i + s * h_zz_r;
  }

#pragma omp simd
  for (size_t j = std::max(nComplexBands, first); j < last; j++)
  {
    int jr = j << 1;
    int ji = jr + 1;

    const ST kX    = k0[j];
    const ST kY    = k1[j];
    const ST kZ    = k2[j];
    const ST val_r = myV[jr];
    const ST val_i = myV[ji];

    //phase
    ST s, c;
    omptarget::sincos(-(x * kX + y * kY + z * kZ), &s, &c);

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

    const size_t psiIndex = first_spo + nComplexBands + j;

    psi[psiIndex]     = c * val_r - s * val_i;
    dpsi[psiIndex][0] = c * gX_r - s * gX_i;
    dpsi[psiIndex][1] = c * gY_r - s * gY_i;
    dpsi[psiIndex][2] = c * gZ_r - s * gZ_i;

    const ST h_xx_r =
        v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g00, g01, g02) + kX * (gX_i + dX_i);
    const ST h_xy_r =
        v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g10, g11, g12) + kX * (gY_i + dY_i);
    const ST h_xz_r =
        v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g20, g21, g22) + kX * (gZ_i + dZ_i);
    const ST h_yx_r =
        v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g10, g11, g12, g00, g01, g02) + kY * (gX_i + dX_i);
    const ST h_yy_r =
        v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g10, g11, g12, g10, g11, g12) + kY * (gY_i + dY_i);
    const ST h_yz_r =
        v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g10, g11, g12, g20, g21, g22) + kY * (gZ_i + dZ_i);
    const ST h_zx_r =
        v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g20, g21, g22, g00, g01, g02) + kZ * (gX_i + dX_i);
    const ST h_zy_r =
        v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g20, g21, g22, g10, g11, g12) + kZ * (gY_i + dY_i);
    const ST h_zz_r =
        v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g20, g21, g22, g20, g21, g22) + kZ * (gZ_i + dZ_i);

    const ST h_xx_i =
        v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g00, g01, g02) - kX * (gX_r + dX_r);
    const ST h_xy_i =
        v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g10, g11, g12) - kX * (gY_r + dY_r);
    const ST h_xz_i =
        v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g20, g21, g22) - kX * (gZ_r + dZ_r);
    const ST h_yx_i =
        v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g10, g11, g12, g00, g01, g02) - kY * (gX_r + dX_r);
    const ST h_yy_i =
        v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g10, g11, g12, g10, g11, g12) - kY * (gY_r + dY_r);
    const ST h_yz_i =
        v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g10, g11, g12, g20, g21, g22) - kY * (gZ_r + dZ_r);
    const ST h_zx_i =
        v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g20, g21, g22, g00, g01, g02) - kZ * (gX_r + dX_r);
    const ST h_zy_i =
        v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g20, g21, g22, g10, g11, g12) - kZ * (gY_r + dY_r);
    const ST h_zz_i =
        v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g20, g21, g22, g20, g21, g22) - kZ * (gZ_r + dZ_r);

    grad_grad_psi[psiIndex][0] = c * h_xx_r - s * h_xx_i;
    grad_grad_psi[psiIndex][1] = c * h_xy_r - s * h_xy_i;
    grad_grad_psi[psiIndex][2] = c * h_xz_r - s * h_xz_i;
    grad_grad_psi[psiIndex][3] = c * h_yx_r - s * h_yx_i;
    grad_grad_psi[psiIndex][4] = c * h_yy_r - s * h_yy_i;
    grad_grad_psi[psiIndex][5] = c * h_yz_r - s * h_yz_i;
    grad_grad_psi[psiIndex][6] = c * h_zx_r - s * h_zx_i;
    grad_grad_psi[psiIndex][7] = c * h_zy_r - s * h_zy_i;
    grad_grad_psi[psiIndex][8] = c * h_zz_r - s * h_zz_i;
  }
}

template<typename ST>
void SplineC2ROMPTarget<ST>::evaluateVGH(const ParticleSet& P,
                                         const int iat,
                                         ValueVector& psi,
                                         GradVector& dpsi,
                                         HessVector& grad_grad_psi)
{
  const PointType& r = P.activeR(iat);
  PointType ru(PrimLattice.toUnit_floor(r));
#pragma omp parallel
  {
    int first, last;
    FairDivideAligned(myV.size(), getAlignment<ST>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

    spline2::evaluate3d_vgh(SplineInst->getSplinePtr(), ru, myV, myG, myH, first, last);
    assign_vgh(r, psi, dpsi, grad_grad_psi, first / 2, last / 2);
  }
}

template<typename ST>
void SplineC2ROMPTarget<ST>::assign_vghgh(const PointType& r,
                                          ValueVector& psi,
                                          GradVector& dpsi,
                                          HessVector& grad_grad_psi,
                                          GGGVector& grad_grad_grad_psi,
                                          int first,
                                          int last) const
{
  // protect last
  last = last < 0 ? kPoints.size() : (last > kPoints.size() ? kPoints.size() : last);

  const ST g00 = PrimLattice.G(0), g01 = PrimLattice.G(1), g02 = PrimLattice.G(2), g10 = PrimLattice.G(3),
           g11 = PrimLattice.G(4), g12 = PrimLattice.G(5), g20 = PrimLattice.G(6), g21 = PrimLattice.G(7),
           g22 = PrimLattice.G(8);
  const ST x = r[0], y = r[1], z = r[2];

  const ST* restrict k0 = myKcart->data(0);
  const ST* restrict k1 = myKcart->data(1);
  const ST* restrict k2 = myKcart->data(2);

  const ST* restrict g0  = myG.data(0);
  const ST* restrict g1  = myG.data(1);
  const ST* restrict g2  = myG.data(2);
  const ST* restrict h00 = myH.data(0);
  const ST* restrict h01 = myH.data(1);
  const ST* restrict h02 = myH.data(2);
  const ST* restrict h11 = myH.data(3);
  const ST* restrict h12 = myH.data(4);
  const ST* restrict h22 = myH.data(5);

  const ST* restrict gh000 = mygH.data(0);
  const ST* restrict gh001 = mygH.data(1);
  const ST* restrict gh002 = mygH.data(2);
  const ST* restrict gh011 = mygH.data(3);
  const ST* restrict gh012 = mygH.data(4);
  const ST* restrict gh022 = mygH.data(5);
  const ST* restrict gh111 = mygH.data(6);
  const ST* restrict gh112 = mygH.data(7);
  const ST* restrict gh122 = mygH.data(8);
  const ST* restrict gh222 = mygH.data(9);

//SIMD doesn't work quite right yet.  Comment out until further debugging.
#pragma omp simd
  for (size_t j = first; j < std::min(nComplexBands, last); j++)
  {
    int jr = j << 1;
    int ji = jr + 1;

    const ST kX    = k0[j];
    const ST kY    = k1[j];
    const ST kZ    = k2[j];
    const ST val_r = myV[jr];
    const ST val_i = myV[ji];

    //phase
    ST s, c;
    omptarget::sincos(-(x * kX + y * kY + z * kZ), &s, &c);

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

    const size_t psiIndex = first_spo + jr;
    psi[psiIndex]         = c * val_r - s * val_i;
    dpsi[psiIndex][0]     = c * gX_r - s * gX_i;
    dpsi[psiIndex][1]     = c * gY_r - s * gY_i;
    dpsi[psiIndex][2]     = c * gZ_r - s * gZ_i;

    psi[psiIndex + 1]     = c * val_i + s * val_r;
    dpsi[psiIndex + 1][0] = c * gX_i + s * gX_r;
    dpsi[psiIndex + 1][1] = c * gY_i + s * gY_r;
    dpsi[psiIndex + 1][2] = c * gZ_i + s * gZ_r;

    //intermediates for computation of hessian. \partial_i \partial_j phi in cartesian coordinates.
    const ST f_xx_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g00, g01, g02);
    const ST f_xy_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g10, g11, g12);
    const ST f_xz_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g20, g21, g22);
    const ST f_yy_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g10, g11, g12, g10, g11, g12);
    const ST f_yz_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g10, g11, g12, g20, g21, g22);
    const ST f_zz_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g20, g21, g22, g20, g21, g22);

    const ST f_xx_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g00, g01, g02);
    const ST f_xy_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g10, g11, g12);
    const ST f_xz_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g20, g21, g22);
    const ST f_yy_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g10, g11, g12, g10, g11, g12);
    const ST f_yz_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g10, g11, g12, g20, g21, g22);
    const ST f_zz_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g20, g21, g22, g20, g21, g22);

    const ST h_xx_r = f_xx_r + 2 * kX * dX_i - kX * kX * val_r;
    const ST h_xy_r = f_xy_r + (kX * dY_i + kY * dX_i) - kX * kY * val_r;
    const ST h_xz_r = f_xz_r + (kX * dZ_i + kZ * dX_i) - kX * kZ * val_r;
    const ST h_yy_r = f_yy_r + 2 * kY * dY_i - kY * kY * val_r;
    const ST h_yz_r = f_yz_r + (kY * dZ_i + kZ * dY_i) - kY * kZ * val_r;
    const ST h_zz_r = f_zz_r + 2 * kZ * dZ_i - kZ * kZ * val_r;

    const ST h_xx_i = f_xx_i - 2 * kX * dX_r - kX * kX * val_i;
    const ST h_xy_i = f_xy_i - (kX * dY_r + kY * dX_r) - kX * kY * val_i;
    const ST h_xz_i = f_xz_i - (kX * dZ_r + kZ * dX_r) - kX * kZ * val_i;
    const ST h_yy_i = f_yy_i - 2 * kY * dY_r - kY * kY * val_i;
    const ST h_yz_i = f_yz_i - (kZ * dY_r + kY * dZ_r) - kZ * kY * val_i;
    const ST h_zz_i = f_zz_i - 2 * kZ * dZ_r - kZ * kZ * val_i;

    grad_grad_psi[psiIndex][0] = c * h_xx_r - s * h_xx_i;
    grad_grad_psi[psiIndex][1] = c * h_xy_r - s * h_xy_i;
    grad_grad_psi[psiIndex][2] = c * h_xz_r - s * h_xz_i;
    grad_grad_psi[psiIndex][3] = c * h_xy_r - s * h_xy_i;
    grad_grad_psi[psiIndex][4] = c * h_yy_r - s * h_yy_i;
    grad_grad_psi[psiIndex][5] = c * h_yz_r - s * h_yz_i;
    grad_grad_psi[psiIndex][6] = c * h_xz_r - s * h_xz_i;
    grad_grad_psi[psiIndex][7] = c * h_yz_r - s * h_yz_i;
    grad_grad_psi[psiIndex][8] = c * h_zz_r - s * h_zz_i;

    grad_grad_psi[psiIndex + 1][0] = c * h_xx_i + s * h_xx_r;
    grad_grad_psi[psiIndex + 1][1] = c * h_xy_i + s * h_xy_r;
    grad_grad_psi[psiIndex + 1][2] = c * h_xz_i + s * h_xz_r;
    grad_grad_psi[psiIndex + 1][3] = c * h_xy_i + s * h_xy_r;
    grad_grad_psi[psiIndex + 1][4] = c * h_yy_i + s * h_yy_r;
    grad_grad_psi[psiIndex + 1][5] = c * h_yz_i + s * h_yz_r;
    grad_grad_psi[psiIndex + 1][6] = c * h_xz_i + s * h_xz_r;
    grad_grad_psi[psiIndex + 1][7] = c * h_yz_i + s * h_yz_r;
    grad_grad_psi[psiIndex + 1][8] = c * h_zz_i + s * h_zz_r;

    //These are the real and imaginary components of the third SPO derivative.  _xxx denotes
    // third derivative w.r.t. x, _xyz, a derivative with resepect to x,y, and z, and so on.

    const ST f3_xxx_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                    gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g00, g01, g02, g00, g01, g02);
    const ST f3_xxy_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                    gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g00, g01, g02, g10, g11, g12);
    const ST f3_xxz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                    gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g00, g01, g02, g20, g21, g22);
    const ST f3_xyy_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                    gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g10, g11, g12, g10, g11, g12);
    const ST f3_xyz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                    gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g10, g11, g12, g20, g21, g22);
    const ST f3_xzz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                    gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g20, g21, g22, g20, g21, g22);
    const ST f3_yyy_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                    gh112[jr], gh122[jr], gh222[jr], g10, g11, g12, g10, g11, g12, g10, g11, g12);
    const ST f3_yyz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                    gh112[jr], gh122[jr], gh222[jr], g10, g11, g12, g10, g11, g12, g20, g21, g22);
    const ST f3_yzz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                    gh112[jr], gh122[jr], gh222[jr], g10, g11, g12, g20, g21, g22, g20, g21, g22);
    const ST f3_zzz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                    gh112[jr], gh122[jr], gh222[jr], g20, g21, g22, g20, g21, g22, g20, g21, g22);

    const ST f3_xxx_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                    gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g00, g01, g02, g00, g01, g02);
    const ST f3_xxy_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                    gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g00, g01, g02, g10, g11, g12);
    const ST f3_xxz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                    gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g00, g01, g02, g20, g21, g22);
    const ST f3_xyy_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                    gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g10, g11, g12, g10, g11, g12);
    const ST f3_xyz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                    gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g10, g11, g12, g20, g21, g22);
    const ST f3_xzz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                    gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g20, g21, g22, g20, g21, g22);
    const ST f3_yyy_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                    gh112[ji], gh122[ji], gh222[ji], g10, g11, g12, g10, g11, g12, g10, g11, g12);
    const ST f3_yyz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                    gh112[ji], gh122[ji], gh222[ji], g10, g11, g12, g10, g11, g12, g20, g21, g22);
    const ST f3_yzz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                    gh112[ji], gh122[ji], gh222[ji], g10, g11, g12, g20, g21, g22, g20, g21, g22);
    const ST f3_zzz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                    gh112[ji], gh122[ji], gh222[ji], g20, g21, g22, g20, g21, g22, g20, g21, g22);

    //Here is where we build up the components of the physical hessian gradient, namely, d^3/dx^3(e^{-ik*r}\phi(r)
    const ST gh_xxx_r = f3_xxx_r + 3 * kX * f_xx_i - 3 * kX * kX * dX_r - kX * kX * kX * val_i;
    const ST gh_xxx_i = f3_xxx_i - 3 * kX * f_xx_r - 3 * kX * kX * dX_i + kX * kX * kX * val_r;
    const ST gh_xxy_r =
        f3_xxy_r + (kY * f_xx_i + 2 * kX * f_xy_i) - (kX * kX * dY_r + 2 * kX * kY * dX_r) - kX * kX * kY * val_i;
    const ST gh_xxy_i =
        f3_xxy_i - (kY * f_xx_r + 2 * kX * f_xy_r) - (kX * kX * dY_i + 2 * kX * kY * dX_i) + kX * kX * kY * val_r;
    const ST gh_xxz_r =
        f3_xxz_r + (kZ * f_xx_i + 2 * kX * f_xz_i) - (kX * kX * dZ_r + 2 * kX * kZ * dX_r) - kX * kX * kZ * val_i;
    const ST gh_xxz_i =
        f3_xxz_i - (kZ * f_xx_r + 2 * kX * f_xz_r) - (kX * kX * dZ_i + 2 * kX * kZ * dX_i) + kX * kX * kZ * val_r;
    const ST gh_xyy_r =
        f3_xyy_r + (2 * kY * f_xy_i + kX * f_yy_i) - (2 * kX * kY * dY_r + kY * kY * dX_r) - kX * kY * kY * val_i;
    const ST gh_xyy_i =
        f3_xyy_i - (2 * kY * f_xy_r + kX * f_yy_r) - (2 * kX * kY * dY_i + kY * kY * dX_i) + kX * kY * kY * val_r;
    const ST gh_xyz_r = f3_xyz_r + (kX * f_yz_i + kY * f_xz_i + kZ * f_xy_i) -
        (kX * kY * dZ_r + kY * kZ * dX_r + kZ * kX * dY_r) - kX * kY * kZ * val_i;
    const ST gh_xyz_i = f3_xyz_i - (kX * f_yz_r + kY * f_xz_r + kZ * f_xy_r) -
        (kX * kY * dZ_i + kY * kZ * dX_i + kZ * kX * dY_i) + kX * kY * kZ * val_r;
    const ST gh_xzz_r =
        f3_xzz_r + (2 * kZ * f_xz_i + kX * f_zz_i) - (2 * kX * kZ * dZ_r + kZ * kZ * dX_r) - kX * kZ * kZ * val_i;
    const ST gh_xzz_i =
        f3_xzz_i - (2 * kZ * f_xz_r + kX * f_zz_r) - (2 * kX * kZ * dZ_i + kZ * kZ * dX_i) + kX * kZ * kZ * val_r;
    const ST gh_yyy_r = f3_yyy_r + 3 * kY * f_yy_i - 3 * kY * kY * dY_r - kY * kY * kY * val_i;
    const ST gh_yyy_i = f3_yyy_i - 3 * kY * f_yy_r - 3 * kY * kY * dY_i + kY * kY * kY * val_r;
    const ST gh_yyz_r =
        f3_yyz_r + (kZ * f_yy_i + 2 * kY * f_yz_i) - (kY * kY * dZ_r + 2 * kY * kZ * dY_r) - kY * kY * kZ * val_i;
    const ST gh_yyz_i =
        f3_yyz_i - (kZ * f_yy_r + 2 * kY * f_yz_r) - (kY * kY * dZ_i + 2 * kY * kZ * dY_i) + kY * kY * kZ * val_r;
    const ST gh_yzz_r =
        f3_yzz_r + (2 * kZ * f_yz_i + kY * f_zz_i) - (2 * kY * kZ * dZ_r + kZ * kZ * dY_r) - kY * kZ * kZ * val_i;
    const ST gh_yzz_i =
        f3_yzz_i - (2 * kZ * f_yz_r + kY * f_zz_r) - (2 * kY * kZ * dZ_i + kZ * kZ * dY_i) + kY * kZ * kZ * val_r;
    const ST gh_zzz_r = f3_zzz_r + 3 * kZ * f_zz_i - 3 * kZ * kZ * dZ_r - kZ * kZ * kZ * val_i;
    const ST gh_zzz_i = f3_zzz_i - 3 * kZ * f_zz_r - 3 * kZ * kZ * dZ_i + kZ * kZ * kZ * val_r;

    grad_grad_grad_psi[psiIndex][0][0] = c * gh_xxx_r - s * gh_xxx_i;
    grad_grad_grad_psi[psiIndex][0][1] = c * gh_xxy_r - s * gh_xxy_i;
    grad_grad_grad_psi[psiIndex][0][2] = c * gh_xxz_r - s * gh_xxz_i;
    grad_grad_grad_psi[psiIndex][0][3] = c * gh_xxy_r - s * gh_xxy_i;
    grad_grad_grad_psi[psiIndex][0][4] = c * gh_xyy_r - s * gh_xyy_i;
    grad_grad_grad_psi[psiIndex][0][5] = c * gh_xyz_r - s * gh_xyz_i;
    grad_grad_grad_psi[psiIndex][0][6] = c * gh_xxz_r - s * gh_xxz_i;
    grad_grad_grad_psi[psiIndex][0][7] = c * gh_xyz_r - s * gh_xyz_i;
    grad_grad_grad_psi[psiIndex][0][8] = c * gh_xzz_r - s * gh_xzz_i;

    grad_grad_grad_psi[psiIndex][1][0] = c * gh_xxy_r - s * gh_xxy_i;
    grad_grad_grad_psi[psiIndex][1][1] = c * gh_xyy_r - s * gh_xyy_i;
    grad_grad_grad_psi[psiIndex][1][2] = c * gh_xyz_r - s * gh_xyz_i;
    grad_grad_grad_psi[psiIndex][1][3] = c * gh_xyy_r - s * gh_xyy_i;
    grad_grad_grad_psi[psiIndex][1][4] = c * gh_yyy_r - s * gh_yyy_i;
    grad_grad_grad_psi[psiIndex][1][5] = c * gh_yyz_r - s * gh_yyz_i;
    grad_grad_grad_psi[psiIndex][1][6] = c * gh_xyz_r - s * gh_xyz_i;
    grad_grad_grad_psi[psiIndex][1][7] = c * gh_yyz_r - s * gh_yyz_i;
    grad_grad_grad_psi[psiIndex][1][8] = c * gh_yzz_r - s * gh_yzz_i;

    grad_grad_grad_psi[psiIndex][2][0] = c * gh_xxz_r - s * gh_xxz_i;
    grad_grad_grad_psi[psiIndex][2][1] = c * gh_xyz_r - s * gh_xyz_i;
    grad_grad_grad_psi[psiIndex][2][2] = c * gh_xzz_r - s * gh_xzz_i;
    grad_grad_grad_psi[psiIndex][2][3] = c * gh_xyz_r - s * gh_xyz_i;
    grad_grad_grad_psi[psiIndex][2][4] = c * gh_yyz_r - s * gh_yyz_i;
    grad_grad_grad_psi[psiIndex][2][5] = c * gh_yzz_r - s * gh_yzz_i;
    grad_grad_grad_psi[psiIndex][2][6] = c * gh_xzz_r - s * gh_xzz_i;
    grad_grad_grad_psi[psiIndex][2][7] = c * gh_yzz_r - s * gh_yzz_i;
    grad_grad_grad_psi[psiIndex][2][8] = c * gh_zzz_r - s * gh_zzz_i;

    grad_grad_grad_psi[psiIndex + 1][0][0] = c * gh_xxx_i + s * gh_xxx_r;
    grad_grad_grad_psi[psiIndex + 1][0][1] = c * gh_xxy_i + s * gh_xxy_r;
    grad_grad_grad_psi[psiIndex + 1][0][2] = c * gh_xxz_i + s * gh_xxz_r;
    grad_grad_grad_psi[psiIndex + 1][0][3] = c * gh_xxy_i + s * gh_xxy_r;
    grad_grad_grad_psi[psiIndex + 1][0][4] = c * gh_xyy_i + s * gh_xyy_r;
    grad_grad_grad_psi[psiIndex + 1][0][5] = c * gh_xyz_i + s * gh_xyz_r;
    grad_grad_grad_psi[psiIndex + 1][0][6] = c * gh_xxz_i + s * gh_xxz_r;
    grad_grad_grad_psi[psiIndex + 1][0][7] = c * gh_xyz_i + s * gh_xyz_r;
    grad_grad_grad_psi[psiIndex + 1][0][8] = c * gh_xzz_i + s * gh_xzz_r;

    grad_grad_grad_psi[psiIndex + 1][1][0] = c * gh_xxy_i + s * gh_xxy_r;
    grad_grad_grad_psi[psiIndex + 1][1][1] = c * gh_xyy_i + s * gh_xyy_r;
    grad_grad_grad_psi[psiIndex + 1][1][2] = c * gh_xyz_i + s * gh_xyz_r;
    grad_grad_grad_psi[psiIndex + 1][1][3] = c * gh_xyy_i + s * gh_xyy_r;
    grad_grad_grad_psi[psiIndex + 1][1][4] = c * gh_yyy_i + s * gh_yyy_r;
    grad_grad_grad_psi[psiIndex + 1][1][5] = c * gh_yyz_i + s * gh_yyz_r;
    grad_grad_grad_psi[psiIndex + 1][1][6] = c * gh_xyz_i + s * gh_xyz_r;
    grad_grad_grad_psi[psiIndex + 1][1][7] = c * gh_yyz_i + s * gh_yyz_r;
    grad_grad_grad_psi[psiIndex + 1][1][8] = c * gh_yzz_i + s * gh_yzz_r;

    grad_grad_grad_psi[psiIndex + 1][2][0] = c * gh_xxz_i + s * gh_xxz_r;
    grad_grad_grad_psi[psiIndex + 1][2][1] = c * gh_xyz_i + s * gh_xyz_r;
    grad_grad_grad_psi[psiIndex + 1][2][2] = c * gh_xzz_i + s * gh_xzz_r;
    grad_grad_grad_psi[psiIndex + 1][2][3] = c * gh_xyz_i + s * gh_xyz_r;
    grad_grad_grad_psi[psiIndex + 1][2][4] = c * gh_yyz_i + s * gh_yyz_r;
    grad_grad_grad_psi[psiIndex + 1][2][5] = c * gh_yzz_i + s * gh_yzz_r;
    grad_grad_grad_psi[psiIndex + 1][2][6] = c * gh_xzz_i + s * gh_xzz_r;
    grad_grad_grad_psi[psiIndex + 1][2][7] = c * gh_yzz_i + s * gh_yzz_r;
    grad_grad_grad_psi[psiIndex + 1][2][8] = c * gh_zzz_i + s * gh_zzz_r;
  }
#pragma omp simd
  for (size_t j = std::max(nComplexBands, first); j < last; j++)
  {
    int jr = j << 1;
    int ji = jr + 1;

    const ST kX    = k0[j];
    const ST kY    = k1[j];
    const ST kZ    = k2[j];
    const ST val_r = myV[jr];
    const ST val_i = myV[ji];

    //phase
    ST s, c;
    omptarget::sincos(-(x * kX + y * kY + z * kZ), &s, &c);

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

    const size_t psiIndex = first_spo + nComplexBands + j;
    psi[psiIndex]         = c * val_r - s * val_i;
    dpsi[psiIndex][0]     = c * gX_r - s * gX_i;
    dpsi[psiIndex][1]     = c * gY_r - s * gY_i;
    dpsi[psiIndex][2]     = c * gZ_r - s * gZ_i;

    //intermediates for computation of hessian. \partial_i \partial_j phi in cartesian coordinates.
    const ST f_xx_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g00, g01, g02);
    const ST f_xy_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g10, g11, g12);
    const ST f_xz_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g00, g01, g02, g20, g21, g22);
    const ST f_yy_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g10, g11, g12, g10, g11, g12);
    const ST f_yz_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g10, g11, g12, g20, g21, g22);
    const ST f_zz_r = v_m_v(h00[jr], h01[jr], h02[jr], h11[jr], h12[jr], h22[jr], g20, g21, g22, g20, g21, g22);

    const ST f_xx_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g00, g01, g02);
    const ST f_xy_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g10, g11, g12);
    const ST f_xz_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g00, g01, g02, g20, g21, g22);
    const ST f_yy_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g10, g11, g12, g10, g11, g12);
    const ST f_yz_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g10, g11, g12, g20, g21, g22);
    const ST f_zz_i = v_m_v(h00[ji], h01[ji], h02[ji], h11[ji], h12[ji], h22[ji], g20, g21, g22, g20, g21, g22);

    const ST h_xx_r = f_xx_r + 2 * kX * dX_i - kX * kX * val_r;
    const ST h_xy_r = f_xy_r + (kX * dY_i + kY * dX_i) - kX * kY * val_r;
    const ST h_xz_r = f_xz_r + (kX * dZ_i + kZ * dX_i) - kX * kZ * val_r;
    const ST h_yy_r = f_yy_r + 2 * kY * dY_i - kY * kY * val_r;
    const ST h_yz_r = f_yz_r + (kY * dZ_i + kZ * dY_i) - kY * kZ * val_r;
    const ST h_zz_r = f_zz_r + 2 * kZ * dZ_i - kZ * kZ * val_r;

    const ST h_xx_i = f_xx_i - 2 * kX * dX_r - kX * kX * val_i;
    const ST h_xy_i = f_xy_i - (kX * dY_r + kY * dX_r) - kX * kY * val_i;
    const ST h_xz_i = f_xz_i - (kX * dZ_r + kZ * dX_r) - kX * kZ * val_i;
    const ST h_yy_i = f_yy_i - 2 * kY * dY_r - kY * kY * val_i;
    const ST h_yz_i = f_yz_i - (kZ * dY_r + kY * dZ_r) - kZ * kY * val_i;
    const ST h_zz_i = f_zz_i - 2 * kZ * dZ_r - kZ * kZ * val_i;

    grad_grad_psi[psiIndex][0] = c * h_xx_r - s * h_xx_i;
    grad_grad_psi[psiIndex][1] = c * h_xy_r - s * h_xy_i;
    grad_grad_psi[psiIndex][2] = c * h_xz_r - s * h_xz_i;
    grad_grad_psi[psiIndex][3] = c * h_xy_r - s * h_xy_i;
    grad_grad_psi[psiIndex][4] = c * h_yy_r - s * h_yy_i;
    grad_grad_psi[psiIndex][5] = c * h_yz_r - s * h_yz_i;
    grad_grad_psi[psiIndex][6] = c * h_xz_r - s * h_xz_i;
    grad_grad_psi[psiIndex][7] = c * h_yz_r - s * h_yz_i;
    grad_grad_psi[psiIndex][8] = c * h_zz_r - s * h_zz_i;

    //These are the real and imaginary components of the third SPO derivative.  _xxx denotes
    // third derivative w.r.t. x, _xyz, a derivative with resepect to x,y, and z, and so on.

    const ST f3_xxx_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                    gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g00, g01, g02, g00, g01, g02);
    const ST f3_xxy_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                    gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g00, g01, g02, g10, g11, g12);
    const ST f3_xxz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                    gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g00, g01, g02, g20, g21, g22);
    const ST f3_xyy_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                    gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g10, g11, g12, g10, g11, g12);
    const ST f3_xyz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                    gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g10, g11, g12, g20, g21, g22);
    const ST f3_xzz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                    gh112[jr], gh122[jr], gh222[jr], g00, g01, g02, g20, g21, g22, g20, g21, g22);
    const ST f3_yyy_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                    gh112[jr], gh122[jr], gh222[jr], g10, g11, g12, g10, g11, g12, g10, g11, g12);
    const ST f3_yyz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                    gh112[jr], gh122[jr], gh222[jr], g10, g11, g12, g10, g11, g12, g20, g21, g22);
    const ST f3_yzz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                    gh112[jr], gh122[jr], gh222[jr], g10, g11, g12, g20, g21, g22, g20, g21, g22);
    const ST f3_zzz_r = t3_contract(gh000[jr], gh001[jr], gh002[jr], gh011[jr], gh012[jr], gh022[jr], gh111[jr],
                                    gh112[jr], gh122[jr], gh222[jr], g20, g21, g22, g20, g21, g22, g20, g21, g22);

    const ST f3_xxx_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                    gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g00, g01, g02, g00, g01, g02);
    const ST f3_xxy_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                    gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g00, g01, g02, g10, g11, g12);
    const ST f3_xxz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                    gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g00, g01, g02, g20, g21, g22);
    const ST f3_xyy_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                    gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g10, g11, g12, g10, g11, g12);
    const ST f3_xyz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                    gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g10, g11, g12, g20, g21, g22);
    const ST f3_xzz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                    gh112[ji], gh122[ji], gh222[ji], g00, g01, g02, g20, g21, g22, g20, g21, g22);
    const ST f3_yyy_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                    gh112[ji], gh122[ji], gh222[ji], g10, g11, g12, g10, g11, g12, g10, g11, g12);
    const ST f3_yyz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                    gh112[ji], gh122[ji], gh222[ji], g10, g11, g12, g10, g11, g12, g20, g21, g22);
    const ST f3_yzz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                    gh112[ji], gh122[ji], gh222[ji], g10, g11, g12, g20, g21, g22, g20, g21, g22);
    const ST f3_zzz_i = t3_contract(gh000[ji], gh001[ji], gh002[ji], gh011[ji], gh012[ji], gh022[ji], gh111[ji],
                                    gh112[ji], gh122[ji], gh222[ji], g20, g21, g22, g20, g21, g22, g20, g21, g22);

    //Here is where we build up the components of the physical hessian gradient, namely, d^3/dx^3(e^{-ik*r}\phi(r)
    const ST gh_xxx_r = f3_xxx_r + 3 * kX * f_xx_i - 3 * kX * kX * dX_r - kX * kX * kX * val_i;
    const ST gh_xxx_i = f3_xxx_i - 3 * kX * f_xx_r - 3 * kX * kX * dX_i + kX * kX * kX * val_r;
    const ST gh_xxy_r =
        f3_xxy_r + (kY * f_xx_i + 2 * kX * f_xy_i) - (kX * kX * dY_r + 2 * kX * kY * dX_r) - kX * kX * kY * val_i;
    const ST gh_xxy_i =
        f3_xxy_i - (kY * f_xx_r + 2 * kX * f_xy_r) - (kX * kX * dY_i + 2 * kX * kY * dX_i) + kX * kX * kY * val_r;
    const ST gh_xxz_r =
        f3_xxz_r + (kZ * f_xx_i + 2 * kX * f_xz_i) - (kX * kX * dZ_r + 2 * kX * kZ * dX_r) - kX * kX * kZ * val_i;
    const ST gh_xxz_i =
        f3_xxz_i - (kZ * f_xx_r + 2 * kX * f_xz_r) - (kX * kX * dZ_i + 2 * kX * kZ * dX_i) + kX * kX * kZ * val_r;
    const ST gh_xyy_r =
        f3_xyy_r + (2 * kY * f_xy_i + kX * f_yy_i) - (2 * kX * kY * dY_r + kY * kY * dX_r) - kX * kY * kY * val_i;
    const ST gh_xyy_i =
        f3_xyy_i - (2 * kY * f_xy_r + kX * f_yy_r) - (2 * kX * kY * dY_i + kY * kY * dX_i) + kX * kY * kY * val_r;
    const ST gh_xyz_r = f3_xyz_r + (kX * f_yz_i + kY * f_xz_i + kZ * f_xy_i) -
        (kX * kY * dZ_r + kY * kZ * dX_r + kZ * kX * dY_r) - kX * kY * kZ * val_i;
    const ST gh_xyz_i = f3_xyz_i - (kX * f_yz_r + kY * f_xz_r + kZ * f_xy_r) -
        (kX * kY * dZ_i + kY * kZ * dX_i + kZ * kX * dY_i) + kX * kY * kZ * val_r;
    const ST gh_xzz_r =
        f3_xzz_r + (2 * kZ * f_xz_i + kX * f_zz_i) - (2 * kX * kZ * dZ_r + kZ * kZ * dX_r) - kX * kZ * kZ * val_i;
    const ST gh_xzz_i =
        f3_xzz_i - (2 * kZ * f_xz_r + kX * f_zz_r) - (2 * kX * kZ * dZ_i + kZ * kZ * dX_i) + kX * kZ * kZ * val_r;
    const ST gh_yyy_r = f3_yyy_r + 3 * kY * f_yy_i - 3 * kY * kY * dY_r - kY * kY * kY * val_i;
    const ST gh_yyy_i = f3_yyy_i - 3 * kY * f_yy_r - 3 * kY * kY * dY_i + kY * kY * kY * val_r;
    const ST gh_yyz_r =
        f3_yyz_r + (kZ * f_yy_i + 2 * kY * f_yz_i) - (kY * kY * dZ_r + 2 * kY * kZ * dY_r) - kY * kY * kZ * val_i;
    const ST gh_yyz_i =
        f3_yyz_i - (kZ * f_yy_r + 2 * kY * f_yz_r) - (kY * kY * dZ_i + 2 * kY * kZ * dY_i) + kY * kY * kZ * val_r;
    const ST gh_yzz_r =
        f3_yzz_r + (2 * kZ * f_yz_i + kY * f_zz_i) - (2 * kY * kZ * dZ_r + kZ * kZ * dY_r) - kY * kZ * kZ * val_i;
    const ST gh_yzz_i =
        f3_yzz_i - (2 * kZ * f_yz_r + kY * f_zz_r) - (2 * kY * kZ * dZ_i + kZ * kZ * dY_i) + kY * kZ * kZ * val_r;
    const ST gh_zzz_r = f3_zzz_r + 3 * kZ * f_zz_i - 3 * kZ * kZ * dZ_r - kZ * kZ * kZ * val_i;
    const ST gh_zzz_i = f3_zzz_i - 3 * kZ * f_zz_r - 3 * kZ * kZ * dZ_i + kZ * kZ * kZ * val_r;
    //[x][xx] //These are the unique entries
    grad_grad_grad_psi[psiIndex][0][0] = c * gh_xxx_r - s * gh_xxx_i;
    grad_grad_grad_psi[psiIndex][0][1] = c * gh_xxy_r - s * gh_xxy_i;
    grad_grad_grad_psi[psiIndex][0][2] = c * gh_xxz_r - s * gh_xxz_i;
    grad_grad_grad_psi[psiIndex][0][3] = c * gh_xxy_r - s * gh_xxy_i;
    grad_grad_grad_psi[psiIndex][0][4] = c * gh_xyy_r - s * gh_xyy_i;
    grad_grad_grad_psi[psiIndex][0][5] = c * gh_xyz_r - s * gh_xyz_i;
    grad_grad_grad_psi[psiIndex][0][6] = c * gh_xxz_r - s * gh_xxz_i;
    grad_grad_grad_psi[psiIndex][0][7] = c * gh_xyz_r - s * gh_xyz_i;
    grad_grad_grad_psi[psiIndex][0][8] = c * gh_xzz_r - s * gh_xzz_i;

    grad_grad_grad_psi[psiIndex][1][0] = c * gh_xxy_r - s * gh_xxy_i;
    grad_grad_grad_psi[psiIndex][1][1] = c * gh_xyy_r - s * gh_xyy_i;
    grad_grad_grad_psi[psiIndex][1][2] = c * gh_xyz_r - s * gh_xyz_i;
    grad_grad_grad_psi[psiIndex][1][3] = c * gh_xyy_r - s * gh_xyy_i;
    grad_grad_grad_psi[psiIndex][1][4] = c * gh_yyy_r - s * gh_yyy_i;
    grad_grad_grad_psi[psiIndex][1][5] = c * gh_yyz_r - s * gh_yyz_i;
    grad_grad_grad_psi[psiIndex][1][6] = c * gh_xyz_r - s * gh_xyz_i;
    grad_grad_grad_psi[psiIndex][1][7] = c * gh_yyz_r - s * gh_yyz_i;
    grad_grad_grad_psi[psiIndex][1][8] = c * gh_yzz_r - s * gh_yzz_i;

    grad_grad_grad_psi[psiIndex][2][0] = c * gh_xxz_r - s * gh_xxz_i;
    grad_grad_grad_psi[psiIndex][2][1] = c * gh_xyz_r - s * gh_xyz_i;
    grad_grad_grad_psi[psiIndex][2][2] = c * gh_xzz_r - s * gh_xzz_i;
    grad_grad_grad_psi[psiIndex][2][3] = c * gh_xyz_r - s * gh_xyz_i;
    grad_grad_grad_psi[psiIndex][2][4] = c * gh_yyz_r - s * gh_yyz_i;
    grad_grad_grad_psi[psiIndex][2][5] = c * gh_yzz_r - s * gh_yzz_i;
    grad_grad_grad_psi[psiIndex][2][6] = c * gh_xzz_r - s * gh_xzz_i;
    grad_grad_grad_psi[psiIndex][2][7] = c * gh_yzz_r - s * gh_yzz_i;
    grad_grad_grad_psi[psiIndex][2][8] = c * gh_zzz_r - s * gh_zzz_i;
  }
}

template<typename ST>
void SplineC2ROMPTarget<ST>::evaluateVGHGH(const ParticleSet& P,
                                           const int iat,
                                           ValueVector& psi,
                                           GradVector& dpsi,
                                           HessVector& grad_grad_psi,
                                           GGGVector& grad_grad_grad_psi)
{
  const PointType& r = P.activeR(iat);
  PointType ru(PrimLattice.toUnit_floor(r));
#pragma omp parallel
  {
    int first, last;
    FairDivideAligned(myV.size(), getAlignment<ST>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

    spline2::evaluate3d_vghgh(SplineInst->getSplinePtr(), ru, myV, myG, myH, mygH, first, last);
    assign_vghgh(r, psi, dpsi, grad_grad_psi, grad_grad_grad_psi, first / 2, last / 2);
  }
}

template<typename ST>
void SplineC2ROMPTarget<ST>::evaluate_notranspose(const ParticleSet& P,
                                                  int first,
                                                  int last,
                                                  ValueMatrix& logdet,
                                                  GradMatrix& dlogdet,
                                                  ValueMatrix& d2logdet)
{
  // chunk the [first, last) loop into blocks to save temporary memory usage
  const int block_size = 16;

  // reference vectors refer to the rows of matrices
  std::vector<ValueVector> multi_psi_v;
  std::vector<GradVector> multi_dpsi_v;
  std::vector<ValueVector> multi_d2psi_v;
  RefVector<ValueVector> psi_v_list;
  RefVector<GradVector> dpsi_v_list;
  RefVector<ValueVector> d2psi_v_list;

  multi_psi_v.reserve(block_size);
  multi_dpsi_v.reserve(block_size);
  multi_d2psi_v.reserve(block_size);
  psi_v_list.reserve(block_size);
  dpsi_v_list.reserve(block_size);
  d2psi_v_list.reserve(block_size);

  for (int iat = first, i = 0; iat < last; iat += block_size, i += block_size)
  {
    const int actual_block_size = std::min(last - iat, block_size);
    multi_pos_copy.resize(actual_block_size * 6);
    multi_psi_v.clear();
    multi_dpsi_v.clear();
    multi_d2psi_v.clear();
    psi_v_list.clear();
    dpsi_v_list.clear();
    d2psi_v_list.clear();

    for (int ipos = 0; ipos < actual_block_size; ++ipos)
    {
      // pack particle positions
      const PointType& r = P.activeR(iat + ipos);
      PointType ru(PrimLattice.toUnit_floor(r));
      multi_pos_copy[ipos * 6]     = r[0];
      multi_pos_copy[ipos * 6 + 1] = r[1];
      multi_pos_copy[ipos * 6 + 2] = r[2];
      multi_pos_copy[ipos * 6 + 3] = ru[0];
      multi_pos_copy[ipos * 6 + 4] = ru[1];
      multi_pos_copy[ipos * 6 + 5] = ru[2];

      multi_psi_v.emplace_back(logdet[i + ipos], OrbitalSetSize);
      multi_dpsi_v.emplace_back(dlogdet[i + ipos], OrbitalSetSize);
      multi_d2psi_v.emplace_back(d2logdet[i + ipos], OrbitalSetSize);

      psi_v_list.push_back(multi_psi_v[ipos]);
      dpsi_v_list.push_back(multi_dpsi_v[ipos]);
      d2psi_v_list.push_back(multi_d2psi_v[ipos]);
    }

    evaluateVGLMultiPos(multi_pos_copy, offload_scratch, results_scratch, psi_v_list, dpsi_v_list, d2psi_v_list);
  }
}

template class SplineC2ROMPTarget<float>;
template class SplineC2ROMPTarget<double>;

} // namespace qmcplusplus
