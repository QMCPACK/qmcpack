//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "SplineR2R.h"
#include "Concurrency/OpenMP.h"
#include "spline2/MultiBspline.hpp"
#include "spline2/MultiBsplineOffload.hpp"
#include "spline2/MultiBsplineEval.hpp"
#include "spline2/MultiBsplineEval_OMPoffload.hpp"
#include "QMCWaveFunctions/BsplineFactory/contraction_helper.hpp"
#include "Platforms/CPU/BLAS.hpp"
#include "CPU/SIMD/inner_product.hpp"
#include "OMPTarget/OMPTargetMath.hpp"
#include <cstdint>

namespace qmcplusplus
{

template<typename T>
inline static std::unique_ptr<MultiBsplineBase<T>> create_MultiBsplineDerived(const bool use_offload)
{
  if (use_offload)
    return std::make_unique<MultiBsplineOffload<T>>();
  else
    return std::make_unique<MultiBspline<T>>();
}

template<typename ST>
SplineR2R<ST>::SplineR2R(const std::string& my_name, bool use_offload)
    : BsplineSet(my_name),
      use_offload_(use_offload),
      offload_timer_(createGlobalTimer("SplineC2ROMPTarget::offload", timer_level_fine)),
      SplineInst(create_MultiBsplineDerived<ST>(use_offload))
{}

template<typename ST>
SplineR2R<ST>::SplineR2R(const SplineR2R& in) = default;

template<typename ST>
inline void SplineR2R<ST>::set_spline(SingleSplineType* spline_r,
                                      SingleSplineType* spline_i,
                                      int twist,
                                      int ispline,
                                      int level)
{
  copy_spline<double, ST>(*spline_r, *SplineInst->getSplinePtr(), ispline);
}

template<typename ST>
void SplineR2R<ST>::finalizeConstruction()
{
  if (use_offload_)
  {
    SplineInst->finalize();
    // transfer static data to GPU
    GGt_offload           = std::make_shared<OffloadVector<ST>>(9);
    PrimLattice_G_offload = std::make_shared<OffloadVector<ST>>(9);
    for (std::uint32_t i = 0; i < 9; i++)
    {
      (*GGt_offload)[i]           = GGt[i];
      (*PrimLattice_G_offload)[i] = PrimLattice.G[i];
    }
    PrimLattice_G_offload->updateTo();
    GGt_offload->updateTo();
  }
}

template<typename ST>
bool SplineR2R<ST>::read_splines(hdf_archive& h5f)
{
  std::ostringstream o;
  o << "spline_" << MyIndex;
  einspline_engine<SplineType> bigtable(SplineInst->getSplinePtr());
  return h5f.readEntry(bigtable, o.str().c_str()); //"spline_0");
}

template<typename ST>
bool SplineR2R<ST>::write_splines(hdf_archive& h5f)
{
  std::ostringstream o;
  o << "spline_" << MyIndex;
  einspline_engine<SplineType> bigtable(SplineInst->getSplinePtr());
  return h5f.writeEntry(bigtable, o.str().c_str()); //"spline_0");
}

template<typename ST>
void SplineR2R<ST>::storeParamsBeforeRotation()
{
  const auto spline_ptr     = SplineInst->getSplinePtr();
  const auto coefs_tot_size = spline_ptr->coefs_size;
  coef_copy_                = std::make_shared<std::vector<ST>>(coefs_tot_size);

  std::copy_n(spline_ptr->coefs, coefs_tot_size, coef_copy_->begin());
}

/*
  ~~ Notes for rotation ~~
  spl_coefs      = Raw pointer to spline coefficients
  basis_set_size = Number of spline coefs per orbital
  OrbitalSetSize = Number of orbitals (excluding padding)

  spl_coefs has a complicated layout depending on dimensionality of splines.
  Luckily, for our purposes, we can think of spl_coefs as pointing to a
  matrix of size BasisSetSize x (OrbitalSetSize + padding), with the spline
  index adjacent in memory. The orbital index is SIMD aligned and therefore
  may include padding.

  As a result, due to SIMD alignment, Nsplines may be larger than the
  actual number of splined orbitals. This means that in practice rot_mat
  may be smaller than the number of 'columns' in the coefs array!

      SplineR2R spl_coef layout:
             ^         | sp1 | ... | spN | pad |
             |         |=====|=====|=====|=====|
             |         | c11 | ... | c1N | 0   |
      basis_set_size   | c21 | ... | c2N | 0   |
             |         | ... | ... | ... | 0   |
             |         | cM1 | ... | cMN | 0   |
             v         |=====|=====|=====|=====|
                       <------ Nsplines ------>

      SplineC2C spl_coef layout:
             ^         | sp1_r | sp1_i |  ...  | spN_r | spN_i |  pad  |
             |         |=======|=======|=======|=======|=======|=======|
             |         | c11_r | c11_i |  ...  | c1N_r | c1N_i |   0   |
      basis_set_size   | c21_r | c21_i |  ...  | c2N_r | c2N_i |   0   |
             |         |  ...  |  ...  |  ...  |  ...  |  ...  |  ...  |
             |         | cM1_r | cM1_i |  ...  | cMN_r | cMN_i |   0   |
             v         |=======|=======|=======|=======|=======|=======|
                       <------------------ Nsplines ------------------>

  NB: For splines (typically) BasisSetSize >> OrbitalSetSize, so the spl_coefs
  "matrix" is very tall and skinny.
*/
template<typename ST>
void SplineR2R<ST>::applyRotation(const ValueMatrix& rot_mat, bool use_stored_copy)
{
  // SplineInst is a MultiBspline. See src/spline2/MultiBspline.hpp
  const auto spline_ptr = SplineInst->getSplinePtr();
  assert(spline_ptr != nullptr);
  const auto spl_coefs      = spline_ptr->coefs;
  const auto Nsplines       = spline_ptr->num_splines; // May include padding
  const auto coefs_tot_size = spline_ptr->coefs_size;
  const auto BasisSetSize   = coefs_tot_size / Nsplines;
  const auto TrueNOrbs      = rot_mat.size1(); // == Nsplines - padding
  assert(OrbitalSetSize == rot_mat.rows());
  assert(OrbitalSetSize == rot_mat.cols());

  if (!use_stored_copy)
  {
    assert(coef_copy_ != nullptr);
    std::copy_n(spl_coefs, coefs_tot_size, coef_copy_->begin());
  }


  if constexpr (std::is_same_v<ST, RealType>)
  {
    //Here, ST should be equal to ValueType, which will be double for R2R. Using BLAS to make things faster
    BLAS::gemm('N', 'N', OrbitalSetSize, BasisSetSize, OrbitalSetSize, ST(1.0), rot_mat.data(), OrbitalSetSize,
               coef_copy_->data(), Nsplines, ST(0.0), spl_coefs, Nsplines);
  }
  else
  {
    //Here, ST is float but ValueType is double for R2R. Due to issues with type conversions, just doing naive matrix multiplication in this case to not lose precision on rot_mat
    for (IndexType i = 0; i < BasisSetSize; i++)
      for (IndexType j = 0; j < OrbitalSetSize; j++)
      {
        const auto cur_elem = Nsplines * i + j;
        RealType newval{0.};
        for (IndexType k = 0; k < OrbitalSetSize; k++)
        {
          const auto index = i * Nsplines + k;
          newval += (*coef_copy_)[index] * rot_mat[k][j];
        }
        spl_coefs[cur_elem] = newval;
      }
  }
  // update coefficients on GPU from host
  SplineInst->finalize();
}


template<typename ST>
inline void SplineR2R<ST>::assign_v(int bc_sign, const vContainer_type& myV, ValueVector& psi, int first, int last)
    const
{
  // protect last against kPoints.size() and psi.size()
  size_t last_real = std::min(kPoints.size(), psi.size());
  last             = last > last_real ? last_real : last;

  const ST signed_one = (bc_sign & 1) ? -1 : 1;
#pragma omp simd
  for (size_t j = first; j < last; ++j)
    psi[first_spo + j] = signed_one * myV[j];
}

template<typename ST>
void SplineR2R<ST>::evaluateValue(const ParticleSet& P, const int iat, ValueVector& psi)
{
  const PointType& r = P.activeR(iat);
  PointType ru;
  int bc_sign = convertPos(r, ru);

#pragma omp parallel
  {
    int first, last;
    FairDivideAligned(psi.size(), getAlignment<ST>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

    spline2::evaluate3d(SplineInst->getSplinePtr(), ru, myV, first, last);
    assign_v(bc_sign, myV, psi, first, last);
  }
}

template<typename ST>
void SplineR2R<ST>::evaluateDetRatios(const VirtualParticleSet& VP,
                                      ValueVector& psi,
                                      const ValueVector& psiinv,
                                      std::vector<ValueType>& ratios)
{
  const bool need_resize = ratios_private.rows() < VP.getTotalNum();

#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    // initialize thread private ratios
    if (need_resize)
    {
      if (tid == 0) // just like #pragma omp master, but one fewer call to the runtime
        ratios_private.resize(VP.getTotalNum(), omp_get_num_threads());
#pragma omp barrier
    }
    int first, last;
    FairDivideAligned(psi.size(), getAlignment<ST>(), omp_get_num_threads(), tid, first, last);
    const int last_real = kPoints.size() < last ? kPoints.size() : last;

    for (int iat = 0; iat < VP.getTotalNum(); ++iat)
    {
      const PointType& r = VP.activeR(iat);
      PointType ru;
      int bc_sign = convertPos(r, ru);

      spline2::evaluate3d(SplineInst->getSplinePtr(), ru, myV, first, last);
      assign_v(bc_sign, myV, psi, first, last_real);
      ratios_private[iat][tid] = simd::dot(psi.data() + first, psiinv.data() + first, last_real - first);
    }
  }

  // do the reduction manually
  for (int iat = 0; iat < VP.getTotalNum(); ++iat)
  {
    ratios[iat] = ValueType(0);
    for (int tid = 0; tid < ratios_private.cols(); tid++)
      ratios[iat] += ratios_private[iat][tid];
  }
}

template<typename ST>
void SplineR2R<ST>::mw_evaluateDetRatios(const RefVectorWithLeader<SPOSet>& spo_list,
                                         const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                                         const RefVector<ValueVector>& psi_list,
                                         const std::vector<const ValueType*>& invRow_ptr_list,
                                         std::vector<std::vector<ValueType>>& ratios_list) const
{
  if (!use_offload_)
  {
    BsplineSet::mw_evaluateDetRatios(spo_list, vp_list, psi_list, invRow_ptr_list, ratios_list);
    return;
  }

  assert(this == &spo_list.getLeader());
  auto& phi_leader                = spo_list.getCastedLeader<SplineR2R>();
  auto& mw_mem                    = phi_leader.mw_mem_handle_.getResource();
  auto& det_ratios_buffer_H2D     = mw_mem.det_ratios_buffer_H2D;
  auto& mw_ratios_private         = mw_mem.mw_ratios_private;
  auto& mw_offload_scratch        = mw_mem.mw_offload_scratch;
  auto& mw_results_scratch        = mw_mem.mw_results_scratch;
  const size_t nw                 = spo_list.size();
  const size_t requested_orb_size = psi_list.size() ? psi_list[0].get().size() : phi_leader.size();

  size_t mw_nVP = 0;
  for (const VirtualParticleSet& VP : vp_list)
    mw_nVP += VP.getTotalNum();

  const size_t packed_size = nw * sizeof(ValueType*) + mw_nVP * (4 * sizeof(ST) + sizeof(int));
  det_ratios_buffer_H2D.resize(packed_size);

  // pack invRow_ptr_list to det_ratios_buffer_H2D
  Vector<const ValueType*> ptr_buffer(reinterpret_cast<const ValueType**>(det_ratios_buffer_H2D.data()), nw);
  for (size_t iw = 0; iw < nw; iw++)
    ptr_buffer[iw] = invRow_ptr_list[iw];

  // pack particle positions
  auto* pos_ptr = reinterpret_cast<ST*>(det_ratios_buffer_H2D.data() + nw * sizeof(ValueType*));
  auto* ref_id_ptr =
      reinterpret_cast<int*>(det_ratios_buffer_H2D.data() + nw * sizeof(ValueType*) + mw_nVP * 4 * sizeof(ST));
  size_t iVP = 0;
  for (size_t iw = 0; iw < nw; iw++)
  {
    const VirtualParticleSet& VP = vp_list[iw];
    assert(ratios_list[iw].size() == VP.getTotalNum());
    for (size_t iat = 0; iat < VP.getTotalNum(); ++iat, ++iVP)
    {
      ref_id_ptr[iVP]    = iw;
      const PointType& r = VP.activeR(iat);
      PointType ru;
      const int bc_sign = phi_leader.convertPos(r, ru);

      pos_ptr[0] = (bc_sign & 1) ? ST(-1) : ST(1);
      pos_ptr[1] = ru[0];
      pos_ptr[2] = ru[1];
      pos_ptr[3] = ru[2];
      pos_ptr += 4;
    }
  }

  const size_t ChunkSizePerTeam = 512;
  const int NumTeams            = (myV.size() + ChunkSizePerTeam - 1) / ChunkSizePerTeam;
  mw_ratios_private.resize(mw_nVP, NumTeams);
  const auto spline_padded_size = myV.size();
  mw_offload_scratch.resize(spline_padded_size * mw_nVP);

  // Ye: need to extract sizes and pointers before entering target region
  const auto* spline_ptr       = SplineInst->getSplinePtr();
  auto* offload_scratch_ptr    = mw_offload_scratch.data();
  auto* buffer_H2D_ptr         = det_ratios_buffer_H2D.data();
  auto* ratios_private_ptr     = mw_ratios_private.data();
  const size_t first_spo_local = first_spo;

  {
    ScopedTimer offload(offload_timer_);
    PRAGMA_OFFLOAD("omp target teams distribute collapse(2) num_teams(NumTeams*mw_nVP) \
                map(always, to: buffer_H2D_ptr[0:det_ratios_buffer_H2D.size()]) \
                map(always, from: ratios_private_ptr[0:NumTeams*mw_nVP])")
    for (int iat = 0; iat < mw_nVP; iat++)
      for (int team_id = 0; team_id < NumTeams; team_id++)
      {
        const size_t first = ChunkSizePerTeam * team_id;
        const size_t last  = omptarget::min(first + ChunkSizePerTeam, spline_padded_size);

        auto* restrict offload_scratch_iat_ptr = offload_scratch_ptr + spline_padded_size * iat;
        auto* ref_id_ptr = reinterpret_cast<int*>(buffer_H2D_ptr + nw * sizeof(ValueType*) + mw_nVP * 4 * sizeof(ST));
        auto* restrict psiinv_ptr  = reinterpret_cast<const ValueType**>(buffer_H2D_ptr)[ref_id_ptr[iat]];
        auto* restrict pos_scratch = reinterpret_cast<ST*>(buffer_H2D_ptr + nw * sizeof(ValueType*));

        int ix, iy, iz;
        ST a[4], b[4], c[4];
        spline2::computeLocationAndFractional(spline_ptr, pos_scratch[iat * 4 + 1], pos_scratch[iat * 4 + 2],
                                              pos_scratch[iat * 4 + 3], ix, iy, iz, a, b, c);

        PRAGMA_OFFLOAD("omp parallel for")
        for (int index = 0; index < last - first; index++)
          spline2offload::evaluate_v_impl_v2(spline_ptr, ix, iy, iz, first + index, a, b, c,
                                             offload_scratch_iat_ptr + first + index);
        ValueType sum(0);
        PRAGMA_OFFLOAD("omp parallel for simd reduction(+:sum)")
        for (int index = first; index < requested_orb_size; index++)
        {
          sum += psiinv_ptr[index] * offload_scratch_iat_ptr[spline_padded_size * SoAFields3D::VAL + index] *
              pos_scratch[iat * 4];
        }
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
      ratios[iat] = ValueType(0);
      for (int tid = 0; tid < NumTeams; ++tid)
        ratios[iat] += mw_ratios_private[iVP][tid];
    }
  }
}

template<typename ST>
inline void SplineR2R<ST>::assign_vgl(int bc_sign,
                                      ValueVector& psi,
                                      GradVector& dpsi,
                                      ValueVector& d2psi,
                                      int first,
                                      int last) const
{
  // protect last against kPoints.size() and psi.size()
  size_t last_real = std::min(kPoints.size(), psi.size());
  last             = last > last_real ? last_real : last;

  const ST signed_one = (bc_sign & 1) ? -1 : 1;
  const ST g00 = PrimLattice.G(0), g01 = PrimLattice.G(1), g02 = PrimLattice.G(2), g10 = PrimLattice.G(3),
           g11 = PrimLattice.G(4), g12 = PrimLattice.G(5), g20 = PrimLattice.G(6), g21 = PrimLattice.G(7),
           g22      = PrimLattice.G(8);
  const ST symGG[6] = {GGt[0], GGt[1] + GGt[3], GGt[2] + GGt[6], GGt[4], GGt[5] + GGt[7], GGt[8]};

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
  for (size_t j = first; j < last; ++j)
  {
    const size_t psiIndex = first_spo + j;
    psi[psiIndex]         = signed_one * myV[j];
    dpsi[psiIndex][0]     = signed_one * (g00 * g0[j] + g01 * g1[j] + g02 * g2[j]);
    dpsi[psiIndex][1]     = signed_one * (g10 * g0[j] + g11 * g1[j] + g12 * g2[j]);
    dpsi[psiIndex][2]     = signed_one * (g20 * g0[j] + g21 * g1[j] + g22 * g2[j]);
    d2psi[psiIndex]       = signed_one * SymTrace(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], symGG);
  }
}

/** assign_vgl_from_l can be used when myL is precomputed and myV,myG,myL in cartesian
   */
template<typename ST>
inline void SplineR2R<ST>::assign_vgl_from_l(int bc_sign, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi)
{
  const ST signed_one   = (bc_sign & 1) ? -1 : 1;
  const ST* restrict g0 = myG.data(0);
  const ST* restrict g1 = myG.data(1);
  const ST* restrict g2 = myG.data(2);

  const size_t last_real = last_spo > psi.size() ? psi.size() : last_spo;
#pragma omp simd
  for (int psiIndex = first_spo; psiIndex < last_real; ++psiIndex)
  {
    const size_t j    = psiIndex - first_spo;
    psi[psiIndex]     = signed_one * myV[j];
    dpsi[psiIndex][0] = signed_one * g0[j];
    dpsi[psiIndex][1] = signed_one * g1[j];
    dpsi[psiIndex][2] = signed_one * g2[j];
    d2psi[psiIndex]   = signed_one * myL[j];
  }
}

template<typename ST>
void SplineR2R<ST>::evaluateVGL(const ParticleSet& P,
                                const int iat,
                                ValueVector& psi,
                                GradVector& dpsi,
                                ValueVector& d2psi)
{
  const PointType& r = P.activeR(iat);
  PointType ru;
  int bc_sign = convertPos(r, ru);

#pragma omp parallel
  {
    int first, last;
    FairDivideAligned(psi.size(), getAlignment<ST>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

    spline2::evaluate3d_vgh(SplineInst->getSplinePtr(), ru, myV, myG, myH, first, last);
    assign_vgl(bc_sign, psi, dpsi, d2psi, first, last);
  }
}

template<typename ST>
void SplineR2R<ST>::mw_evaluateVGLandDetRatioGrads(const RefVectorWithLeader<SPOSet>& spo_list,
                                                   const RefVectorWithLeader<ParticleSet>& P_list,
                                                   int iat,
                                                   const std::vector<const ValueType*>& invRow_ptr_list,
                                                   OffloadMWVGLArray& phi_vgl_v,
                                                   std::vector<ValueType>& ratios,
                                                   std::vector<GradType>& grads) const
{
  if (!use_offload_)
  {
    BsplineSet::mw_evaluateVGLandDetRatioGrads(spo_list, P_list, iat, invRow_ptr_list, phi_vgl_v, ratios, grads);
    return;
  }

  assert(this == &spo_list.getLeader());
  auto& phi_leader         = spo_list.getCastedLeader<SplineR2R>();
  auto& mw_mem             = phi_leader.mw_mem_handle_.getResource();
  auto& buffer_H2D         = mw_mem.buffer_H2D;
  auto& rg_private         = mw_mem.rg_private;
  auto& mw_offload_scratch = mw_mem.mw_offload_scratch;
  const int nwalkers       = spo_list.size();
  buffer_H2D.resize(nwalkers, sizeof(ST) * 4 + sizeof(ValueType*));

  // pack particle positions and invRow pointers.
  for (int iw = 0; iw < nwalkers; ++iw)
  {
    const PointType& r = P_list[iw].activeR(iat);
    PointType ru;
    const int bc_sign = phi_leader.convertPos(r, ru);
    Vector<ST> pos_copy(reinterpret_cast<ST*>(buffer_H2D[iw]), 4);

    pos_copy[0] = (bc_sign & 1) ? ST(-1) : ST(1);
    pos_copy[1] = ru[0];
    pos_copy[2] = ru[1];
    pos_copy[3] = ru[2];

    auto& invRow_ptr = *reinterpret_cast<const ValueType**>(buffer_H2D[iw] + sizeof(ST) * 4);
    invRow_ptr       = invRow_ptr_list[iw];
  }

  const size_t num_pos          = nwalkers;
  const auto spline_padded_size = myV.size();
  const size_t ChunkSizePerTeam = 512;
  const int NumTeams            = (myV.size() + ChunkSizePerTeam - 1) / ChunkSizePerTeam;

  mw_offload_scratch.resize(spline_padded_size * num_pos * SoAFields3D::NUM_FIELDS);
  // per team ratio and grads
  rg_private.resize(num_pos, NumTeams * 4);

  // Ye: need to extract sizes and pointers before entering target region
  const auto* spline_ptr         = SplineInst->getSplinePtr();
  auto* buffer_H2D_ptr           = buffer_H2D.data();
  auto* offload_scratch_ptr      = mw_offload_scratch.data();
  auto* GGt_ptr                  = GGt_offload->data();
  auto* PrimLattice_G_ptr        = PrimLattice_G_offload->data();
  auto* phi_vgl_ptr              = phi_vgl_v.data();
  auto* rg_private_ptr           = rg_private.data();
  const size_t buffer_H2D_stride = buffer_H2D.cols();
  const size_t first_spo_local   = first_spo;
  const auto requested_orb_size  = phi_vgl_v.size(2);
  const size_t phi_vgl_stride    = num_pos * requested_orb_size;

  {
    ScopedTimer offload(offload_timer_);
    PRAGMA_OFFLOAD("omp target teams distribute collapse(2) num_teams(NumTeams*num_pos) \
                    map(always, to: buffer_H2D_ptr[:buffer_H2D.size()])")
    for (int iw = 0; iw < num_pos; iw++)
      for (int team_id = 0; team_id < NumTeams; team_id++)
      {
        const size_t first = ChunkSizePerTeam * team_id;
        const size_t last  = omptarget::min(first + ChunkSizePerTeam, spline_padded_size);

        auto* restrict offload_scratch_iw_ptr = offload_scratch_ptr + spline_padded_size * iw * SoAFields3D::NUM_FIELDS;
        const auto* restrict pos_iw_ptr       = reinterpret_cast<ST*>(buffer_H2D_ptr + buffer_H2D_stride * iw);
        int ix, iy, iz;
        ST a[4], b[4], c[4], da[4], db[4], dc[4], d2a[4], d2b[4], d2c[4];
        spline2::computeLocationAndFractional(spline_ptr, pos_iw_ptr[1], pos_iw_ptr[2], pos_iw_ptr[3], ix, iy, iz, a, b,
                                              c, da, db, dc, d2a, d2b, d2c);
        const ST symGGt[6] = {GGt_ptr[0], GGt_ptr[1] + GGt_ptr[3], GGt_ptr[2] + GGt_ptr[6],
                              GGt_ptr[4], GGt_ptr[5] + GGt_ptr[7], GGt_ptr[8]};

        PRAGMA_OFFLOAD("omp parallel for")
        for (int index = 0; index < last - first; index++)
        {
          spline2offload::evaluate_vgh_impl_v2(spline_ptr, ix, iy, iz, first + index, a, b, c, da, db, dc, d2a, d2b,
                                               d2c, offload_scratch_iw_ptr + first + index, spline_padded_size);
          const int output_index = first + index;
          offload_scratch_iw_ptr[spline_padded_size * SoAFields3D::LAPL + output_index] =
              SymTrace(offload_scratch_iw_ptr[spline_padded_size * SoAFields3D::HESS00 + output_index],
                       offload_scratch_iw_ptr[spline_padded_size * SoAFields3D::HESS01 + output_index],
                       offload_scratch_iw_ptr[spline_padded_size * SoAFields3D::HESS02 + output_index],
                       offload_scratch_iw_ptr[spline_padded_size * SoAFields3D::HESS11 + output_index],
                       offload_scratch_iw_ptr[spline_padded_size * SoAFields3D::HESS12 + output_index],
                       offload_scratch_iw_ptr[spline_padded_size * SoAFields3D::HESS22 + output_index], symGGt);
        }
      }
  }

  {
    ScopedTimer offload(offload_timer_);
    PRAGMA_OFFLOAD("omp target teams distribute collapse(2) num_teams(NumTeams*num_pos) \
                    map(always, to: buffer_H2D_ptr[:buffer_H2D.size()]) \
                    map(always, from: rg_private_ptr[:NumTeams*num_pos*4])")
    for (int iw = 0; iw < num_pos; iw++)
      for (int team_id = 0; team_id < NumTeams; team_id++)
      {
        const size_t first       = ChunkSizePerTeam * team_id;
        const size_t last        = omptarget::min(first + ChunkSizePerTeam, spline_padded_size);
        const size_t reduce_last = omptarget::min(last, static_cast<size_t>(requested_orb_size));

        auto* restrict offload_scratch_iw_ptr = offload_scratch_ptr + spline_padded_size * iw * SoAFields3D::NUM_FIELDS;
        const auto* restrict pos_iw_ptr       = reinterpret_cast<ST*>(buffer_H2D_ptr + buffer_H2D_stride * iw);
        const auto* restrict invRow_iw_ptr =
            *reinterpret_cast<ValueType**>(buffer_H2D_ptr + buffer_H2D_stride * iw + sizeof(ST) * 4);
        const ST G[9] = {PrimLattice_G_ptr[0], PrimLattice_G_ptr[1], PrimLattice_G_ptr[2],
                         PrimLattice_G_ptr[3], PrimLattice_G_ptr[4], PrimLattice_G_ptr[5],
                         PrimLattice_G_ptr[6], PrimLattice_G_ptr[7], PrimLattice_G_ptr[8]};

        ValueType ratio(0), grad_x(0), grad_y(0), grad_z(0);
        PRAGMA_OFFLOAD("omp parallel for reduction(+: ratio, grad_x, grad_y, grad_z)")
        for (int index = first; index < requested_orb_size; index++)
        {
          const ST* restrict val   = offload_scratch_iw_ptr + spline_padded_size * SoAFields3D::VAL;
          const ST* restrict g0    = offload_scratch_iw_ptr + spline_padded_size * SoAFields3D::GRAD0;
          const ST* restrict g1    = offload_scratch_iw_ptr + spline_padded_size * SoAFields3D::GRAD1;
          const ST* restrict g2    = offload_scratch_iw_ptr + spline_padded_size * SoAFields3D::GRAD2;
          const ST* restrict lcart = offload_scratch_iw_ptr + spline_padded_size * SoAFields3D::LAPL;

          ValueType* restrict out_phi    = phi_vgl_ptr + iw * requested_orb_size;
          ValueType* restrict out_dphi_x = out_phi + phi_vgl_stride;
          ValueType* restrict out_dphi_y = out_dphi_x + phi_vgl_stride;
          ValueType* restrict out_dphi_z = out_dphi_y + phi_vgl_stride;
          ValueType* restrict out_d2phi  = out_dphi_z + phi_vgl_stride;

          const size_t psiIndex = first_spo_local + index;
          out_phi[psiIndex]     = pos_iw_ptr[0] * val[index];
          out_dphi_x[psiIndex]  = pos_iw_ptr[0] * (G[0] * g0[index] + G[1] * g1[index] + G[2] * g2[index]);
          out_dphi_y[psiIndex]  = pos_iw_ptr[0] * (G[3] * g0[index] + G[4] * g1[index] + G[5] * g2[index]);
          out_dphi_z[psiIndex]  = pos_iw_ptr[0] * (G[6] * g0[index] + G[7] * g1[index] + G[8] * g2[index]);
          out_d2phi[psiIndex]   = pos_iw_ptr[0] * lcart[index];

          ratio += out_phi[psiIndex] * invRow_iw_ptr[psiIndex];
          grad_x += out_dphi_x[psiIndex] * invRow_iw_ptr[psiIndex];
          grad_y += out_dphi_y[psiIndex] * invRow_iw_ptr[psiIndex];
          grad_z += out_dphi_z[psiIndex] * invRow_iw_ptr[psiIndex];
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
void SplineR2R<ST>::assign_vgh(int bc_sign,
                               ValueVector& psi,
                               GradVector& dpsi,
                               HessVector& grad_grad_psi,
                               int first,
                               int last) const
{
  // protect last against kPoints.size() and psi.size()
  const size_t last_real = std::min(kPoints.size(), psi.size());
  last                   = last > last_real ? last_real : last;

  const ST signed_one = (bc_sign & 1) ? -1 : 1;
  const ST g00 = PrimLattice.G(0), g01 = PrimLattice.G(1), g02 = PrimLattice.G(2), g10 = PrimLattice.G(3),
           g11 = PrimLattice.G(4), g12 = PrimLattice.G(5), g20 = PrimLattice.G(6), g21 = PrimLattice.G(7),
           g22 = PrimLattice.G(8);

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
  for (size_t j = first; j < last; ++j)
  {
    //dot(PrimLattice.G,myG[j])
    const ST dX_r = g00 * g0[j] + g01 * g1[j] + g02 * g2[j];
    const ST dY_r = g10 * g0[j] + g11 * g1[j] + g12 * g2[j];
    const ST dZ_r = g20 * g0[j] + g21 * g1[j] + g22 * g2[j];

    const size_t psiIndex = j + first_spo;
    psi[psiIndex]         = signed_one * myV[j];
    dpsi[psiIndex][0]     = signed_one * dX_r;
    dpsi[psiIndex][1]     = signed_one * dY_r;
    dpsi[psiIndex][2]     = signed_one * dZ_r;

    const ST h_xx_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g00, g01, g02, g00, g01, g02);
    const ST h_xy_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g00, g01, g02, g10, g11, g12);
    const ST h_xz_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g00, g01, g02, g20, g21, g22);
    const ST h_yx_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g10, g11, g12, g00, g01, g02);
    const ST h_yy_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g10, g11, g12, g10, g11, g12);
    const ST h_yz_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g10, g11, g12, g20, g21, g22);
    const ST h_zx_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g20, g21, g22, g00, g01, g02);
    const ST h_zy_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g20, g21, g22, g10, g11, g12);
    const ST h_zz_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g20, g21, g22, g20, g21, g22);

    grad_grad_psi[psiIndex][0] = signed_one * h_xx_r;
    grad_grad_psi[psiIndex][1] = signed_one * h_xy_r;
    grad_grad_psi[psiIndex][2] = signed_one * h_xz_r;
    grad_grad_psi[psiIndex][3] = signed_one * h_yx_r;
    grad_grad_psi[psiIndex][4] = signed_one * h_yy_r;
    grad_grad_psi[psiIndex][5] = signed_one * h_yz_r;
    grad_grad_psi[psiIndex][6] = signed_one * h_zx_r;
    grad_grad_psi[psiIndex][7] = signed_one * h_zy_r;
    grad_grad_psi[psiIndex][8] = signed_one * h_zz_r;
  }
}

template<typename ST>
void SplineR2R<ST>::evaluateVGH(const ParticleSet& P,
                                const int iat,
                                ValueVector& psi,
                                GradVector& dpsi,
                                HessVector& grad_grad_psi)
{
  const PointType& r = P.activeR(iat);
  PointType ru;
  int bc_sign = convertPos(r, ru);

#pragma omp parallel
  {
    int first, last;
    FairDivideAligned(psi.size(), getAlignment<ST>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

    spline2::evaluate3d_vgh(SplineInst->getSplinePtr(), ru, myV, myG, myH, first, last);
    assign_vgh(bc_sign, psi, dpsi, grad_grad_psi, first, last);
  }
}

template<typename ST>
void SplineR2R<ST>::assign_vghgh(int bc_sign,
                                 ValueVector& psi,
                                 GradVector& dpsi,
                                 HessVector& grad_grad_psi,
                                 GGGVector& grad_grad_grad_psi,
                                 int first,
                                 int last) const
{
  // protect last against kPoints.size() and psi.size()
  const size_t last_real = std::min(kPoints.size(), psi.size());
  last                   = last < 0 ? last_real : (last > last_real ? last_real : last);

  const ST signed_one = (bc_sign & 1) ? -1 : 1;
  const ST g00 = PrimLattice.G(0), g01 = PrimLattice.G(1), g02 = PrimLattice.G(2), g10 = PrimLattice.G(3),
           g11 = PrimLattice.G(4), g12 = PrimLattice.G(5), g20 = PrimLattice.G(6), g21 = PrimLattice.G(7),
           g22 = PrimLattice.G(8);

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
  //#pragma omp simd
  for (size_t j = first; j < last; ++j)
  {
    const ST val_r = myV[j];


    //dot(PrimLattice.G,myG[j])
    const ST dX_r = g00 * g0[j] + g01 * g1[j] + g02 * g2[j];
    const ST dY_r = g10 * g0[j] + g11 * g1[j] + g12 * g2[j];
    const ST dZ_r = g20 * g0[j] + g21 * g1[j] + g22 * g2[j];

    const size_t psiIndex = j + first_spo;
    psi[psiIndex]         = signed_one * val_r;
    dpsi[psiIndex][0]     = signed_one * dX_r;
    dpsi[psiIndex][1]     = signed_one * dY_r;
    dpsi[psiIndex][2]     = signed_one * dZ_r;

    //intermediates for computation of hessian. \partial_i \partial_j phi in cartesian coordinates.
    const ST f_xx_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g00, g01, g02, g00, g01, g02);
    const ST f_xy_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g00, g01, g02, g10, g11, g12);
    const ST f_xz_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g00, g01, g02, g20, g21, g22);
    const ST f_yy_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g10, g11, g12, g10, g11, g12);
    const ST f_yz_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g10, g11, g12, g20, g21, g22);
    const ST f_zz_r = v_m_v(h00[j], h01[j], h02[j], h11[j], h12[j], h22[j], g20, g21, g22, g20, g21, g22);

    /*    const ST h_xx_r=f_xx_r;
      const ST h_xy_r=f_xy_r+(kX*dY_i+kY*dX_i)-kX*kY*val_r;
      const ST h_xz_r=f_xz_r+(kX*dZ_i+kZ*dX_i)-kX*kZ*val_r;
      const ST h_yy_r=f_yy_r+2*kY*dY_i-kY*kY*val_r;
      const ST h_yz_r=f_yz_r+(kY*dZ_i+kZ*dY_i)-kY*kZ*val_r;
      const ST h_zz_r=f_zz_r+2*kZ*dZ_i-kZ*kZ*val_r; */

    grad_grad_psi[psiIndex][0] = f_xx_r * signed_one;
    grad_grad_psi[psiIndex][1] = f_xy_r * signed_one;
    grad_grad_psi[psiIndex][2] = f_xz_r * signed_one;
    grad_grad_psi[psiIndex][4] = f_yy_r * signed_one;
    grad_grad_psi[psiIndex][5] = f_yz_r * signed_one;
    grad_grad_psi[psiIndex][8] = f_zz_r * signed_one;

    //symmetry:
    grad_grad_psi[psiIndex][3] = grad_grad_psi[psiIndex][1];
    grad_grad_psi[psiIndex][6] = grad_grad_psi[psiIndex][2];
    grad_grad_psi[psiIndex][7] = grad_grad_psi[psiIndex][5];
    //These are the real and imaginary components of the third SPO derivative.  _xxx denotes
    // third derivative w.r.t. x, _xyz, a derivative with resepect to x,y, and z, and so on.

    const ST f3_xxx_r = t3_contract(gh000[j], gh001[j], gh002[j], gh011[j], gh012[j], gh022[j], gh111[j], gh112[j],
                                    gh122[j], gh222[j], g00, g01, g02, g00, g01, g02, g00, g01, g02);
    const ST f3_xxy_r = t3_contract(gh000[j], gh001[j], gh002[j], gh011[j], gh012[j], gh022[j], gh111[j], gh112[j],
                                    gh122[j], gh222[j], g00, g01, g02, g00, g01, g02, g10, g11, g12);
    const ST f3_xxz_r = t3_contract(gh000[j], gh001[j], gh002[j], gh011[j], gh012[j], gh022[j], gh111[j], gh112[j],
                                    gh122[j], gh222[j], g00, g01, g02, g00, g01, g02, g20, g21, g22);
    const ST f3_xyy_r = t3_contract(gh000[j], gh001[j], gh002[j], gh011[j], gh012[j], gh022[j], gh111[j], gh112[j],
                                    gh122[j], gh222[j], g00, g01, g02, g10, g11, g12, g10, g11, g12);
    const ST f3_xyz_r = t3_contract(gh000[j], gh001[j], gh002[j], gh011[j], gh012[j], gh022[j], gh111[j], gh112[j],
                                    gh122[j], gh222[j], g00, g01, g02, g10, g11, g12, g20, g21, g22);
    const ST f3_xzz_r = t3_contract(gh000[j], gh001[j], gh002[j], gh011[j], gh012[j], gh022[j], gh111[j], gh112[j],
                                    gh122[j], gh222[j], g00, g01, g02, g20, g21, g22, g20, g21, g22);
    const ST f3_yyy_r = t3_contract(gh000[j], gh001[j], gh002[j], gh011[j], gh012[j], gh022[j], gh111[j], gh112[j],
                                    gh122[j], gh222[j], g10, g11, g12, g10, g11, g12, g10, g11, g12);
    const ST f3_yyz_r = t3_contract(gh000[j], gh001[j], gh002[j], gh011[j], gh012[j], gh022[j], gh111[j], gh112[j],
                                    gh122[j], gh222[j], g10, g11, g12, g10, g11, g12, g20, g21, g22);
    const ST f3_yzz_r = t3_contract(gh000[j], gh001[j], gh002[j], gh011[j], gh012[j], gh022[j], gh111[j], gh112[j],
                                    gh122[j], gh222[j], g10, g11, g12, g20, g21, g22, g20, g21, g22);
    const ST f3_zzz_r = t3_contract(gh000[j], gh001[j], gh002[j], gh011[j], gh012[j], gh022[j], gh111[j], gh112[j],
                                    gh122[j], gh222[j], g20, g21, g22, g20, g21, g22, g20, g21, g22);

    //Here is where we build up the components of the physical hessian gradient, namely, d^3/dx^3(e^{-ik*r}\phi(r)
    /*     const ST gh_xxx_r= f3_xxx_r + 3*kX*f_xx_i - 3*kX*kX*dX_r - kX*kX*kX*val_i;
      const ST gh_xxy_r= f3_xxy_r +(kY*f_xx_i+2*kX*f_xy_i) - (kX*kX*dY_r+2*kX*kY*dX_r)-kX*kX*kY*val_i;
      const ST gh_xxz_r= f3_xxz_r +(kZ*f_xx_i+2*kX*f_xz_i) - (kX*kX*dZ_r+2*kX*kZ*dX_r)-kX*kX*kZ*val_i;
      const ST gh_xyy_r= f3_xyy_r +(2*kY*f_xy_i+kX*f_yy_i) - (2*kX*kY*dY_r+kY*kY*dX_r)-kX*kY*kY*val_i;
      const ST gh_xyz_r= f3_xyz_r +(kX*f_yz_i+kY*f_xz_i+kZ*f_xy_i)-(kX*kY*dZ_r+kY*kZ*dX_r+kZ*kX*dY_r) - kX*kY*kZ*val_i;
      const ST gh_xzz_r= f3_xzz_r +(2*kZ*f_xz_i+kX*f_zz_i) - (2*kX*kZ*dZ_r+kZ*kZ*dX_r)-kX*kZ*kZ*val_i;
      const ST gh_yyy_r= f3_yyy_r + 3*kY*f_yy_i - 3*kY*kY*dY_r - kY*kY*kY*val_i;
      const ST gh_yyz_r= f3_yyz_r +(kZ*f_yy_i+2*kY*f_yz_i) - (kY*kY*dZ_r+2*kY*kZ*dY_r)-kY*kY*kZ*val_i;
      const ST gh_yzz_r= f3_yzz_r +(2*kZ*f_yz_i+kY*f_zz_i) - (2*kY*kZ*dZ_r+kZ*kZ*dY_r)-kY*kZ*kZ*val_i;
      const ST gh_zzz_r= f3_zzz_r + 3*kZ*f_zz_i - 3*kZ*kZ*dZ_r - kZ*kZ*kZ*val_i;*/
    //[x][xx] //These are the unique entries
    grad_grad_grad_psi[psiIndex][0][0] = signed_one * f3_xxx_r;
    grad_grad_grad_psi[psiIndex][0][1] = signed_one * f3_xxy_r;
    grad_grad_grad_psi[psiIndex][0][2] = signed_one * f3_xxz_r;
    grad_grad_grad_psi[psiIndex][0][4] = signed_one * f3_xyy_r;
    grad_grad_grad_psi[psiIndex][0][5] = signed_one * f3_xyz_r;
    grad_grad_grad_psi[psiIndex][0][8] = signed_one * f3_xzz_r;

    //filling in the symmetric terms.  Filling out the xij terms
    grad_grad_grad_psi[psiIndex][0][3] = grad_grad_grad_psi[psiIndex][0][1];
    grad_grad_grad_psi[psiIndex][0][6] = grad_grad_grad_psi[psiIndex][0][2];
    grad_grad_grad_psi[psiIndex][0][7] = grad_grad_grad_psi[psiIndex][0][5];

    //Now for everything that's a permutation of the above:
    grad_grad_grad_psi[psiIndex][1][0] = grad_grad_grad_psi[psiIndex][0][1];
    grad_grad_grad_psi[psiIndex][1][1] = grad_grad_grad_psi[psiIndex][0][4];
    grad_grad_grad_psi[psiIndex][1][2] = grad_grad_grad_psi[psiIndex][0][5];
    grad_grad_grad_psi[psiIndex][1][3] = grad_grad_grad_psi[psiIndex][0][4];
    grad_grad_grad_psi[psiIndex][1][6] = grad_grad_grad_psi[psiIndex][0][5];

    grad_grad_grad_psi[psiIndex][2][0] = grad_grad_grad_psi[psiIndex][0][2];
    grad_grad_grad_psi[psiIndex][2][1] = grad_grad_grad_psi[psiIndex][0][5];
    grad_grad_grad_psi[psiIndex][2][2] = grad_grad_grad_psi[psiIndex][0][8];
    grad_grad_grad_psi[psiIndex][2][3] = grad_grad_grad_psi[psiIndex][0][5];
    grad_grad_grad_psi[psiIndex][2][6] = grad_grad_grad_psi[psiIndex][0][8];

    grad_grad_grad_psi[psiIndex][1][4] = signed_one * f3_yyy_r;
    grad_grad_grad_psi[psiIndex][1][5] = signed_one * f3_yyz_r;
    grad_grad_grad_psi[psiIndex][1][8] = signed_one * f3_yzz_r;

    grad_grad_grad_psi[psiIndex][1][7] = grad_grad_grad_psi[psiIndex][1][5];
    grad_grad_grad_psi[psiIndex][2][4] = grad_grad_grad_psi[psiIndex][1][5];
    grad_grad_grad_psi[psiIndex][2][5] = grad_grad_grad_psi[psiIndex][1][8];
    grad_grad_grad_psi[psiIndex][2][7] = grad_grad_grad_psi[psiIndex][1][8];

    grad_grad_grad_psi[psiIndex][2][8] = signed_one * f3_zzz_r;
  }
}

template<typename ST>
void SplineR2R<ST>::evaluateVGHGH(const ParticleSet& P,
                                  const int iat,
                                  ValueVector& psi,
                                  GradVector& dpsi,
                                  HessVector& grad_grad_psi,
                                  GGGVector& grad_grad_grad_psi)
{
  const PointType& r = P.activeR(iat);
  PointType ru;
  int bc_sign = convertPos(r, ru);

#pragma omp parallel
  {
    int first, last;
    FairDivideAligned(psi.size(), getAlignment<ST>(), omp_get_num_threads(), omp_get_thread_num(), first, last);

    spline2::evaluate3d_vghgh(SplineInst->getSplinePtr(), ru, myV, myG, myH, mygH, first, last);
    assign_vghgh(bc_sign, psi, dpsi, grad_grad_psi, grad_grad_grad_psi, first, last);
  }
}

template class SplineR2R<float>;
template class SplineR2R<double>;

} // namespace qmcplusplus
