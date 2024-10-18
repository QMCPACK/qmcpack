//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////


/** @file HybridRepCplx.h
 *
 * hold HybridRepCplx
 */
#ifndef QMCPLUSPLUS_HYBRIDREP_CPLX_H
#define QMCPLUSPLUS_HYBRIDREP_CPLX_H

#include "QMCWaveFunctions/BsplineFactory/HybridRepCenterOrbitals.h"
#include "CPU/SIMD/inner_product.hpp"
#include "BsplineSet.h"

namespace qmcplusplus
{
/** hybrid representation orbitals combining B-spline orbitals on a grid and atomic centered orbitals.
 * @tparam SPLINEBASE B-spline orbital class.
 *
 * Only works with SPLINEBASE class containing complex splines
 */
template<typename SPLINEBASE>
class HybridRepCplx : public SPLINEBASE, private HybridRepCenterOrbitals<typename SPLINEBASE::DataType>
{
public:
  using SplineBase       = SPLINEBASE;
  using HYBRIDBASE       = HybridRepCenterOrbitals<typename SPLINEBASE::DataType>;
  using ST               = typename SPLINEBASE::DataType;
  using DataType         = typename SPLINEBASE::DataType;
  using PointType        = typename SPLINEBASE::PointType;
  using SingleSplineType = typename SPLINEBASE::SingleSplineType;
  using RealType         = typename SPLINEBASE::RealType;
  // types for evaluation results
  using typename SPLINEBASE::GGGVector;
  using typename SPLINEBASE::GradMatrix;
  using typename SPLINEBASE::GradType;
  using typename SPLINEBASE::GradVector;
  using typename SPLINEBASE::HessVector;
  using typename SPLINEBASE::OffloadMWVGLArray;
  using typename SPLINEBASE::ValueMatrix;
  using typename SPLINEBASE::ValueType;
  using typename SPLINEBASE::ValueVector;

private:
  using typename HYBRIDBASE::Region;

  ValueVector psi_AO, d2psi_AO;
  GradVector dpsi_AO;
  Matrix<ST, aligned_allocator<ST>> multi_myV;
  typename HYBRIDBASE::LocationSmoothingInfo info;

  using SPLINEBASE::myG;
  using SPLINEBASE::myH;
  using SPLINEBASE::myL;
  using SPLINEBASE::myV;

public:
  HybridRepCplx(const std::string& my_name) : SPLINEBASE(my_name) {}

  bool isRotationSupported() const override { return SPLINEBASE::isRotationSupported(); }
  void storeParamsBeforeRotation() override
  {
    SPLINEBASE::storeParamsBeforeRotation();
    HYBRIDBASE::storeParamsBeforeRotation();
  }
  void applyRotation(const ValueMatrix& rot_mat, bool use_stored_copy) override
  {
    SPLINEBASE::applyRotation(rot_mat, use_stored_copy);
    HYBRIDBASE::applyRotation(rot_mat, use_stored_copy);
  }

  std::string getClassName() const final { return "Hybrid" + SPLINEBASE::getClassName(); }
  std::string getKeyword() const final { return "Hybrid" + SPLINEBASE::getKeyword(); }
  bool isOMPoffload() const final { return false; }

  std::unique_ptr<SPOSet> makeClone() const override { return std::make_unique<HybridRepCplx>(*this); }

  void bcast_tables(Communicate* comm)
  {
    SPLINEBASE::bcast_tables(comm);
    HYBRIDBASE::bcast_tables(comm);
  }

  void gather_tables(Communicate* comm)
  {
    SPLINEBASE::gather_tables(comm);
    HYBRIDBASE::gather_atomic_tables(comm, SPLINEBASE::offset);
  }

  bool read_splines(hdf_archive& h5f) { return HYBRIDBASE::read_splines(h5f) && SPLINEBASE::read_splines(h5f); }

  bool write_splines(hdf_archive& h5f) { return HYBRIDBASE::write_splines(h5f) && SPLINEBASE::write_splines(h5f); }

  void evaluateValue(const ParticleSet& P, const int iat, ValueVector& psi) override
  {
    HYBRIDBASE::evaluate_v(P, iat, myV, info);
    if (info.region == Region::INTER)
      SPLINEBASE::evaluateValue(P, iat, psi);
    else if (info.region == Region::INSIDE)
      SPLINEBASE::assign_v(P.activeR(iat), myV, psi, 0, myV.size() / 2);
    else
    {
      psi_AO.resize(psi.size());
      SPLINEBASE::assign_v(P.activeR(iat), myV, psi_AO, 0, myV.size() / 2);
      SPLINEBASE::evaluateValue(P, iat, psi);
      HYBRIDBASE::interpolate_buffer_v(psi, psi_AO, info.f);
    }
  }


  void evaluateDetRatios(const VirtualParticleSet& VP,
                         ValueVector& psi,
                         const ValueVector& psiinv,
                         std::vector<ValueType>& ratios) override
  {
    if (VP.isOnSphere() && HYBRIDBASE::is_VP_batching_safe(VP))
    {
      // resize scratch space
      psi_AO.resize(psi.size());
      if (multi_myV.rows() < VP.getTotalNum())
        multi_myV.resize(VP.getTotalNum(), myV.size());
      HYBRIDBASE::evaluateValuesC2X(VP, multi_myV, info);
      for (int iat = 0; iat < VP.getTotalNum(); ++iat)
      {
        if (info.region == Region::INTER)
          SPLINEBASE::evaluateValue(VP, iat, psi);
        else if (info.region == Region::INSIDE)
        {
          Vector<ST, aligned_allocator<ST>> myV_one(multi_myV[iat], myV.size());
          SPLINEBASE::assign_v(VP.R[iat], myV_one, psi, 0, myV.size() / 2);
        }
        else
        {
          Vector<ST, aligned_allocator<ST>> myV_one(multi_myV[iat], myV.size());
          SPLINEBASE::assign_v(VP.R[iat], myV_one, psi_AO, 0, myV.size() / 2);
          SPLINEBASE::evaluateValue(VP, iat, psi);
          HYBRIDBASE::interpolate_buffer_v(psi, psi_AO, info.f);
        }
        ratios[iat] = simd::dot(psi.data(), psiinv.data(), psi.size());
      }
    }
    else
    {
      for (int iat = 0; iat < VP.getTotalNum(); ++iat)
      {
        evaluateValue(VP, iat, psi);
        ratios[iat] = simd::dot(psi.data(), psiinv.data(), psi.size());
      }
    }
  }

  void mw_evaluateDetRatios(const RefVectorWithLeader<SPOSet>& spo_list,
                            const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                            const RefVector<ValueVector>& psi_list,
                            const std::vector<const ValueType*>& invRow_ptr_list,
                            std::vector<std::vector<ValueType>>& ratios_list) const final
  {
    BsplineSet::mw_evaluateDetRatios(spo_list, vp_list, psi_list, invRow_ptr_list, ratios_list);
  }

  void evaluateVGL(const ParticleSet& P, const int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override
  {
    HYBRIDBASE::evaluate_vgl(P, iat, myV, myG, myL, info);
    if (info.region == Region::INTER)
      SPLINEBASE::evaluateVGL(P, iat, psi, dpsi, d2psi);
    else if (info.region == Region::INSIDE)
      SPLINEBASE::assign_vgl_from_l(P.activeR(iat), psi, dpsi, d2psi);
    else
    {
      psi_AO.resize(psi.size());
      dpsi_AO.resize(psi.size());
      d2psi_AO.resize(psi.size());
      SPLINEBASE::assign_vgl_from_l(P.activeR(iat), psi_AO, dpsi_AO, d2psi_AO);
      SPLINEBASE::evaluateVGL(P, iat, psi, dpsi, d2psi);
      HYBRIDBASE::interpolate_buffer_vgl(psi, dpsi, d2psi, psi_AO, dpsi_AO, d2psi_AO, info);
    }
  }

  void mw_evaluateVGL(const RefVectorWithLeader<SPOSet>& sa_list,
                      const RefVectorWithLeader<ParticleSet>& P_list,
                      int iat,
                      const RefVector<ValueVector>& psi_v_list,
                      const RefVector<GradVector>& dpsi_v_list,
                      const RefVector<ValueVector>& d2psi_v_list) const final
  {
    BsplineSet::mw_evaluateVGL(sa_list, P_list, iat, psi_v_list, dpsi_v_list, d2psi_v_list);
  }

  void mw_evaluateVGLandDetRatioGrads(const RefVectorWithLeader<SPOSet>& spo_list,
                                      const RefVectorWithLeader<ParticleSet>& P_list,
                                      int iat,
                                      const std::vector<const ValueType*>& invRow_ptr_list,
                                      OffloadMWVGLArray& phi_vgl_v,
                                      std::vector<ValueType>& ratios,
                                      std::vector<GradType>& grads) const final
  {
    BsplineSet::mw_evaluateVGLandDetRatioGrads(spo_list, P_list, iat, invRow_ptr_list, phi_vgl_v, ratios, grads);
  }

  void evaluateVGH(const ParticleSet& P,
                   const int iat,
                   ValueVector& psi,
                   GradVector& dpsi,
                   HessVector& grad_grad_psi) override
  {
    APP_ABORT("HybridRepCplx::evaluate_vgh not implemented!");
    HYBRIDBASE::evaluate_vgh(P, iat, myV, myG, myH, info);
    if (info.region == Region::INTER)
      SPLINEBASE::evaluateVGH(P, iat, psi, dpsi, grad_grad_psi);
    else
      SPLINEBASE::assign_vgh(P.activeR(iat), psi, dpsi, grad_grad_psi, 0, myV.size() / 2);
  }

  void evaluateVGHGH(const ParticleSet& P,
                     const int iat,
                     ValueVector& psi,
                     GradVector& dpsi,
                     HessVector& grad_grad_psi,
                     GGGVector& grad_grad_grad_psi) override
  {
    APP_ABORT("HybridRepCplx::evaluate_vghgh not implemented!");
  }

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) final
  {
    // bypass SPLINEBASE::evaluate_notranspose
    BsplineSet::evaluate_notranspose(P, first, last, logdet, dlogdet, d2logdet);
  }

  template<class BSPLINESPO>
  friend class HybridRepSetReader;
  template<class BSPLINESPO>
  friend class SplineSetReader;
  friend struct BsplineReader;
};

} // namespace qmcplusplus
#endif
