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


/** @file HybridRepReal.h
 *
 * hold HybridRepReal
 */
#ifndef QMCPLUSPLUS_HYBRIDREP_REAL_H
#define QMCPLUSPLUS_HYBRIDREP_REAL_H

#include "QMCWaveFunctions/BsplineFactory/HybridRepCenterOrbitals.h"
namespace qmcplusplus
{
/** hybrid representation orbitals combining B-spline orbitals on a grid and atomic centered orbitals.
 * @tparam SPLINEBASE B-spline orbital class.
 *
 * Only works with SPLINEBASE class containing real splines
 */
template<typename SPLINEBASE>
class HybridRepReal : public SPLINEBASE, private HybridRepCenterOrbitals<typename SPLINEBASE::DataType>
{
public:
  using HYBRIDBASE       = HybridRepCenterOrbitals<typename SPLINEBASE::DataType>;
  using ST               = typename SPLINEBASE::DataType;
  using PointType        = typename SPLINEBASE::PointType;
  using SingleSplineType = typename SPLINEBASE::SingleSplineType;
  using RealType         = typename SPLINEBASE::RealType;
  // types for evaluation results
  using typename SPLINEBASE::GGGVector;
  using typename SPLINEBASE::GradVector;
  using typename SPLINEBASE::HessVector;
  using typename SPLINEBASE::ValueType;
  using typename SPLINEBASE::ValueVector;

private:
  ValueVector psi_AO, d2psi_AO;
  GradVector dpsi_AO;
  Matrix<ST, aligned_allocator<ST>> multi_myV;

  using SPLINEBASE::HalfG;
  using SPLINEBASE::myG;
  using SPLINEBASE::myH;
  using SPLINEBASE::myL;
  using SPLINEBASE::myV;
  using SPLINEBASE::PrimLattice;

public:
  HybridRepReal()
  {
    this->className = "Hybrid" + this->className;
    this->KeyWord   = "Hybrid" + this->KeyWord;
  }

  std::unique_ptr<SPOSet> makeClone() const override { return std::make_unique<HybridRepReal>(*this); }

  inline void resizeStorage(size_t n, size_t nvals)
  {
    SPLINEBASE::resizeStorage(n, nvals);
    HYBRIDBASE::resizeStorage(myV.size());
  }

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

  inline void flush_zero()
  {
    //SPLINEBASE::flush_zero();
    HYBRIDBASE::flush_zero();
  }

  bool read_splines(hdf_archive& h5f) { return HYBRIDBASE::read_splines(h5f) && SPLINEBASE::read_splines(h5f); }

  bool write_splines(hdf_archive& h5f) { return HYBRIDBASE::write_splines(h5f) && SPLINEBASE::write_splines(h5f); }

  void evaluateValue(const ParticleSet& P, const int iat, ValueVector& psi) override
  {
    const RealType smooth_factor = HYBRIDBASE::evaluate_v(P, iat, myV);
    const RealType cone(1);
    if (smooth_factor < 0)
    {
      SPLINEBASE::evaluateValue(P, iat, psi);
    }
    else if (smooth_factor == cone)
    {
      const PointType& r = P.activeR(iat);
      int bc_sign        = HYBRIDBASE::get_bc_sign(r, PrimLattice, HalfG);
      SPLINEBASE::assign_v(bc_sign, myV, psi, 0, myV.size());
    }
    else
    {
      const PointType& r = P.activeR(iat);
      psi_AO.resize(psi.size());
      int bc_sign = HYBRIDBASE::get_bc_sign(r, PrimLattice, HalfG);
      SPLINEBASE::assign_v(bc_sign, myV, psi_AO, 0, myV.size());
      SPLINEBASE::evaluateValue(P, iat, psi);
      HYBRIDBASE::interpolate_buffer_v(psi, psi_AO);
    }
  }

  void evaluateDetRatios(const VirtualParticleSet& VP,
                         ValueVector& psi,
                         const ValueVector& psiinv,
                         std::vector<ValueType>& ratios) override
  {
    if (VP.isOnSphere() && HYBRIDBASE::is_batched_safe(VP))
    {
      // resize scratch space
      psi_AO.resize(psi.size());
      if (multi_myV.rows() < VP.getTotalNum())
        multi_myV.resize(VP.getTotalNum(), myV.size());
      std::vector<int> bc_signs(VP.getTotalNum());
      const RealType smooth_factor = HYBRIDBASE::evaluateValuesR2R(VP, PrimLattice, HalfG, multi_myV, bc_signs);
      const RealType cone(1);
      for (int iat = 0; iat < VP.getTotalNum(); ++iat)
      {
        if (smooth_factor < 0)
          SPLINEBASE::evaluateValue(VP, iat, psi);
        else if (smooth_factor == cone)
        {
          Vector<ST, aligned_allocator<ST>> myV_one(multi_myV[iat], myV.size());
          SPLINEBASE::assign_v(bc_signs[iat], myV_one, psi, 0, myV.size());
        }
        else
        {
          Vector<ST, aligned_allocator<ST>> myV_one(multi_myV[iat], myV.size());
          SPLINEBASE::assign_v(bc_signs[iat], myV_one, psi_AO, 0, myV.size());
          SPLINEBASE::evaluateValue(VP, iat, psi);
          HYBRIDBASE::interpolate_buffer_v(psi, psi_AO);
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

  void evaluateVGL(const ParticleSet& P, const int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) override
  {
    const RealType smooth_factor = HYBRIDBASE::evaluate_vgl(P, iat, myV, myG, myL);
    const RealType cone(1);
    if (smooth_factor < 0)
    {
      SPLINEBASE::evaluateVGL(P, iat, psi, dpsi, d2psi);
    }
    else if (smooth_factor == cone)
    {
      const PointType& r = P.activeR(iat);
      int bc_sign        = HYBRIDBASE::get_bc_sign(r, PrimLattice, HalfG);
      SPLINEBASE::assign_vgl_from_l(bc_sign, psi, dpsi, d2psi);
    }
    else
    {
      const PointType& r = P.activeR(iat);
      psi_AO.resize(psi.size());
      dpsi_AO.resize(psi.size());
      d2psi_AO.resize(psi.size());
      int bc_sign = HYBRIDBASE::get_bc_sign(r, PrimLattice, HalfG);
      SPLINEBASE::assign_vgl_from_l(bc_sign, psi_AO, dpsi_AO, d2psi_AO);
      SPLINEBASE::evaluateVGL(P, iat, psi, dpsi, d2psi);
      HYBRIDBASE::interpolate_buffer_vgl(psi, dpsi, d2psi, psi_AO, dpsi_AO, d2psi_AO);
    }
  }

  void evaluateVGH(const ParticleSet& P,
                   const int iat,
                   ValueVector& psi,
                   GradVector& dpsi,
                   HessVector& grad_grad_psi) override
  {
    APP_ABORT("HybridRepReal::evaluateVGH not implemented!");
    if (HYBRIDBASE::evaluate_vgh(P, iat, myV, myG, myH))
    {
      const PointType& r = P.activeR(iat);
      int bc_sign        = HYBRIDBASE::get_bc_sign(r, PrimLattice, HalfG);
      SPLINEBASE::assign_vgh(bc_sign, psi, dpsi, grad_grad_psi, 0, myV.size());
    }
    else
      SPLINEBASE::evaluateVGH(P, iat, psi, dpsi, grad_grad_psi);
  }

  void evaluateVGHGH(const ParticleSet& P,
                     const int iat,
                     ValueVector& psi,
                     GradVector& dpsi,
                     HessVector& grad_grad_psi,
                     GGGVector& grad_grad_grad_psi) override
  {
    APP_ABORT("HybridRepCplx::evaluateVGHGH not implemented!");
  }

  template<class BSPLINESPO>
  friend class HybridRepSetReader;
  template<class BSPLINESPO>
  friend struct SplineSetReader;
  friend struct BsplineReaderBase;
};

} // namespace qmcplusplus
#endif
