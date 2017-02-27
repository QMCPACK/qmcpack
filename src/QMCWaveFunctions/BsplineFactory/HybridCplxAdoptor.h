//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory 
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory

//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file HybridCplxAdoptor.h
 *
 * Adoptor classes to handle complex hybrid orbitals with arbitrary precision
 */
#ifndef QMCPLUSPLUS_HYBRID_CPLX_SOA_ADOPTOR_H
#define QMCPLUSPLUS_HYBRID_CPLX_SOA_ADOPTOR_H

#include <QMCWaveFunctions/BsplineFactory/HybridAdoptorBase.h>
namespace qmcplusplus
{

/** adoptor class to match 
 *
 */
template<typename BaseAdoptor>
struct HybridCplxSoA: public BaseAdoptor, public HybridAdoptorBase<typename BaseAdoptor::DataType>
{
  using HybridBase       = HybridAdoptorBase<typename BaseAdoptor::DataType>;
  using ST               = typename BaseAdoptor::DataType;
  using PointType        = typename BaseAdoptor::PointType;
  using SingleSplineType = typename BaseAdoptor::SingleSplineType;

  using BaseAdoptor::PrimLattice;
  using BaseAdoptor::myV;
  using BaseAdoptor::myG;
  using BaseAdoptor::myL;
  using BaseAdoptor::myH;
  using BaseAdoptor::SplineInst;

  HybridCplxSoA(): BaseAdoptor()
  {
    this->is_complex=true;
    this->is_soa_ready=true;
    this->AdoptorName="Hybrid"+this->AdoptorName;
    this->KeyWord="Hybrid"+this->KeyWord;
  }

  bool read_splines(hdf_archive& h5f)
  {
    return BaseAdoptor::read_splines(h5f) && HybridBase::read_splines(h5f);
  }

  bool write_splines(hdf_archive& h5f)
  {
    return BaseAdoptor::write_splines(h5f) && HybridBase::write_splines(h5f);
  }

  template<typename VV>
  inline void evaluate_v(const ParticleSet& P, const int iat, VV& psi)
  {
    const PointType& r=P.R[iat];
    if(!HybridBase::evaluate_v(P,myV))
    {
      PointType ru(PrimLattice.toUnit_floor(r));
      SplineInst->evaluate(ru,myV);
    }
    BaseAdoptor::assign_v(r,psi);
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, VV& d2psi)
  {
    const PointType& r=P.R[iat];
    if(!HybridBase::evaluate_vgh(P,myV,myG,myH))
    {
      PointType ru(PrimLattice.toUnit_floor(r));
      SplineInst->evaluate_vgh(ru,myV,myG,myH);
    }
    BaseAdoptor::assign_vgl(r,psi,dpsi,d2psi);
  }

  /** evaluate VGL using VectorSoaContainer
   * @param r position
   * @param psi value container
   * @param dpsi gradient-laplacian container
   */
  template<typename VGL>
  inline void evaluate_vgl_combo(const ParticleSet& P, const int iat, VGL& vgl)
  {
    const PointType& r=P.R[iat];
    if(!HybridBase::evaluate_vgh(P,myV,myG,myH))
    {
      PointType ru(PrimLattice.toUnit_floor(r));
      SplineInst->evaluate_vgh(ru,myV,myG,myH);
    }
    BaseAdoptor::assign_vgl_soa(r,vgl);
  }

  template<typename VV, typename GV, typename GGV>
  void evaluate_vgh(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    const PointType& r=P.R[iat];
    if(!HybridBase::evaluate_vgh(P,myV,myG,myH))
    {
      PointType ru(PrimLattice.toUnit_floor(r));
      SplineInst->evaluate_vgh(ru,myV,myG,myH);
    }
    BaseAdoptor::assign_vgh(r,psi,dpsi,grad_grad_psi);
  }
};

}
#endif
