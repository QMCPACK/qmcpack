//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
//////////////////////////////////////////////////////////////////////////////////////


/** @file HybridRealAdoptor.h
 *
 * Adoptor classes to handle real hybrid orbitals with arbitrary precision
 */
#ifndef QMCPLUSPLUS_HYBRID_REAL_SOA_ADOPTOR_H
#define QMCPLUSPLUS_HYBRID_REAL_SOA_ADOPTOR_H

#include <QMCWaveFunctions/BsplineFactory/HybridAdoptorBase.h>
namespace qmcplusplus
{

/** adoptor class to match
 *
 */
template<typename BaseAdoptor>
struct HybridRealSoA: public BaseAdoptor, public HybridAdoptorBase<typename BaseAdoptor::DataType>
{
  using HybridBase       = HybridAdoptorBase<typename BaseAdoptor::DataType>;
  using PointType        = typename BaseAdoptor::PointType;
  using SingleSplineType = typename BaseAdoptor::SingleSplineType;
  using RealType         = typename SPOSetBase::RealType;
  using ValueType        = typename SPOSetBase::ValueType;

  typename OrbitalSetTraits<ValueType>::ValueVector_t psi_AO, d2psi_AO;
  typename OrbitalSetTraits<ValueType>::GradVector_t dpsi_AO;

  using BaseAdoptor::myV;
  using BaseAdoptor::myG;
  using BaseAdoptor::myL;
  using BaseAdoptor::myH;
  using BaseAdoptor::HalfG;
  using BaseAdoptor::PrimLattice;
  using HybridBase::dist_r;
  using HybridBase::dist_dr;
  using HybridBase::df_dr;
  using HybridBase::d2f_dr2;

  HybridRealSoA(): BaseAdoptor()
  {
    this->AdoptorName="Hybrid"+this->AdoptorName;
    this->KeyWord="Hybrid"+this->KeyWord;
  }

  inline void resizeStorage(size_t n, size_t nvals)
  {
    BaseAdoptor::resizeStorage(n,nvals);
    HybridBase::resizeStorage(myV.size());
  }

  void bcast_tables(Communicate* comm)
  {
    BaseAdoptor::bcast_tables(comm);
    HybridBase::bcast_tables(comm);
  }

  void reduce_tables(Communicate* comm)
  {
    //BaseAdoptor::reduce_tables(comm);
    BaseAdoptor::gather_tables(comm);
    HybridBase::reduce_atomic_tables(comm);
  }

  inline void flush_zero()
  {
    BaseAdoptor::flush_zero();
    HybridBase::flush_zero();
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
    const RealType smooth_factor=HybridBase::evaluate_v(P,iat,myV);
    const RealType cone(1);
    if(smooth_factor<0)
    {
      BaseAdoptor::evaluate_v(P,iat,psi);
    }
    else if (smooth_factor==cone)
    {
      const PointType& r=P.R[iat];
      int bc_sign=HybridBase::get_bc_sign(r, PrimLattice, HalfG);
      BaseAdoptor::assign_v(bc_sign,psi);
    }
    else
    {
      const PointType& r=P.R[iat];
      psi_AO.resize(psi.size());
      int bc_sign=HybridBase::get_bc_sign(r, PrimLattice, HalfG);
      BaseAdoptor::assign_v(bc_sign,psi_AO);
      BaseAdoptor::evaluate_v(P,iat,psi);
      for(size_t i=0; i<psi.size(); i++)
        psi[i] = psi_AO[i]*smooth_factor + psi[i]*(cone-smooth_factor);
    }
  }

  template<typename VV, typename GV>
  inline void evaluate_vgl(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, VV& d2psi)
  {
    const RealType smooth_factor=HybridBase::evaluate_vgl(P,iat,myV,myG,myL);
    const RealType cone(1);
    if(smooth_factor<0)
    {
      BaseAdoptor::evaluate_vgl(P,iat,psi,dpsi,d2psi);
    }
    else if(smooth_factor==cone)
    {
      const PointType& r=P.R[iat];
      int bc_sign=HybridBase::get_bc_sign(r, PrimLattice, HalfG);
      BaseAdoptor::assign_vgl_from_l(bc_sign,psi,dpsi,d2psi);
    }
    else
    {
      const PointType& r=P.R[iat];
      const RealType ctwo(2);
      const RealType rinv(1.0/dist_r);
      psi_AO.resize(psi.size());
      dpsi_AO.resize(psi.size());
      d2psi_AO.resize(psi.size());
      int bc_sign=HybridBase::get_bc_sign(r, PrimLattice, HalfG);
      BaseAdoptor::assign_vgl_from_l(bc_sign,psi_AO,dpsi_AO,d2psi_AO);
      BaseAdoptor::evaluate_vgl(P,iat,psi,dpsi,d2psi);
      for(size_t i=0; i<psi.size(); i++)
      {
        d2psi[i] = d2psi_AO[i]*smooth_factor + d2psi[i]*(cone-smooth_factor)
                 + df_dr * rinv * ctwo * dot(dpsi[i]-dpsi_AO[i], dist_dr)
                 + (psi_AO[i]-psi[i]) * (d2f_dr2 + ctwo * rinv *df_dr);
         dpsi[i] =  dpsi_AO[i]*smooth_factor +  dpsi[i]*(cone-smooth_factor)
                 + df_dr * rinv * dist_dr * (psi[i]-psi_AO[i]);
          psi[i] =   psi_AO[i]*smooth_factor +   psi[i]*(cone-smooth_factor);
      }
    }
  }

  /** evaluate VGL using VectorSoaContainer
   * @param r position
   * @param psi value container
   * @param dpsi gradient-laplacian container
   */
  template<typename VGL>
  inline void evaluate_vgl_combo(const ParticleSet& P, const int iat, VGL& vgl)
  {
    APP_ABORT("HybridRealSoA::evaluate_vgl_combo not implemented!");
    if(HybridBase::evaluate_vgh(P,iat,myV,myG,myH))
    {
      const PointType& r=P.R[iat];
      int bc_sign=HybridBase::get_bc_sign(r, PrimLattice, HalfG);
      BaseAdoptor::assign_vgl_soa(bc_sign,vgl);
    }
    else
      BaseAdoptor::evaluate_vgl_combo(P,iat,vgl);
  }

  template<typename VV, typename GV, typename GGV>
  inline void evaluate_vgh(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    APP_ABORT("HybridRealSoA::evaluate_vgh not implemented!");
    if(HybridBase::evaluate_vgh(P,iat,myV,myG,myH))
    {
      const PointType& r=P.R[iat];
      int bc_sign=HybridBase::get_bc_sign(r, PrimLattice, HalfG);
      BaseAdoptor::assign_vgh(bc_sign,psi,dpsi,grad_grad_psi);
    }
    else
      BaseAdoptor::evaluate_vgh(P,iat,psi,dpsi,grad_grad_psi);
  }
};

}
#endif
