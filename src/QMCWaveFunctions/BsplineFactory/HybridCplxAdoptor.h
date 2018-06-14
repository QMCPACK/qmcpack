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
  using RealType         = typename SPOSetBase::RealType;
  using ValueType        = typename SPOSetBase::ValueType;

  typename OrbitalSetTraits<ValueType>::ValueVector_t psi_AO, d2psi_AO;
  typename OrbitalSetTraits<ValueType>::GradVector_t dpsi_AO;

  using BaseAdoptor::myV;
  using BaseAdoptor::myG;
  using BaseAdoptor::myL;
  using BaseAdoptor::myH;
  using HybridBase::dist_r;
  using HybridBase::dist_dr;
  using HybridBase::df_dr;
  using HybridBase::d2f_dr2;

  HybridCplxSoA(): BaseAdoptor()
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

  void gather_tables(Communicate* comm)
  {
    BaseAdoptor::gather_tables(comm);
    HybridBase::gather_atomic_tables(comm, this->offset_cplx, this->offset_real);
  }

  bool read_splines(hdf_archive& h5f)
  {
    return HybridBase::read_splines(h5f) && BaseAdoptor::read_splines(h5f);
  }

  bool write_splines(hdf_archive& h5f)
  {
    return HybridBase::write_splines(h5f) && BaseAdoptor::write_splines(h5f);
  }

  inline void flush_zero()
  {
    //BaseAdoptor::flush_zero();
    HybridBase::flush_zero();
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
      const PointType& r=P.activeR(iat);
      BaseAdoptor::assign_v(r,myV,psi);
    }
    else
    {
      const PointType& r=P.activeR(iat);
      psi_AO.resize(psi.size());
      BaseAdoptor::assign_v(r,myV,psi_AO);
      BaseAdoptor::evaluate_v(P,iat,psi);
      for(size_t i=0; i<psi.size(); i++)
        psi[i] = psi_AO[i]*smooth_factor + psi[i]*(cone-smooth_factor);
    }
  }


  template<typename VM>
  inline void evaluateValues(VirtualParticleSet& VP, VM& psiM)
  {
    const size_t m=psiM.cols();
    if(VP.isOnSphere())
    {
      Matrix<ST,aligned_allocator<ST> > multi_myV((ST*)VP.SPOMem.data(),VP.getTotalNum(),myV.size());
      const RealType smooth_factor=HybridBase::evaluateValuesC2X(VP,multi_myV);
      const RealType cone(1);
      if(smooth_factor<0)
      {
        for(int iat=0; iat<VP.getTotalNum(); ++iat)
        {
          Vector<SPOSetBase::ValueType> psi(psiM[iat],m);
          BaseAdoptor::evaluate_v(VP,iat,psi);
        }
      }
      else if (smooth_factor==cone)
      {
        for(int iat=0; iat<VP.getTotalNum(); ++iat)
        {
          const PointType& r=VP.R[iat];
          Vector<SPOSetBase::ValueType> psi(psiM[iat],m);
          Vector<ST,aligned_allocator<ST> > myV_one(multi_myV[iat],myV.size());
          BaseAdoptor::assign_v(r,myV_one,psi);
        }
      }
      else
      {
        psi_AO.resize(m);
        for(int iat=0; iat<VP.getTotalNum(); ++iat)
        {
          const PointType& r=VP.R[iat];
          Vector<SPOSetBase::ValueType> psi(psiM[iat],m);
          Vector<ST,aligned_allocator<ST> > myV_one(multi_myV[iat],myV.size());
          BaseAdoptor::assign_v(r,myV_one,psi_AO);
          BaseAdoptor::evaluate_v(VP,iat,psi);
          for(size_t i=0; i<psi.size(); i++)
            psi[i] = psi_AO[i]*smooth_factor + psi[i]*(cone-smooth_factor);
        }
      }
    }
    else
    {
      for(int iat=0; iat<VP.getTotalNum(); ++iat)
      {
        Vector<SPOSetBase::ValueType> psi(psiM[iat],m);
        evaluate_v(VP,iat,psi);
      }
    }
  }

  inline size_t estimateMemory(const int nP)
  {
    return BaseAdoptor::estimateMemory(nP)+myV.size()*sizeof(ST)/sizeof(ValueType)*nP;
  }

  template<typename TT>
  inline TT evaluate_dot(const ParticleSet& P, int iat, const TT* restrict arow, ST* scratch)
  {
    Vector<ST> vtmp(scratch,myV.size());
    if(HybridBase::evaluate_v(P,iat,vtmp))
      return BaseAdoptor::evaluate_dot(P,iat,arow,scratch,false);
    else
      return BaseAdoptor::evaluate_dot(P,iat,arow,scratch,true);
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
      const PointType& r=P.activeR(iat);
      BaseAdoptor::assign_vgl_from_l(r,psi,dpsi,d2psi);
    }
    else
    {
      const PointType& r=P.activeR(iat);
      const RealType ctwo(2);
      const RealType rinv(1.0/dist_r);
      psi_AO.resize(psi.size());
      dpsi_AO.resize(psi.size());
      d2psi_AO.resize(psi.size());
      BaseAdoptor::assign_vgl_from_l(r,psi_AO,dpsi_AO,d2psi_AO);
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
    APP_ABORT("HybridCplxSoA::evaluate_vgl_combo not implemented!");
    if(HybridBase::evaluate_vgh(P,iat,myV,myG,myH))
    {
      const PointType& r=P.activeR(iat);
      BaseAdoptor::assign_vgl_soa(r,vgl);
    }
    else
      BaseAdoptor::evaluate_vgl_combo(P,iat,vgl);
  }

  template<typename VV, typename GV, typename GGV>
  inline void evaluate_vgh(const ParticleSet& P, const int iat, VV& psi, GV& dpsi, GGV& grad_grad_psi)
  {
    APP_ABORT("HybridCplxSoA::evaluate_vgh not implemented!");
    if(HybridBase::evaluate_vgh(P,iat,myV,myG,myH))
    {
      const PointType& r=P.activeR(iat);
      BaseAdoptor::assign_vgh(r,psi,dpsi,grad_grad_psi);
    }
    else
      BaseAdoptor::evaluate_vgh(P,iat,psi,dpsi,grad_grad_psi);
  }
};

}
#endif
