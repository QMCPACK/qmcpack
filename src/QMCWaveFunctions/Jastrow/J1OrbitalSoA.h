//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_ONEBODYJASTROW_OPTIMIZED_SOA_H
#define QMCPLUSPLUS_ONEBODYJASTROW_OPTIMIZED_SOA_H
#include "Configuration.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/Jastrow/DiffOneBodyJastrowOrbital.h"
#include <qmc_common.h>
#include <simd/allocator.hpp>
#include <simd/algorithm.hpp>
#include  <map>
#include  <numeric>

namespace qmcplusplus
{

/** @ingroup OrbitalComponent
 *  @brief Specialization for one-body Jastrow function using multiple functors
 */
template<class FT>
struct  J1OrbitalSoA : public OrbitalBase
{
  ///alias FuncType
  using FuncType=FT;
  ///type of each component U, dU, d2U;
  using valT=typename FT::real_type;
  ///element position type
  using posT=TinyVector<valT,OHMMS_DIM>;
  ///use the same container 
  using RowContainer=DistanceTableData::RowContainer;
  ///table index
  int myTableID;
  ///number of ions
  int Nions;
  ///number of electrons
  int Nelec;
  ///number of groups
  int NumGroups;
  ///reference to the sources (ions)
  const ParticleSet& Ions;

  valT LogValue;
  valT curAt;
  valT curLap;
  posT curGrad;

  ///\f$Vat[i] = sum_(j) u_{i,j}\f$
  aligned_vector<RealType> Vat;
  aligned_vector<valT> U, dU, d2U;
  aligned_vector<valT> DistCompressed;
  aligned_vector<int> DistIndice;
  Vector<posT> Grad;
  Vector<valT> Lap;
  ///Container for \f$F[ig*NumGroups+jg]\f$
  std::vector<FT*> F;

  J1OrbitalSoA(const ParticleSet& ions, ParticleSet& els) : Ions(ions)
  {
    initalize(els);
    myTableID=els.addTable(ions,DT_SOA);
    OrbitalName = "J1OrbitalSoA";
  }

  J1OrbitalSoA(const J1OrbitalSoA& rhs)=delete;

  ~J1OrbitalSoA() 
  { 
    for(int i=0; i<F.size(); ++i)
      if(F[i] != nullptr) delete F[i];
  }

  /* initialize storage */
  void initalize(ParticleSet& els)
  {
    Nions=Ions.getTotalNum();
    NumGroups=Ions.getSpeciesSet().getTotalNum();
    F.resize(std::max(NumGroups,4),nullptr);
    if(NumGroups>1  && !Ions.IsGrouped) 
    {
      NumGroups=0;
    }
    Nelec=els.getTotalNum();
    Vat.resize(Nelec);
    Grad.resize(Nelec);
    Lap.resize(Nelec);

    U.resize(Nions);
    dU.resize(Nions);
    d2U.resize(Nions);
    DistCompressed.resize(Nions);
    DistIndice.resize(Nions);
  }

  void addFunc(int source_type, FT* afunc, int target_type=-1)
  {
    if(F[source_type]!=nullptr) delete F[source_type];
    F[source_type]=afunc;
  }

  void recompute(ParticleSet& P)
  {
    const DistanceTableData& d_ie(*(P.DistTables[myTableID]));
    for(int iat=0; iat<Nelec; ++iat)
    {
      computeU3(P,iat,d_ie.Distances[iat]);
      Vat[iat]=simd::accumulate_n(U.data(),Nions,valT());
      Lap[iat]=accumulateGL(dU.data(),d2U.data(),d_ie.Displacements[iat],Grad[iat]);
    }
  }

  RealType evaluateLog(ParticleSet& P,
                       ParticleSet::ParticleGradient_t& G,
                       ParticleSet::ParticleLaplacian_t& L)
  {
    evaluateGL(P,G,L,true);
    return LogValue;
  }

  ValueType evaluate(ParticleSet& P,
                     ParticleSet::ParticleGradient_t& G,
                     ParticleSet::ParticleLaplacian_t& L)
  {
    return std::exp(evaluateLog(P,G,L));
  }
  
  ValueType ratio(ParticleSet& P, int iat)
  {
    UpdateMode=ORB_PBYP_RATIO;
    curAt = valT(0);
    const valT* restrict dist=P.DistTables[myTableID]->Temp_r.data();
    if(NumGroups>0)
    {
      for(int jg=0; jg<NumGroups; ++jg)
      {
        if(F[jg]!=nullptr)
          curAt += F[jg]->evaluateV(-1, Ions.first(jg), Ions.last(jg), dist, DistCompressed.data());
      }
    }
    else
    {
      for(int c=0; c<Nions; ++c)
      {
        int gid=Ions.GroupID[c];
        if(F[gid]!=nullptr) curAt += F[gid]->evaluate(dist[c]);
      }
    }

    if(!P.Ready4Measure)
    {//need to compute per atom
      computeU3(P,iat,P.DistTables[myTableID]->Distances[iat]);
      Lap[iat]=accumulateGL(dU.data(),d2U.data(),P.DistTables[myTableID]->Displacements[iat],Grad[iat]);
      Vat[iat]=simd::accumulate_n(U.data(),Nions,valT());
    }

    return std::exp(Vat[iat]-curAt);
  }

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
  {
    const valT* restrict dist=P.DistTables[myTableID]->Temp_r.data();
    curAt = valT(0);
    if(NumGroups>0)
    {
      for(int jg=0; jg<NumGroups; ++jg)
      {
        if(F[jg]!=nullptr)
          curAt += F[jg]->evaluateV(-1, Ions.first(jg), Ions.last(jg), dist, DistCompressed.data());
      }
    }
    else
    {
      for(int c=0; c<Nions; ++c)
      {
        int gid=Ions.GroupID[c];
        if(F[gid]!=nullptr) curAt += F[gid]->evaluate(dist[c]);
      }
    }

    for(int i=0; i<Nelec; ++i)
      ratios[i]=std::exp(Vat[i]-curAt);
  }

  inline void evaluateGL(ParticleSet& P,
                         ParticleSet::ParticleGradient_t& G,
                         ParticleSet::ParticleLaplacian_t& L,
                         bool fromscratch=false)
  {
    if(fromscratch) recompute(P);

    for(size_t iat=0; iat<Nelec; ++iat) G[iat]+=Grad[iat];
    for(size_t iat=0; iat<Nelec; ++iat) L[iat]-=Lap[iat];
    LogValue=-simd::accumulate_n(Vat.data(), Nelec, valT());
  }

  /** compute gradient and lap
   * @return lap
   */
  inline valT accumulateGL(const valT* restrict du, const valT* restrict d2u,
      const RowContainer& displ, posT& grad) const
  {
    valT lap(0);
    constexpr valT lapfac=OHMMS_DIM-RealType(1);
//#pragma omp simd reduction(+:lap)
    for(int jat=0; jat<Nions; ++jat)
      lap+=d2u[jat]+lapfac*du[jat];
    for(int idim=0; idim<OHMMS_DIM; ++idim)
    {
      const valT* restrict dX=displ.data(idim);
      valT s=valT();
//#pragma omp simd reduction(+:s)
      for(int jat=0; jat<Nions; ++jat) s+=du[jat]*dX[jat];
      grad[idim]=s;
    }
    return lap;
  }

  /** compute U, dU and d2U 
   * @param P quantum particleset
   * @param iat the moving particle
   * @param dist starting address of the distances of the ions wrt the iat-th particle
   */
  inline void computeU3(ParticleSet& P, int iat, const valT* dist)
  {
    if(NumGroups>0)
    {//ions are grouped
      CONSTEXPR valT czero(0);
      std::fill_n(U.data(),Nions,czero);
      std::fill_n(dU.data(),Nions,czero);
      std::fill_n(d2U.data(),Nions,czero);

      for(int jg=0; jg<NumGroups; ++jg)
      {
        if(F[jg]==nullptr) continue;
        F[jg]->evaluateVGL(-1, Ions.first(jg), Ions.last(jg), dist,
            U.data(), dU.data(), d2U.data(), DistCompressed.data(), DistIndice.data());
      }
    }
    else
    {
      for(int c=0; c<Nions; ++c)
      {
        int gid=Ions.GroupID[c];
        if(F[gid]!=nullptr)
        {
          U[c]= F[gid]->evaluate(dist[c],dU[c],d2U[c]);
          dU[c]/=dist[c];
        }
      }
    }
  }

  /** compute the gradient during particle-by-particle update
   * @param P quantum particleset
   * @param iat particle index
   */
  GradType evalGrad(ParticleSet& P, int iat)
  {
    computeU3(P,iat,P.DistTables[myTableID]->Distances[iat]);
    Lap[iat]=accumulateGL(dU.data(),d2U.data(),P.DistTables[myTableID]->Displacements[iat],Grad[iat]);
    Vat[iat]=simd::accumulate_n(U.data(),Nions,valT());
    return GradType(Grad[iat]);
  }

  /** compute the gradient during particle-by-particle update
   * @param P quantum particleset
   * @param iat particle index
   *
   * Using Temp_r. curAt, curGrad and curLap are computed.
   */
  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
  {
    UpdateMode=ORB_PBYP_PARTIAL;

    computeU3(P,iat,P.DistTables[myTableID]->Temp_r.data());
    curLap=accumulateGL(dU.data(),d2U.data(),P.DistTables[myTableID]->Temp_dr,curGrad);
    curAt=simd::accumulate_n(U.data(),Nions,valT());
    grad_iat+=curGrad;
    return std::exp(Vat[iat]-curAt);
  }

  /** Rejected move. Nothing to do */
  inline void restore(int iat) { }

  /** Accpted move. Update Vat[iat],Grad[iat] and Lap[iat] */
  void acceptMove(ParticleSet& P, int iat)
  {

    if(UpdateMode == ORB_PBYP_RATIO)
    {
      computeU3(P,iat,P.DistTables[myTableID]->Temp_r.data());
      curLap=accumulateGL(dU.data(),d2U.data(),P.DistTables[myTableID]->Temp_dr,curGrad);
    }

    LogValue += Vat[iat]-curAt;
    Vat[iat]  = curAt;
    Grad[iat] = curGrad;
    Lap[iat]  = curLap;
  }


  inline void registerData(ParticleSet& P, WFBufferType& buf) { }

  inline RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false)
  {
    evaluateGL(P, P.G, P.L, false);
    return LogValue;
  }

  inline void copyFromBuffer(ParticleSet& P, WFBufferType& buf) { }

  OrbitalBasePtr makeClone(ParticleSet& tqp) const
  {
    J1OrbitalSoA<FT>* j1copy=new J1OrbitalSoA<FT>(Ions,tqp);
    j1copy->Optimizable=Optimizable;
    for (size_t i=0, n=F.size(); i<n; ++i)
    {
      if (F[i] != nullptr) j1copy->addFunc(i,new FT(*F[i]));
    }
    if (dPsi)
    {
      j1copy->dPsi =  dPsi->makeClone(tqp);
    }
    return j1copy;
  }

  /**@{ OrbitalBase virtual functions that are not essential for the development */
  void resetTargetParticleSet(ParticleSet& P){}
  void reportStatus(std::ostream& os)
  {
    for (size_t i=0,n=F.size(); i<n; ++i)
    {
      if(F[i] != nullptr) F[i]->myVars.print(os);
    }
  }

  void checkInVariables(opt_variables_type& active)
  {
    myVars.clear();
    for (size_t i=0,n=F.size(); i<n; ++i)
    {
      if(F[i] != nullptr)
      {
        F[i]->checkInVariables(active);
        F[i]->checkInVariables(myVars);
      }
    }
  }
  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
    Optimizable=myVars.is_optimizable();
    for (size_t i=0,n=F.size(); i<n; ++i)
      if (F[i] != nullptr) F[i]->checkOutVariables(active);
    if (dPsi)
      dPsi->checkOutVariables(active);
  }

  void resetParameters(const opt_variables_type& active)
  {
    if (!Optimizable)
      return;
    for (size_t i=0,n=F.size(); i<n; ++i)
      if (F[i] != nullptr) F[i]->resetParameters(active);

    for (int i=0; i<myVars.size(); ++i)
    {
      int ii=myVars.Index[i];
      if (ii>=0)
        myVars[i]= active[ii];
    }
    if (dPsi)
      dPsi->resetParameters(active);
  }
  /**@} */

};




}
#endif
