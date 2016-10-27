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
 *
 * Each pair-type can have distinct function \f$u(r_{ij})\f$.
 * For electrons, distinct pair correlation functions are used
 * for spins up-up/down-down and up-down/down-up.
 *
 * Based on TwoBodyJatrowOrbital.h with these considerations
 * - DistanceTableData using SoA containers
 * - support mixed precision: FT::real_type != OHMMS_PRECISION
 * - loops over the groups: elminated PairID
 * - support simd function
 * - double the loop counts
 * - UseBuffer is introduced. Default is false.
 * - Memory use is 3*N*N+3N. 
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
  ///number of groups
  int NumGroups;
  ///task id
  int TaskID;
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
  ParticleAttrib<posT> Grad;
  ParticleAttrib<valT> Lap;
  ///Container for \f$F[ig*NumGroups+jg]\f$
  std::vector<FT*> F;

  J1OrbitalSoA(ParticleSet& ions, ParticleSet& els, int tid) : Ions(ions),TaskID(tid)
  {
    init(ions,els);
    myTableID=els.addTable(ions,DT_SOA);
  }

  J1OrbitalSoA(const J1OrbitalSoA& rhs)=delete;

  ~J1OrbitalSoA() 
  { 
    for(int i=0; i<F.size(); ++i)
      if(F[i] != nullptr) delete F[i];
  }

  /* initialize storage */
  void init(ParticleSet& ions, ParticleSet& els)
  {
    Nions=ions.getTotalNum();
    NumGroups=ions.getSpeciesSet().getTotalNum();
    F.resize(std::max(NumGroups,4),nullptr);
    if(NumGroups>1  && !ions.IsGrouped) 
    {
      NumGroups=0;
    }
    const int N=els.getTotalNum();
    Vat.resize(N); 
    Grad.resize(N);
    Lap.resize(N);

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

  RealType evaluateLog(ParticleSet& P,
                       ParticleSet::ParticleGradient_t& G,
                       ParticleSet::ParticleLaplacian_t& L)
  {
    evaluateLogAndStore(P,G,L);
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
    curAt = valT(0);
    if(P.Ready4Measure)
    {
      const valT* restrict dist=P.DistTables[myTableID]->Temp_r.data();
      if(NumGroups>0)
      {
        for(int jg=0; jg<NumGroups; ++jg)
        {
          if(F[jg]!=nullptr) 
            curAt += F[jg]->evaluateV(Ions.first(jg), Ions.last(jg), dist, DistCompressed.data() );
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
    }
    else
    {
      evaluateU3(P,iat,P.DistTables[myTableID]->Temp_r.data());
      curLap=valT();
      accumulateGL(dU.data(),d2U.data(),P.DistTables[myTableID]->Temp_dr,curGrad,curLap);
      curAt=simd::accumulate_n(U.data(),Nions,valT());
    }
    return std::exp(Vat[iat]-curAt);
  }

  inline void evaluateGL(ParticleSet& P)
  {
    const int n=P.getTotalNum();
    for(int iat=0; iat<n; ++iat) P.G[iat]+=Grad[iat];
    for(int iat=0; iat<n; ++iat) P.L[iat]-=Lap[iat];
  }

  /** compute gradient and lap
   */
  inline void accumulateGL(const valT* restrict du, const valT* restrict d2u,
      const RowContainer& displ, posT& grad, valT& lap) const
  {
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
  }

  /** compute U, dU and d2U 
   * @param P quantum particleset
   * @param iat the moving particle
   * @param dist starting address of the distances of the ions wrt the iat-th particle
   */
  inline void evaluateU3(ParticleSet& P, int iat, const valT* dist)
  {
    if(NumGroups>0)
    {//ions are grouped
      for(int jg=0; jg<NumGroups; ++jg)
      {
        if(F[jg]==nullptr) continue;
        F[jg]->evaluateVGL( Ions.first(jg), Ions.last(jg), dist, 
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
   *
   * Using Temp_r. Vat[iat], Grad[iat] and Lap[iat] are computed.
   */
  GradType evalGrad(ParticleSet& P, int iat)
  {
    evaluateU3(P,iat,P.DistTables[myTableID]->Temp_r.data());
    Grad[iat]=valT();
    Lap[iat]=valT();
    accumulateGL(dU.data(),d2U.data(),P.DistTables[myTableID]->Temp_dr,Grad[iat],Lap[iat]);
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
    evaluateU3(P,iat,P.DistTables[myTableID]->Temp_r.data());
    curLap=valT();
    accumulateGL(dU.data(),d2U.data(),P.DistTables[myTableID]->Temp_dr,curGrad,curLap);
    curAt=simd::accumulate_n(U.data(),Nions,valT());
    grad_iat+=curGrad;
    return std::exp(Vat[iat]-curAt);
  }

  /** Rejected move. Nothing to do */
  inline void restore(int iat) { }

  /** Accpted move. Update Vat[iat],Grad[iat] and Lap[iat] */
  void acceptMove(ParticleSet& P, int iat)
  {
    LogValue += Vat[iat]-curAt;
    Vat[iat]  = curAt;
    Grad[iat] = curGrad;
    Lap[iat]  = curLap;
  }


  inline void evaluateLogAndStore(ParticleSet& P, 
      ParticleSet::ParticleGradient_t& dG, ParticleSet::ParticleLaplacian_t& dL )
  {
    const int n=P.getTotalNum();
    const DistanceTableData& d_ie(*(P.DistTables[myTableID]));
    LogValue=valT();
    for(int iat=0; iat<n; ++iat)
    {
      evaluateU3(P,iat,d_ie.Distances[iat]);
      LogValue-=Vat[iat]=simd::accumulate_n(U.data(),Nions,valT());
      Grad[iat]=valT();
      Lap[iat]=valT();
      accumulateGL(dU.data(),d2U.data(),d_ie.Displacements[iat],Grad[iat],Lap[iat]);
      dG[iat]+=Grad[iat];
      dL[iat]-=Lap[iat];
    }
  }

  inline RealType registerData(ParticleSet& P, PooledData<RealType>& buf)
  {
    evaluateLogAndStore(P,P.G,P.L);
    return LogValue;
  }

  inline RealType updateBuffer(ParticleSet& P, PooledData<RealType>& buf, bool fromscratch=false)
  {
    evaluateGL(P);
    constexpr RealType mone(-1);
    LogValue=mone*simd::accumulate_n(Vat.data(), P.getTotalNum(), valT());
    return LogValue;
  }

  inline void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf) { }

  inline RealType evaluateLog(ParticleSet& P, PooledData<RealType>& buf)
  {
    evaluateLogAndStore(P,P.G,P.L);
    return LogValue;
  }

  /**@{ OrbitalBase virtual functions that are not essential for the development */
  void resetTargetParticleSet(ParticleSet& P){}
  void reportStatus(std::ostream& os){}
  void checkInVariables(opt_variables_type& active){}
  void checkOutVariables(const opt_variables_type& active){}
  void resetParameters(const opt_variables_type& active){}
  /**@} */

  /** must be removed */
  ValueType ratio(ParticleSet& P, int iat,
                  ParticleSet::ParticleGradient_t& dG,
                  ParticleSet::ParticleLaplacian_t& dL)
  {
    APP_ABORT("OrbitalBase::ratio(P,iat,dG,dL) shuold not Used")
    return 1;
  }

  /** must be removed */
  inline void update(ParticleSet& P,
                     ParticleSet::ParticleGradient_t& dG,
                     ParticleSet::ParticleLaplacian_t& dL,
                     int iat)
  {
    APP_ABORT("J1OrbitalSoA::update must not be used");
  }




};




}
#endif
/***************************************************************************
 * $RCSfile$   $Author: yeluo $
 * $Revision: 6449 $   $Date: 2015-04-18 00:22:48 -0500 (Sat, 18 Apr 2015) $
 * $Id: TwoBodyJastrowOrbital.h 6449 2015-04-18 05:22:48Z yeluo $
 ***************************************************************************/

