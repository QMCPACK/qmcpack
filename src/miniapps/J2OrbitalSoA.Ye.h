//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp. 
//                    Amrita Mathuriya, amrita.mathuriya@intel.com, Intel Corp.
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_TWOBODYJASTROW_OPTIMIZED_SOA_H
#define QMCPLUSPLUS_TWOBODYJASTROW_OPTIMIZED_SOA_H
#include "Configuration.h"
#if QMC_BUILD_LEVEL<5
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/Jastrow/DiffTwoBodyJastrowOrbital.h"
#include <qmc_common.h>
#endif
#include "Particle/DistanceTableData.h"
#include <simd/allocator.hpp>
#include <simd/algorithm.hpp>
#include <map>
#include <numeric>

namespace qmcplusplus
{

/** @ingroup OrbitalComponent
 *  @brief Specialization for two-body Jastrow function using multiple functors
 *
 * Each pair-type can have distinct function \f$u(r_{ij})\f$.
 * For electrons, distinct pair correlation functions are used
 * for spins up-up/down-down and up-down/down-up.
 *
 * Based on J2OrbitalSoA.h with these considerations
 * - DistanceTableData using SoA containers
 * - support mixed precision: FT::real_type != OHMMS_PRECISION
 * - loops over the groups: elminated PairID
 * - support simd function
 * - double the loop counts
 * - UseBuffer is introduced. Default is false.
 * - Memory use is O(N). 
 */
template<class FT>
struct  J2OrbitalSoA : public OrbitalBase
{
  ///alias FuncType
  using FuncType=FT;
  ///type of each component U, dU, d2U;
  using valT=typename FT::real_type;
  ///element position type
  using posT=TinyVector<valT,OHMMS_DIM>;
  ///use the same container 
  using RowContainer=DistanceTableData::RowContainer;

  ///number of particles
  int N;
  ///number of groups of the target particleset
  int NumGroups;
  ///task id
  int TaskID;
  ///Used to compute correction
  bool FirstTime;
  ///diff value
  RealType DiffVal;
  ///log value
  RealType LogValue;
  ///Correction
  RealType KEcorr;
  ///\f$Uat[i] = sum_(j) u_{i,j}\f$
  ParticleAttrib<valT> Uat;
  ///\f$dUat[i] = sum_(j) du_{i,j}\f$
  ParticleAttrib<posT> dUat;
  valT *FirstAddressOfdU, *LastAddressOfdU;
  ///\f$d2Uat[i] = sum_(j) d2u_{i,j}\f$
  ParticleAttrib<valT> d2Uat;
  valT cur_Uat;
  aligned_vector<valT> cur_u, cur_du, cur_d2u;
  aligned_vector<valT> old_u, old_du, old_d2u;
  aligned_vector<valT> DistCompressed;
  aligned_vector<int> DistIndice;
  ///Container for \f$F[ig*NumGroups+jg]\f$
  std::vector<FT*> F;
  ///Uniquue J2 set for cleanup
  std::map<std::string,FT*> J2Unique;

  J2OrbitalSoA(ParticleSet& p, int tid);
  J2OrbitalSoA(const J2OrbitalSoA& rhs)=delete;
  ~J2OrbitalSoA();

  /* initialize storage */
  void init(ParticleSet& p);

  /** add functor for (ia,ib) pair */
  void addFunc(int ia, int ib, FT* j);


  void resetTargetParticleSet(ParticleSet& P)
  {
    if(dPsi)
      dPsi->resetTargetParticleSet(P);
  }

  /** check in an optimizable parameter
   * @param o a super set of optimizable variables
   */
  void checkInVariables(opt_variables_type& active)
  {
  }

  /** check out optimizable variables
   */
  void checkOutVariables(const opt_variables_type& active)
  {
  }

  ///reset the value of all the unique Two-Body Jastrow functions
  void resetParameters(const opt_variables_type& active)
  {
  }

  /** print the state, e.g., optimizables */
  void reportStatus(std::ostream& os)
  {
  }
  RealType ChiesaKEcorrection() { return RealType();}
  /**@} */

  OrbitalBasePtr makeClone(ParticleSet& tqp) const;

  RealType evaluateLog(ParticleSet& P,
                       ParticleSet::ParticleGradient_t& G,
                       ParticleSet::ParticleLaplacian_t& L);

  /** recompute internal data assuming distance table is fully ready */
  void recompute(ParticleSet& P);

  void evaluateLogAndStore(ParticleSet& P,
                       ParticleSet::ParticleGradient_t& G,
                       ParticleSet::ParticleLaplacian_t& L)
  {
     RealType x=evaluateLog(P,G,L);
  }

  ValueType evaluate(ParticleSet& P,
                     ParticleSet::ParticleGradient_t& G,
                     ParticleSet::ParticleLaplacian_t& L)
  {
    evaluateLog(P,G,L);
    return std::exp(LogValue);
  }

  ValueType ratio(ParticleSet& P, int iat);
  GradType evalGrad(ParticleSet& P, int iat);
  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);
  void acceptMove(ParticleSet& P, int iat);
  inline void restore(int iat) {}

  /** compute G and L after the sweep
   */
  void evaluateGL(ParticleSet& P,
                     ParticleSet::ParticleGradient_t& G,
                     ParticleSet::ParticleLaplacian_t& L, bool fromscratch=false);

  inline RealType registerData(ParticleSet& P, PooledData<RealType>& buf)
  {
    evaluateLog(P,P.G,P.L);
    buf.add(Uat.begin(), Uat.end());
    buf.add(FirstAddressOfdU,LastAddressOfdU);
    buf.add(d2Uat.begin(), d2Uat.end());
    return LogValue;
  }

  inline void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
  {
    buf.get(Uat.begin(), Uat.end());
    buf.get(FirstAddressOfdU,LastAddressOfdU);
    buf.get(d2Uat.begin(), d2Uat.end());
  }

  RealType updateBuffer(ParticleSet& P, PooledData<RealType>& buf, bool fromscratch=false)
  {
    evaluateGL(P, P.G, P.L, false);
    buf.put(Uat.begin(), Uat.end());
    buf.put(FirstAddressOfdU,LastAddressOfdU);
    buf.put(d2Uat.begin(), d2Uat.end());
    return LogValue;
  }

  inline RealType evaluateLog(ParticleSet& P, PooledData<RealType>& buf)
  {
    buf.put(Uat.begin(), Uat.end());
    buf.put(FirstAddressOfdU,LastAddressOfdU);
    buf.put(d2Uat.begin(), d2Uat.end());
    return LogValue;
  }

  //to be removed from QMCPACK: these are not used anymore with PbyPFast
  void update(ParticleSet& P,
              ParticleSet::ParticleGradient_t& dG,
              ParticleSet::ParticleLaplacian_t& dL,
              int iat) {}
  ValueType ratio(ParticleSet& P, int iat,
      ParticleSet::ParticleGradient_t& dG,
      ParticleSet::ParticleLaplacian_t& dL){ return ValueType(1);}

  
  /*@{ internal compute engines*/
  inline void computeU3(ParticleSet& P, int iat, const RealType* restrict dist,
      RealType* restrict u, RealType* restrict du, RealType* restrict d2u);

  /** compute gradient
   */
  inline posT accumulateG(const valT* restrict du, const RowContainer& displ) const
  {
    posT grad;
    for(int idim=0; idim<OHMMS_DIM; ++idim)
    {
      const valT* restrict dX=displ.data(idim);
      valT s=valT();

      ASSUME_ALIGNED(du);
      ASSUME_ALIGNED(dX);
      #pragma omp simd reduction(+:s)
      for(int jat=0; jat<N; ++jat) s+=du[jat]*dX[jat];
      grad[idim]=s;
    }
    return grad;
  }

  /** compute gradient and lap
   */
  inline void accumulateGL(const valT* restrict du, const valT* restrict d2u,
      const RowContainer& displ, posT& grad, valT& lap) const
  {
    constexpr valT lapfac=OHMMS_DIM-RealType(1);
    lap=valT(0);
//#pragma omp simd reduction(+:lap)
    for(int jat=0; jat<N; ++jat)
      lap+=d2u[jat]+lapfac*du[jat];
    for(int idim=0; idim<OHMMS_DIM; ++idim)
    {
      const valT* restrict dX=displ.data(idim);
      valT s=valT();
//#pragma omp simd reduction(+:s)
      for(int jat=0; jat<N; ++jat) s+=du[jat]*dX[jat];
      grad[idim]=s;
    }
  }

};

template<typename FT>
J2OrbitalSoA<FT>::J2OrbitalSoA(ParticleSet& p, int tid) : TaskID(tid)
{
  init(p);
  FirstTime =true;
  KEcorr=0.0;
  //OrbitalName = "J2OrbitalSoA";
}

template<typename FT>
J2OrbitalSoA<FT>::~J2OrbitalSoA()
{ 
  auto it=J2Unique.begin();
  while(it != J2Unique.end())
  {
    delete ((*it).second);
    ++it;
  }
}//need to clean up J2Unique 

template<typename FT>
void J2OrbitalSoA<FT>::init(ParticleSet& p)
{
  N=p.getTotalNum();
  NumGroups=p.groups();

  Uat.resize(N); 
  dUat.resize(N);
  FirstAddressOfdU = &(dUat[0][0]);
  LastAddressOfdU = FirstAddressOfdU + dUat.size()*OHMMS_DIM;
  d2Uat.resize(N);
  cur_u.resize(N);
  cur_du.resize(N);
  cur_d2u.resize(N);
  old_u.resize(N);
  old_du.resize(N);
  old_d2u.resize(N);
  F.resize(NumGroups*NumGroups,nullptr);
  DistCompressed.resize(N);
  DistIndice.resize(N);
}

template<typename FT>
void J2OrbitalSoA<FT>::addFunc(int ia, int ib, FT* j)
{
  if(ia==ib)
  {
    if(ia==0)//first time, assign everything
    {
      int ij=0;
      for(int ig=0; ig<NumGroups; ++ig)
        for(int jg=0; jg<NumGroups; ++jg, ++ij)
          if(F[ij]==nullptr) F[ij]=j;
    }
  }
  else
  {
    F[ia*NumGroups+ib]=j;
    if(ia<ib)
      F[ib*NumGroups+ia]=j;
  }
  std::stringstream aname;
  aname<<ia<<ib;
  J2Unique[aname.str()]=j;
  //ChiesaKEcorrection();
  FirstTime = false;
}

template<typename FT>
OrbitalBasePtr J2OrbitalSoA<FT>::makeClone(ParticleSet& tqp) const
{
  J2OrbitalSoA<FT>* j2copy=new J2OrbitalSoA<FT>(tqp,-1);
  if (dPsi)
    j2copy->dPsi = dPsi->makeClone(tqp);
  std::map<const FT*,FT*> fcmap;
  for(int ig=0; ig<NumGroups; ++ig)
    for(int jg=ig; jg<NumGroups; ++jg)
    {
      int ij=ig*NumGroups+jg;
      if(F[ij]==0)
        continue;
      typename std::map<const FT*,FT*>::iterator fit=fcmap.find(F[ij]);
      if(fit == fcmap.end())
      {
        FT* fc=new FT(*F[ij]);
        j2copy->addFunc(ig,jg,fc);
        //if (dPsi) (j2copy->dPsi)->addFunc(aname.str(),ig,jg,fc);
        fcmap[F[ij]]=fc;
      }
    }
  j2copy->Optimizable = Optimizable;
  return j2copy;
}

/** intenal function to compute \f$\sum_j u(r_j), du/dr, d2u/dr2\f$
 * @param P particleset
 * @param iat particle index
 * @param dist starting distance
 * @param u starting value
 * @param du starting first deriv
 * @param d2u starting second deriv
 */
template<typename FT>
inline void
J2OrbitalSoA<FT>::computeU3(ParticleSet& P, int iat, const RealType* restrict dist,
    RealType* restrict u, RealType* restrict du, RealType* restrict d2u)
{
  const int igt=P.GroupID[iat]*NumGroups;
  for(int jg=0; jg<NumGroups; ++jg)
  {
    const FuncType& f2(*F[igt+jg]);
    int iStart = P.first(jg);
    int iEnd = P.last(jg);
    f2.evaluateVGL(iStart, iEnd, dist, u, du, d2u, DistCompressed.data(), DistIndice.data());
  }
  constexpr valT czero(0);
  u[iat]=czero;
  du[iat]=czero;
  d2u[iat]=czero;
}

template<typename FT>
typename  J2OrbitalSoA<FT>::ValueType 
J2OrbitalSoA<FT>::ratio(ParticleSet& P, int iat)
{
  //only ratio, ready to compute it again
  UpdateMode=ORB_PBYP_RATIO;

  const DistanceTableData* d_table=P.DistTables[0];
  const auto dist=d_table->Temp_r.data();
  cur_Uat=valT(0);
  const int igt=P.GroupID[iat]*NumGroups;
  for(int jg=0; jg<NumGroups; ++jg)
  {
    const FuncType& f2(*F[igt+jg]);
    int iStart = P.first(jg);
    int iEnd = P.last(jg);
    cur_Uat += f2.evaluateV( iStart, iEnd, dist, DistCompressed.data() );
  }

  return std::exp(Uat[iat]-cur_Uat);
}

template<typename FT>
typename  J2OrbitalSoA<FT>::GradType 
J2OrbitalSoA<FT>::evalGrad(ParticleSet& P, int iat)
{
  return GradType(dUat[iat]);
}

template<typename FT>
typename  J2OrbitalSoA<FT>::ValueType
J2OrbitalSoA<FT>::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{

  UpdateMode=ORB_PBYP_PARTIAL;

  computeU3(P,iat,P.DistTables[0]->Temp_r.data(), cur_u.data(),cur_du.data(),cur_d2u.data());
  cur_Uat=simd::accumulate_n(cur_u.data(),N,valT());
  DiffVal=Uat[iat]-cur_Uat;
  grad_iat+=accumulateG(cur_du.data(),P.DistTables[0]->Temp_dr);
  return std::exp(DiffVal);
}

template<typename FT>
void
J2OrbitalSoA<FT>::acceptMove(ParticleSet& P, int iat)
{
  // get the old u, du, d2u
  const DistanceTableData* d_table=P.DistTables[0];
  computeU3(P,iat,d_table->Distances[iat],old_u.data(),old_du.data(),old_d2u.data());
  if(UpdateMode == ORB_PBYP_RATIO)
  {//ratio-only during the move; need to compute derivatives
    const auto dist=d_table->Temp_r.data();
    computeU3(P,iat,dist,cur_u.data(),cur_du.data(),cur_d2u.data());
  }

  valT cur_d2Uat(0);
  posT cur_dUat;
  for(int jat=0; jat<N; jat++)
  {
    posT dg, newg;
    valT dl, newl, du;
    constexpr valT lapfac=OHMMS_DIM-RealType(1);
    du = cur_u[jat] - old_u[jat];
    for(int idim=0; idim<OHMMS_DIM; ++idim)
    {
      newg[idim] = cur_du[jat]*d_table->Temp_dr.data(idim)[jat];
      dg[idim]   = newg[idim] - old_du[jat]*d_table->Displacements[iat].data(idim)[jat];
    }
    newl = cur_d2u[jat] + lapfac*cur_du[jat];
    dl   = old_d2u[jat] + lapfac*old_du[jat] - newl;
    Uat[jat]   += du;
    dUat[jat]  -= dg;
    d2Uat[jat] += dl;
    cur_dUat   += newg;
    cur_d2Uat  -= newl;
  }
  Uat[iat]   = cur_Uat;
  dUat[iat]  = cur_dUat;
  d2Uat[iat] = cur_d2Uat;
}

template<typename FT>
void
J2OrbitalSoA<FT>::recompute(ParticleSet& P)
{
  const DistanceTableData* d_table=P.DistTables[0];
  for(int ig=0; ig<NumGroups; ++ig)
  {
    const int igt=ig*NumGroups;
    for(int iat=P.first(ig),last=P.last(ig); iat<last; ++iat)
    {
      computeU3(P,iat,d_table->Distances[iat],cur_u.data(),cur_du.data(),cur_d2u.data());
      Uat[iat]=simd::accumulate_n(cur_u.data(),N,valT());
      posT grad;
      valT lap;
      accumulateGL(cur_du.data(),cur_d2u.data(),d_table->Displacements[iat],grad,lap);
      dUat[iat]=grad;
      d2Uat[iat]=-lap;
    }
  }
}

template<typename FT>
typename J2OrbitalSoA<FT>::RealType
J2OrbitalSoA<FT>::evaluateLog(ParticleSet& P,
    ParticleSet::ParticleGradient_t& dG,
    ParticleSet::ParticleLaplacian_t& dL)
{
  evaluateGL(P,dG,dL,true);
  return LogValue;
}

template<typename FT>
void
J2OrbitalSoA<FT>::evaluateGL(ParticleSet& P,
    ParticleSet::ParticleGradient_t& G,
    ParticleSet::ParticleLaplacian_t& L, bool fromscratch)
{
  if(fromscratch) recompute(P);
  LogValue=valT(0);
  for(int iat=0; iat<N; ++iat)
  {
    LogValue += Uat[iat];
    G[iat] += dUat[iat];
    L[iat] += d2Uat[iat];
  }

  constexpr valT mhalf(-0.5);
  LogValue=mhalf*LogValue;
}

}
#endif
