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
  size_t N;
  ///number of particles + padded
  size_t N_padded;
  ///number of groups of the target particleset
  size_t NumGroups;
  ///task id
  int TaskID;
  ///Used to compute correction
  bool FirstTime;
  ///diff value
  RealType DiffVal;
  ///Correction
  RealType KEcorr;
  ///\f$Uat[i] = sum_(j) u_{i,j}\f$
  Vector<valT> Uat;
  ///\f$dUat[i] = sum_(j) du_{i,j}\f$
  using gContainer_type=VectorSoaContainer<valT,OHMMS_DIM>;
  gContainer_type dUat;
  ///\f$d2Uat[i] = sum_(j) d2u_{i,j}\f$
  Vector<valT> d2Uat;
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
    myVars.clear();
    typename std::map<std::string,FT*>::iterator it(J2Unique.begin()),it_end(J2Unique.end());
    while(it != it_end)
    {
      (*it).second->checkInVariables(active);
      (*it).second->checkInVariables(myVars);
      ++it;
    }
  }

  /** check out optimizable variables
   */
  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
    Optimizable=myVars.is_optimizable();
    typename std::map<std::string,FT*>::iterator it(J2Unique.begin()),it_end(J2Unique.end());
    while(it != it_end)
    {
      (*it).second->checkOutVariables(active);
      ++it;
    }
    if(dPsi)
      dPsi->checkOutVariables(active);
  }

  ///reset the value of all the unique Two-Body Jastrow functions
  void resetParameters(const opt_variables_type& active)
  {
    if(!Optimizable)
      return;
    typename std::map<std::string,FT*>::iterator it(J2Unique.begin()),it_end(J2Unique.end());
    while(it != it_end)
    {
      (*it).second->resetParameters(active);
      ++it;
    }
    if(dPsi)
      dPsi->resetParameters( active );
    for(int i=0; i<myVars.size(); ++i)
    {
      int ii=myVars.Index[i];
      if(ii>=0)
        myVars[i]= active[ii];
    }
  }

  /** print the state, e.g., optimizables */
  void reportStatus(std::ostream& os)
  {
    typename std::map<std::string,FT*>::iterator it(J2Unique.begin()),it_end(J2Unique.end());
    while(it != it_end)
    {
      (*it).second->myVars.print(os);
      ++it;
    }
    ChiesaKEcorrection();
  }
  RealType ChiesaKEcorrection() { return RealType();}
  /**@} */

  OrbitalBasePtr makeClone(ParticleSet& tqp) const;

  RealType evaluateLog(ParticleSet& P,
                       ParticleSet::ParticleGradient_t& G,
                       ParticleSet::ParticleLaplacian_t& L);

  /** recompute internal data assuming distance table is fully ready */
  void recompute(ParticleSet& P);

  ValueType evaluate(ParticleSet& P,
                     ParticleSet::ParticleGradient_t& G,
                     ParticleSet::ParticleLaplacian_t& L)
  {
    evaluateLog(P,G,L);
    return std::exp(LogValue);
  }

  ValueType ratio(ParticleSet& P, int iat);
  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios);
  GradType evalGrad(ParticleSet& P, int iat);
  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);
  void acceptMove(ParticleSet& P, int iat);
  inline void restore(int iat) {}

  /** compute G and L after the sweep
   */
  void evaluateGL(ParticleSet& P,
                     ParticleSet::ParticleGradient_t& G,
                     ParticleSet::ParticleLaplacian_t& L, bool fromscratch=false);

  inline void registerData(ParticleSet& P, WFBufferType& buf)
  {
    if ( Bytes_in_WFBuffer == 0 )
    {
      Bytes_in_WFBuffer = buf.current();
      buf.add(Uat.begin(), Uat.end());
      buf.add(dUat.data(), dUat.end());
      buf.add(d2Uat.begin(), d2Uat.end());
      Bytes_in_WFBuffer = buf.current()-Bytes_in_WFBuffer;
      // free local space
      Uat.free();
      dUat.free();
      d2Uat.free();
    }
    else
    {
      buf.forward(Bytes_in_WFBuffer);
    }
  }

  inline void copyFromBuffer(ParticleSet& P, WFBufferType& buf)
  {
    Uat.attachReference(buf.lendReference<valT>(N), N);
    dUat.attachReference(N, N_padded, buf.lendReference<valT>(N_padded*OHMMS_DIM));
    d2Uat.attachReference(buf.lendReference<valT>(N), N);
  }

  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false)
  {
    evaluateGL(P, P.G, P.L, false);
    buf.forward(Bytes_in_WFBuffer);
    return LogValue;
  }

  /*@{ internal compute engines*/
  inline void computeU3(ParticleSet& P, int iat, const RealType* restrict dist,
      RealType* restrict u, RealType* restrict du, RealType* restrict d2u, bool triangle=false);

  /** compute gradient
   */
  inline posT accumulateG(const valT* restrict du, const RowContainer& displ) const
  {
    posT grad;
    for(int idim=0; idim<OHMMS_DIM; ++idim)
    {
      const valT* restrict dX=displ.data(idim);
      valT s=valT();

      #pragma omp simd reduction(+:s) aligned(du,dX)
      for(int jat=0; jat<N; ++jat) s+=du[jat]*dX[jat];
      grad[idim]=s;
    }
    return grad;
  }

};

template<typename FT>
J2OrbitalSoA<FT>::J2OrbitalSoA(ParticleSet& p, int tid) : TaskID(tid)
{
  init(p);
  FirstTime =true;
  KEcorr=0.0;
  OrbitalName = "J2OrbitalSoA";
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
  N_padded=getAlignedSize<valT>(N);
  NumGroups=p.groups();

  Uat.resize(N); 
  dUat.resize(N);
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
    else
      F[ia*NumGroups+ib]=j;
  }
  else
  {
    if(N==2)
    {
      // a very special case, 1 up + 1 down
      // uu/dd was prevented by the builder
      for(int ig=0; ig<NumGroups; ++ig)
        for(int jg=0; jg<NumGroups; ++jg)
          F[ig*NumGroups+jg]=j;
    }
    else
    {
      // generic case
      F[ia*NumGroups+ib]=j;
      F[ib*NumGroups+ia]=j;
    }
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
    RealType* restrict u, RealType* restrict du, RealType* restrict d2u, bool triangle)
{
  const int jelmax=triangle?iat:N;
  constexpr valT czero(0);
  std::fill_n(u,  jelmax,czero);
  std::fill_n(du, jelmax,czero);
  std::fill_n(d2u,jelmax,czero);

  const int igt=P.GroupID[iat]*NumGroups;
  for(int jg=0; jg<NumGroups; ++jg)
  {
    const FuncType& f2(*F[igt+jg]);
    int iStart = P.first(jg);
    int iEnd = std::min(jelmax,P.last(jg));
    f2.evaluateVGL(iat, iStart, iEnd, dist, u, du, d2u, DistCompressed.data(), DistIndice.data());
  }
  //u[iat]=czero;
  //du[iat]=czero;
  //d2u[iat]=czero;
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
    cur_Uat += f2.evaluateV(iat, iStart, iEnd, dist, DistCompressed.data());
  }

  return std::exp(Uat[iat]-cur_Uat);
}

template<typename FT>
inline void
J2OrbitalSoA<FT>::evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
{
  const DistanceTableData* d_table=P.DistTables[0];
  const auto dist=d_table->Temp_r.data();

  for(int ig=0; ig<NumGroups; ++ig)
  {
    const int igt=ig*NumGroups;
    valT sumU(0);
    for(int jg=0; jg<NumGroups; ++jg)
    {
      const FuncType& f2(*F[igt+jg]);
      int iStart = P.first(jg);
      int iEnd = P.last(jg);
      sumU += f2.evaluateV(-1, iStart, iEnd, dist, DistCompressed.data());
    }

    for(int i=P.first(ig); i<P.last(ig); ++i)
    {
      // remove self-interaction
      const valT Uself = F[igt+ig]->evaluate(dist[i]);
      ratios[i]=std::exp(Uat[i]+Uself-sumU);
    }
  }
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
  const auto& new_dr=d_table->Temp_dr;
  const auto& old_dr=d_table->Displacements[iat];
  constexpr valT lapfac=OHMMS_DIM-RealType(1);
  #pragma omp simd reduction(+:cur_d2Uat)
  for(int jat=0; jat<N; jat++)
  {
    const valT du   = cur_u[jat] - old_u[jat];
    const valT newl = cur_d2u[jat] + lapfac*cur_du[jat];
    const valT dl   = old_d2u[jat] + lapfac*old_du[jat] - newl;
    Uat[jat]   += du;
    d2Uat[jat] += dl;
    cur_d2Uat  -= newl;
  }
  posT cur_dUat;
  for(int idim=0; idim<OHMMS_DIM; ++idim)
  {
    const valT* restrict new_dX=new_dr.data(idim);
    const valT* restrict old_dX=old_dr.data(idim);
    const valT* restrict cur_du_pt=cur_du.data();
    const valT* restrict old_du_pt=old_du.data();
    valT* restrict save_g=dUat.data(idim);
    valT cur_g=cur_dUat[idim];
    #pragma omp simd reduction(+:cur_g) aligned(old_dX,new_dX,save_g,cur_du_pt,old_du_pt)
    for(int jat=0; jat<N; jat++)
    {
      const valT newg = cur_du_pt[jat] * new_dX[jat];
      const valT dg   = newg - old_du_pt[jat]*old_dX[jat];
      save_g[jat]  -= dg;
      cur_g += newg;
    }
    cur_dUat[idim] = cur_g;
  }
  Uat[iat]   = cur_Uat;
  dUat(iat)  = cur_dUat;
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
      computeU3(P,iat,d_table->Distances[iat],cur_u.data(),cur_du.data(),cur_d2u.data(),true);
      Uat[iat]=simd::accumulate_n(cur_u.data(),iat,valT());
      posT grad;
      valT lap(0);
      const valT* restrict    u = cur_u.data();
      const valT* restrict   du = cur_du.data();
      const valT* restrict  d2u = cur_d2u.data();
      const RowContainer& displ = d_table->Displacements[iat];
      constexpr valT lapfac=OHMMS_DIM-RealType(1);
      #pragma omp simd reduction(+:lap) aligned(du,d2u)
      for(int jat=0; jat<iat; ++jat)
        lap+=d2u[jat]+lapfac*du[jat];
      for(int idim=0; idim<OHMMS_DIM; ++idim)
      {
        const valT* restrict dX=displ.data(idim);
        valT s=valT();
        #pragma omp simd reduction(+:s) aligned(du,dX)
        for(int jat=0; jat<iat; ++jat) s+=du[jat]*dX[jat];
        grad[idim]=s;
      }
      dUat(iat)=grad;
      d2Uat[iat]=-lap;
      // add the contribution from the upper triangle
      #pragma omp simd aligned(u,du,d2u)
      for(int jat=0; jat<iat; jat++)
      {
        Uat[jat] += u[jat];
        d2Uat[jat] -= d2u[jat]+lapfac*du[jat];
      }
      for(int idim=0; idim<OHMMS_DIM; ++idim)
      {
        valT* restrict save_g=dUat.data(idim);
        const valT* restrict dX=displ.data(idim);
        #pragma omp simd aligned(save_g,du,dX)
        for(int jat=0; jat<iat; jat++)
          save_g[jat]-=du[jat]*dX[jat];
      }
    }
  }
}

template<typename FT>
typename J2OrbitalSoA<FT>::RealType
J2OrbitalSoA<FT>::evaluateLog(ParticleSet& P,
    ParticleSet::ParticleGradient_t& G,
    ParticleSet::ParticleLaplacian_t& L)
{
  evaluateGL(P,G,L,true);
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
