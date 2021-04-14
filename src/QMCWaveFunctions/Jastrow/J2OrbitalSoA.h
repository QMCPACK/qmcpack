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

#include <map>
#include <numeric>
#include "Configuration.h"
#if !defined(QMC_BUILD_SANDBOX_ONLY)
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Jastrow/DiffTwoBodyJastrowOrbital.h"
#endif
#include "Particle/DistanceTableData.h"
#include "LongRange/StructFact.h"
#include "CPU/SIMD/aligned_allocator.hpp"
#include "CPU/SIMD/algorithm.hpp"

namespace qmcplusplus
{
// helper class to activate KEcorr during optimizing Jastrow
template<typename RT, class FT>
class J2KECorrection
{
  size_t num_groups_;
  std::vector<size_t> num_elec_in_groups_;
  RT num_elecs_;
  RT vol;
  RT G0mag;
  const std::vector<FT*>& F_;
  bool SK_enabled;

public:
  J2KECorrection(const ParticleSet& targetPtcl, const std::vector<FT*>& F)
      : num_groups_(targetPtcl.groups()),
        num_elecs_(targetPtcl.getTotalNum()),
        vol(targetPtcl.Lattice.Volume),
        F_(F),
        SK_enabled(targetPtcl.SK != nullptr)
  {
    // compute num_elec_in_groups_
    num_elec_in_groups_.reserve(3);
    for (int i = 0; i < num_groups_; i++)
      num_elec_in_groups_.push_back(targetPtcl.last(i) - targetPtcl.first(i));

    if (SK_enabled)
      G0mag = std::sqrt(targetPtcl.SK->KLists.ksq[0]);
  }

  RT computeKEcorr()
  {
    if (!SK_enabled)
      return 0;

    const int numPoints = 1000;
    RT uk               = 0.0;
    RT a                = 1.0;

    for (int i = 0; i < num_groups_; i++)
    {
      int Ni = num_elec_in_groups_[i];
      for (int j = 0; j < num_groups_; j++)
      {
        int Nj = num_elec_in_groups_[j];
        if (F_[i * num_groups_ + j])
        {
          FT& ufunc = *(F_[i * num_groups_ + j]);
          RT radius = ufunc.cutoff_radius;
          RT k      = G0mag;
          RT dr     = radius / (RT)(numPoints - 1);
          for (int ir = 0; ir < numPoints; ir++)
          {
            RT r = dr * (RT)ir;
            RT u = ufunc.evaluate(r);
            uk += 0.5 * 4.0 * M_PI * r * std::sin(k * r) / k * u * dr * (RT)Nj / (RT)(Ni + Nj);
          }
        }
      }
    }
    for (int iter = 0; iter < 20; iter++)
      a = uk / (4.0 * M_PI * (1.0 / (G0mag * G0mag) - 1.0 / (G0mag * G0mag + 1.0 / a)));
    return 4.0 * M_PI * a / (4.0 * vol) * num_elecs_;
  }
};

/** @ingroup WaveFunctionComponent
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
class J2OrbitalSoA : public WaveFunctionComponent
{
public:
  ///alias FuncType
  using FuncType = FT;
  ///type of each component U, dU, d2U;
  using valT = typename FT::real_type;
  ///element position type
  using posT = TinyVector<valT, OHMMS_DIM>;
  ///use the same container
  using DistRow         = DistanceTableData::DistRow;
  using DisplRow        = DistanceTableData::DisplRow;
  using gContainer_type = VectorSoaContainer<valT, OHMMS_DIM>;

  // Ye: leaving this public is bad but currently used by unit tests.
  ///Container for \f$F[ig*NumGroups+jg]\f$.
  std::vector<FT*> F;

protected:
  ///number of particles
  size_t N;
  ///number of particles + padded
  size_t N_padded;
  ///number of groups of the target particleset
  size_t NumGroups;
  ///diff value
  RealType DiffVal;
  ///Correction
  RealType KEcorr;
  ///\f$Uat[i] = sum_(j) u_{i,j}\f$
  Vector<valT> Uat;
  ///\f$dUat[i] = sum_(j) du_{i,j}\f$
  gContainer_type dUat;
  ///\f$d2Uat[i] = sum_(j) d2u_{i,j}\f$
  Vector<valT> d2Uat;
  valT cur_Uat;
  aligned_vector<valT> cur_u, cur_du, cur_d2u;
  aligned_vector<valT> old_u, old_du, old_d2u;
  aligned_vector<valT> DistCompressed;
  aligned_vector<int> DistIndice;
  ///Uniquue J2 set for cleanup
  std::map<std::string, FT*> J2Unique;
  /// e-e table ID
  const int my_table_ID_;
  // helper for compute J2 Chiesa KE correction
  J2KECorrection<RealType, FT> j2_ke_corr_helper;

public:
  J2OrbitalSoA(const std::string& obj_name, ParticleSet& p, int tid);
  J2OrbitalSoA(const J2OrbitalSoA& rhs) = delete;
  ~J2OrbitalSoA();

  /* initialize storage */
  void init(ParticleSet& p);

  /** add functor for (ia,ib) pair */
  void addFunc(int ia, int ib, FT* j);

  /** check in an optimizable parameter
   * @param o a super set of optimizable variables
   */
  void checkInVariables(opt_variables_type& active)
  {
    myVars.clear();
    typename std::map<std::string, FT*>::iterator it(J2Unique.begin()), it_end(J2Unique.end());
    while (it != it_end)
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
    Optimizable = myVars.is_optimizable();
    typename std::map<std::string, FT*>::iterator it(J2Unique.begin()), it_end(J2Unique.end());
    while (it != it_end)
    {
      (*it).second->checkOutVariables(active);
      ++it;
    }
    if (dPsi)
      dPsi->checkOutVariables(active);
  }

  ///reset the value of all the unique Two-Body Jastrow functions
  void resetParameters(const opt_variables_type& active)
  {
    if (!Optimizable)
      return;
    typename std::map<std::string, FT*>::iterator it(J2Unique.begin()), it_end(J2Unique.end());
    while (it != it_end)
    {
      (*it).second->resetParameters(active);
      ++it;
    }
    if (dPsi)
      dPsi->resetParameters(active);
    for (int i = 0; i < myVars.size(); ++i)
    {
      int ii = myVars.Index[i];
      if (ii >= 0)
        myVars[i] = active[ii];
    }
  }


  void finalizeOptimization() { KEcorr = j2_ke_corr_helper.computeKEcorr(); }

  /** print the state, e.g., optimizables */
  void reportStatus(std::ostream& os)
  {
    typename std::map<std::string, FT*>::iterator it(J2Unique.begin()), it_end(J2Unique.end());
    while (it != it_end)
    {
      (*it).second->myVars.print(os);
      ++it;
    }
  }

  WaveFunctionComponentPtr makeClone(ParticleSet& tqp) const;

  LogValueType evaluateLog(const ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L);

  void evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi);

  /** recompute internal data assuming distance table is fully ready */
  void recompute(const ParticleSet& P);

  PsiValueType ratio(ParticleSet& P, int iat);
  void evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios)
  {
    for (int k = 0; k < ratios.size(); ++k)
      ratios[k] =
          std::exp(Uat[VP.refPtcl] - computeU(VP.refPS, VP.refPtcl, VP.getDistTable(my_table_ID_).getDistRow(k)));
  }
  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios);

  GradType evalGrad(ParticleSet& P, int iat);

  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);

  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false);
  inline void restore(int iat) {}

  /** compute G and L after the sweep
   */
  LogValueType evaluateGL(const ParticleSet& P,
                          ParticleSet::ParticleGradient_t& G,
                          ParticleSet::ParticleLaplacian_t& L,
                          bool fromscratch = false);

  inline void registerData(ParticleSet& P, WFBufferType& buf)
  {
    if (Bytes_in_WFBuffer == 0)
    {
      Bytes_in_WFBuffer = buf.current();
      buf.add(Uat.begin(), Uat.end());
      buf.add(dUat.data(), dUat.end());
      buf.add(d2Uat.begin(), d2Uat.end());
      Bytes_in_WFBuffer = buf.current() - Bytes_in_WFBuffer;
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
    dUat.attachReference(N, N_padded, buf.lendReference<valT>(N_padded * OHMMS_DIM));
    d2Uat.attachReference(buf.lendReference<valT>(N), N);
  }

  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false)
  {
    evaluateGL(P, P.G, P.L, false);
    buf.forward(Bytes_in_WFBuffer);
    return LogValue;
  }

  /*@{ internal compute engines*/
  inline valT computeU(const ParticleSet& P, int iat, const DistRow& dist)
  {
    valT curUat(0);
    const int igt = P.GroupID[iat] * NumGroups;
    for (int jg = 0; jg < NumGroups; ++jg)
    {
      const FuncType& f2(*F[igt + jg]);
      int iStart = P.first(jg);
      int iEnd   = P.last(jg);
      curUat += f2.evaluateV(iat, iStart, iEnd, dist.data(), DistCompressed.data());
    }
    return curUat;
  }

  inline void computeU3(const ParticleSet& P,
                        int iat,
                        const DistRow& dist,
                        RealType* restrict u,
                        RealType* restrict du,
                        RealType* restrict d2u,
                        bool triangle = false);

  /** compute gradient
   */
  inline posT accumulateG(const valT* restrict du, const DisplRow& displ) const
  {
    posT grad;
    for (int idim = 0; idim < OHMMS_DIM; ++idim)
    {
      const valT* restrict dX = displ.data(idim);
      valT s                  = valT();

#pragma omp simd reduction(+ : s) aligned(du, dX: QMC_SIMD_ALIGNMENT)
      for (int jat = 0; jat < N; ++jat)
        s += du[jat] * dX[jat];
      grad[idim] = s;
    }
    return grad;
  }
  /**@} */

  RealType ChiesaKEcorrection() { return KEcorr = j2_ke_corr_helper.computeKEcorr(); }

  RealType KECorrection() { return KEcorr; }
};

template<typename FT>
J2OrbitalSoA<FT>::J2OrbitalSoA(const std::string& obj_name, ParticleSet& p, int tid)
    : WaveFunctionComponent("J2OrbitalSoA", obj_name), my_table_ID_(p.addTable(p)), j2_ke_corr_helper(p, F)
{
  if (myName.empty())
    throw std::runtime_error("J2OrbitalSoA object name cannot be empty!");
  init(p);
  KEcorr = 0.0;
}

template<typename FT>
J2OrbitalSoA<FT>::~J2OrbitalSoA()
{
  auto it = J2Unique.begin();
  while (it != J2Unique.end())
  {
    delete ((*it).second);
    ++it;
  }
} //need to clean up J2Unique

template<typename FT>
void J2OrbitalSoA<FT>::init(ParticleSet& p)
{
  N         = p.getTotalNum();
  N_padded  = getAlignedSize<valT>(N);
  NumGroups = p.groups();

  Uat.resize(N);
  dUat.resize(N);
  d2Uat.resize(N);
  cur_u.resize(N);
  cur_du.resize(N);
  cur_d2u.resize(N);
  old_u.resize(N);
  old_du.resize(N);
  old_d2u.resize(N);
  F.resize(NumGroups * NumGroups, nullptr);
  DistCompressed.resize(N);
  DistIndice.resize(N);
}

template<typename FT>
void J2OrbitalSoA<FT>::addFunc(int ia, int ib, FT* j)
{
  if (ia == ib)
  {
    if (ia == 0) //first time, assign everything
    {
      int ij = 0;
      for (int ig = 0; ig < NumGroups; ++ig)
        for (int jg = 0; jg < NumGroups; ++jg, ++ij)
          if (F[ij] == nullptr)
            F[ij] = j;
    }
    else
      F[ia * NumGroups + ib] = j;
  }
  else
  {
    if (N == 2)
    {
      // a very special case, 1 up + 1 down
      // uu/dd was prevented by the builder
      for (int ig = 0; ig < NumGroups; ++ig)
        for (int jg = 0; jg < NumGroups; ++jg)
          F[ig * NumGroups + jg] = j;
    }
    else
    {
      // generic case
      F[ia * NumGroups + ib] = j;
      F[ib * NumGroups + ia] = j;
    }
  }
  std::stringstream aname;
  aname << ia << ib;
  J2Unique[aname.str()] = j;
}

template<typename FT>
WaveFunctionComponentPtr J2OrbitalSoA<FT>::makeClone(ParticleSet& tqp) const
{
  J2OrbitalSoA<FT>* j2copy = new J2OrbitalSoA<FT>(myName, tqp, -1);
  if (dPsi)
    j2copy->dPsi = dPsi->makeClone(tqp);
  std::map<const FT*, FT*> fcmap;
  for (int ig = 0; ig < NumGroups; ++ig)
    for (int jg = ig; jg < NumGroups; ++jg)
    {
      int ij = ig * NumGroups + jg;
      if (F[ij] == 0)
        continue;
      typename std::map<const FT*, FT*>::iterator fit = fcmap.find(F[ij]);
      if (fit == fcmap.end())
      {
        FT* fc = new FT(*F[ij]);
        j2copy->addFunc(ig, jg, fc);
        //if (dPsi) (j2copy->dPsi)->addFunc(aname.str(),ig,jg,fc);
        fcmap[F[ij]] = fc;
      }
    }
  j2copy->KEcorr      = KEcorr;
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
inline void J2OrbitalSoA<FT>::computeU3(const ParticleSet& P,
                                        int iat,
                                        const DistRow& dist,
                                        RealType* restrict u,
                                        RealType* restrict du,
                                        RealType* restrict d2u,
                                        bool triangle)
{
  const int jelmax = triangle ? iat : N;
  constexpr valT czero(0);
  std::fill_n(u, jelmax, czero);
  std::fill_n(du, jelmax, czero);
  std::fill_n(d2u, jelmax, czero);

  const int igt = P.GroupID[iat] * NumGroups;
  for (int jg = 0; jg < NumGroups; ++jg)
  {
    const FuncType& f2(*F[igt + jg]);
    int iStart = P.first(jg);
    int iEnd   = std::min(jelmax, P.last(jg));
    f2.evaluateVGL(iat, iStart, iEnd, dist.data(), u, du, d2u, DistCompressed.data(), DistIndice.data());
  }
  //u[iat]=czero;
  //du[iat]=czero;
  //d2u[iat]=czero;
}

template<typename FT>
typename J2OrbitalSoA<FT>::PsiValueType J2OrbitalSoA<FT>::ratio(ParticleSet& P, int iat)
{
  //only ratio, ready to compute it again
  UpdateMode = ORB_PBYP_RATIO;
  cur_Uat    = computeU(P, iat, P.getDistTable(my_table_ID_).getTempDists());
  return std::exp(static_cast<PsiValueType>(Uat[iat] - cur_Uat));
}

template<typename FT>
inline void J2OrbitalSoA<FT>::evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
{
  const auto& d_table = P.getDistTable(my_table_ID_);
  const auto& dist    = d_table.getTempDists();

  for (int ig = 0; ig < NumGroups; ++ig)
  {
    const int igt = ig * NumGroups;
    valT sumU(0);
    for (int jg = 0; jg < NumGroups; ++jg)
    {
      const FuncType& f2(*F[igt + jg]);
      int iStart = P.first(jg);
      int iEnd   = P.last(jg);
      sumU += f2.evaluateV(-1, iStart, iEnd, dist.data(), DistCompressed.data());
    }

    for (int i = P.first(ig); i < P.last(ig); ++i)
    {
      // remove self-interaction
      const valT Uself = F[igt + ig]->evaluate(dist[i]);
      ratios[i]        = std::exp(Uat[i] + Uself - sumU);
    }
  }
}

template<typename FT>
typename J2OrbitalSoA<FT>::GradType J2OrbitalSoA<FT>::evalGrad(ParticleSet& P, int iat)
{
  return GradType(dUat[iat]);
}

template<typename FT>
typename J2OrbitalSoA<FT>::PsiValueType J2OrbitalSoA<FT>::ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
{
  UpdateMode = ORB_PBYP_PARTIAL;

  computeU3(P, iat, P.getDistTable(my_table_ID_).getTempDists(), cur_u.data(), cur_du.data(), cur_d2u.data());
  cur_Uat = simd::accumulate_n(cur_u.data(), N, valT());
  DiffVal = Uat[iat] - cur_Uat;
  grad_iat += accumulateG(cur_du.data(), P.getDistTable(my_table_ID_).getTempDispls());
  return std::exp(static_cast<PsiValueType>(DiffVal));
}

template<typename FT>
void J2OrbitalSoA<FT>::acceptMove(ParticleSet& P, int iat, bool safe_to_delay)
{
  // get the old u, du, d2u
  const auto& d_table = P.getDistTable(my_table_ID_);
  computeU3(P, iat, d_table.getOldDists(), old_u.data(), old_du.data(), old_d2u.data());
  if (UpdateMode == ORB_PBYP_RATIO)
  { //ratio-only during the move; need to compute derivatives
    const auto& dist = d_table.getTempDists();
    computeU3(P, iat, dist, cur_u.data(), cur_du.data(), cur_d2u.data());
  }

  valT cur_d2Uat(0);
  const auto& new_dr    = d_table.getTempDispls();
  const auto& old_dr    = d_table.getOldDispls();
  constexpr valT lapfac = OHMMS_DIM - RealType(1);
#pragma omp simd reduction(+ : cur_d2Uat)
  for (int jat = 0; jat < N; jat++)
  {
    const valT du   = cur_u[jat] - old_u[jat];
    const valT newl = cur_d2u[jat] + lapfac * cur_du[jat];
    const valT dl   = old_d2u[jat] + lapfac * old_du[jat] - newl;
    Uat[jat] += du;
    d2Uat[jat] += dl;
    cur_d2Uat -= newl;
  }
  posT cur_dUat;
  for (int idim = 0; idim < OHMMS_DIM; ++idim)
  {
    const valT* restrict new_dX    = new_dr.data(idim);
    const valT* restrict old_dX    = old_dr.data(idim);
    const valT* restrict cur_du_pt = cur_du.data();
    const valT* restrict old_du_pt = old_du.data();
    valT* restrict save_g          = dUat.data(idim);
    valT cur_g                     = cur_dUat[idim];
#pragma omp simd reduction(+ : cur_g) aligned(old_dX, new_dX, save_g, cur_du_pt, old_du_pt: QMC_SIMD_ALIGNMENT)
    for (int jat = 0; jat < N; jat++)
    {
      const valT newg = cur_du_pt[jat] * new_dX[jat];
      const valT dg   = newg - old_du_pt[jat] * old_dX[jat];
      save_g[jat] -= dg;
      cur_g += newg;
    }
    cur_dUat[idim] = cur_g;
  }
  LogValue += Uat[iat] - cur_Uat;
  Uat[iat]   = cur_Uat;
  dUat(iat)  = cur_dUat;
  d2Uat[iat] = cur_d2Uat;
}

template<typename FT>
void J2OrbitalSoA<FT>::recompute(const ParticleSet& P)
{
  const auto& d_table = P.getDistTable(my_table_ID_);
  for (int ig = 0; ig < NumGroups; ++ig)
  {
    for (int iat = P.first(ig), last = P.last(ig); iat < last; ++iat)
    {
      computeU3(P, iat, d_table.getDistRow(iat), cur_u.data(), cur_du.data(), cur_d2u.data(), true);
      Uat[iat] = simd::accumulate_n(cur_u.data(), iat, valT());
      posT grad;
      valT lap(0);
      const valT* restrict u   = cur_u.data();
      const valT* restrict du  = cur_du.data();
      const valT* restrict d2u = cur_d2u.data();
      const auto& displ        = d_table.getDisplRow(iat);
      constexpr valT lapfac    = OHMMS_DIM - RealType(1);
#pragma omp simd reduction(+ : lap) aligned(du, d2u: QMC_SIMD_ALIGNMENT)
      for (int jat = 0; jat < iat; ++jat)
        lap += d2u[jat] + lapfac * du[jat];
      for (int idim = 0; idim < OHMMS_DIM; ++idim)
      {
        const valT* restrict dX = displ.data(idim);
        valT s                  = valT();
#pragma omp simd reduction(+ : s) aligned(du, dX: QMC_SIMD_ALIGNMENT)
        for (int jat = 0; jat < iat; ++jat)
          s += du[jat] * dX[jat];
        grad[idim] = s;
      }
      dUat(iat)  = grad;
      d2Uat[iat] = -lap;
// add the contribution from the upper triangle
#pragma omp simd aligned(u, du, d2u: QMC_SIMD_ALIGNMENT)
      for (int jat = 0; jat < iat; jat++)
      {
        Uat[jat] += u[jat];
        d2Uat[jat] -= d2u[jat] + lapfac * du[jat];
      }
      for (int idim = 0; idim < OHMMS_DIM; ++idim)
      {
        valT* restrict save_g   = dUat.data(idim);
        const valT* restrict dX = displ.data(idim);
#pragma omp simd aligned(save_g, du, dX: QMC_SIMD_ALIGNMENT)
        for (int jat = 0; jat < iat; jat++)
          save_g[jat] -= du[jat] * dX[jat];
      }
    }
  }
}

template<typename FT>
typename J2OrbitalSoA<FT>::LogValueType J2OrbitalSoA<FT>::evaluateLog(const ParticleSet& P,
                                                                      ParticleSet::ParticleGradient_t& G,
                                                                      ParticleSet::ParticleLaplacian_t& L)
{
  return evaluateGL(P, G, L, true);
}

template<typename FT>
WaveFunctionComponent::LogValueType J2OrbitalSoA<FT>::evaluateGL(const ParticleSet& P,
                                                                 ParticleSet::ParticleGradient_t& G,
                                                                 ParticleSet::ParticleLaplacian_t& L,
                                                                 bool fromscratch)
{
  if (fromscratch)
    recompute(P);
  LogValue = valT(0);
  for (int iat = 0; iat < N; ++iat)
  {
    LogValue += Uat[iat];
    G[iat] += dUat[iat];
    L[iat] += d2Uat[iat];
  }

  return LogValue = -LogValue * 0.5;
}

template<typename FT>
void J2OrbitalSoA<FT>::evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi)
{
  LogValue = 0.0;
  const DistanceTableData& d_ee(P.getDistTable(my_table_ID_));
  valT dudr, d2udr2;

  Tensor<valT, DIM> ident;
  grad_grad_psi = 0.0;
  ident.diagonal(1.0);

  for (int i = 1; i < N; ++i)
  {
    const auto& dist  = d_ee.getDistRow(i);
    const auto& displ = d_ee.getDisplRow(i);
    auto ig           = P.GroupID[i];
    const int igt     = ig * NumGroups;
    for (int j = 0; j < i; ++j)
    {
      auto r    = dist[j];
      auto rinv = 1.0 / r;
      auto dr   = displ[j];
      auto jg   = P.GroupID[j];
      auto uij  = F[igt + jg]->evaluate(r, dudr, d2udr2);
      LogValue -= uij;
      auto hess = rinv * rinv * outerProduct(dr, dr) * (d2udr2 - dudr * rinv) + ident * dudr * rinv;
      grad_grad_psi[i] -= hess;
      grad_grad_psi[j] -= hess;
    }
  }
}

} // namespace qmcplusplus
#endif
