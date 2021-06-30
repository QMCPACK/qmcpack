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
#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "QMCWaveFunctions/Jastrow/DiffOneBodyJastrowOrbital.h"
#include "Utilities/qmc_common.h"
#include "CPU/SIMD/aligned_allocator.hpp"
#include "CPU/SIMD/algorithm.hpp"
#include <map>
#include <numeric>

namespace qmcplusplus
{
/** @ingroup WaveFunctionComponent
 *  @brief Specialization for one-body Jastrow function using multiple functors
 */
template<class FT>
struct J1OrbitalSoA : public WaveFunctionComponent
{
  ///alias FuncType
  using FuncType = FT;
  ///type of each component U, dU, d2U;
  using valT = typename FT::real_type;
  ///element position type
  using posT = TinyVector<valT, OHMMS_DIM>;
  ///use the same container
  using DistRow  = DistanceTableData::DistRow;
  using DisplRow = DistanceTableData::DisplRow;
  ///table index
  const int myTableID;
  ///number of ions
  int Nions;
  ///number of electrons
  int Nelec;
  ///number of groups
  int NumGroups;
  ///reference to the sources (ions)
  const ParticleSet& Ions;

  valT curAt;
  valT curLap;
  posT curGrad;

  ///\f$Vat[i] = sum_(j) u_{i,j}\f$
  Vector<valT> Vat;
  aligned_vector<valT> U, dU, d2U, d3U;
  aligned_vector<valT> DistCompressed;
  aligned_vector<int> DistIndice;
  Vector<posT> Grad;
  Vector<valT> Lap;
  ///Container for \f$F[ig*NumGroups+jg]\f$
  std::vector<FT*> F;

  J1OrbitalSoA(const std::string& obj_name, const ParticleSet& ions, ParticleSet& els)
      : WaveFunctionComponent("J1OrbitalSoA", obj_name), myTableID(els.addTable(ions)), Ions(ions)
  {
    if (myName.empty())
      throw std::runtime_error("J1OrbitalSoA object name cannot be empty!");
    initialize(els);
  }

  J1OrbitalSoA(const J1OrbitalSoA& rhs) = delete;

  ~J1OrbitalSoA()
  {
    for (int i = 0; i < F.size(); ++i)
      if (F[i] != nullptr)
        delete F[i];
  }

  /* initialize storage */
  void initialize(const ParticleSet& els)
  {
    Nions     = Ions.getTotalNum();
    NumGroups = Ions.getSpeciesSet().getTotalNum();
    F.resize(std::max(NumGroups, 4), nullptr);
    if (NumGroups > 1 && !Ions.IsGrouped)
    {
      NumGroups = 0;
    }
    Nelec = els.getTotalNum();
    Vat.resize(Nelec);
    Grad.resize(Nelec);
    Lap.resize(Nelec);

    U.resize(Nions);
    dU.resize(Nions);
    d2U.resize(Nions);
    d3U.resize(Nions);
    DistCompressed.resize(Nions);
    DistIndice.resize(Nions);
  }

  void addFunc(int source_type, FT* afunc, int target_type = -1)
  {
    if (F[source_type] != nullptr)
      delete F[source_type];
    F[source_type] = afunc;
  }

  void recompute(const ParticleSet& P)
  {
    const DistanceTableData& d_ie(P.getDistTable(myTableID));
    for (int iat = 0; iat < Nelec; ++iat)
    {
      computeU3(P, iat, d_ie.getDistRow(iat));
      Vat[iat] = simd::accumulate_n(U.data(), Nions, valT());
      Lap[iat] = accumulateGL(dU.data(), d2U.data(), d_ie.getDisplRow(iat), Grad[iat]);
    }
  }

  LogValueType evaluateLog(const ParticleSet& P, ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
  {
    return evaluateGL(P, G, L, true);
  }

  void evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi)
  {
    const DistanceTableData& d_ie(P.getDistTable(myTableID));
    valT dudr, d2udr2;

    Tensor<valT, DIM> ident;
    grad_grad_psi = 0.0;
    ident.diagonal(1.0);

    for (int iel = 0; iel < Nelec; ++iel)
    {
      const auto& dist  = d_ie.getDistRow(iel);
      const auto& displ = d_ie.getDisplRow(iel);
      for (int iat = 0; iat < Nions; iat++)
      {
        int gid    = Ions.GroupID[iat];
        auto* func = F[gid];
        if (func != nullptr)
        {
          RealType r    = dist[iat];
          RealType rinv = 1.0 / r;
          PosType dr    = displ[iat];
          func->evaluate(r, dudr, d2udr2);
          grad_grad_psi[iel] -= rinv * rinv * outerProduct(dr, dr) * (d2udr2 - dudr * rinv) + ident * dudr * rinv;
        }
      }
    }
  }

  PsiValueType ratio(ParticleSet& P, int iat)
  {
    UpdateMode = ORB_PBYP_RATIO;
    curAt      = computeU(P.getDistTable(myTableID).getTempDists());
    return std::exp(static_cast<PsiValueType>(Vat[iat] - curAt));
  }

  inline void evaluateRatios(const VirtualParticleSet& VP, std::vector<ValueType>& ratios)
  {
    for (int k = 0; k < ratios.size(); ++k)
      ratios[k] = std::exp(Vat[VP.refPtcl] - computeU(VP.getDistTable(myTableID).getDistRow(k)));
  }

  inline valT computeU(const DistRow& dist)
  {
    valT curVat(0);
    if (NumGroups > 0)
    {
      for (int jg = 0; jg < NumGroups; ++jg)
      {
        if (F[jg] != nullptr)
          curVat += F[jg]->evaluateV(-1, Ions.first(jg), Ions.last(jg), dist.data(), DistCompressed.data());
      }
    }
    else
    {
      for (int c = 0; c < Nions; ++c)
      {
        int gid = Ions.GroupID[c];
        if (F[gid] != nullptr)
          curVat += F[gid]->evaluate(dist[c]);
      }
    }
    return curVat;
  }

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
  {
    const auto& dist = P.getDistTable(myTableID).getTempDists();
    curAt            = valT(0);
    if (NumGroups > 0)
    {
      for (int jg = 0; jg < NumGroups; ++jg)
      {
        if (F[jg] != nullptr)
          curAt += F[jg]->evaluateV(-1, Ions.first(jg), Ions.last(jg), dist.data(), DistCompressed.data());
      }
    }
    else
    {
      for (int c = 0; c < Nions; ++c)
      {
        int gid = Ions.GroupID[c];
        if (F[gid] != nullptr)
          curAt += F[gid]->evaluate(dist[c]);
      }
    }

    for (int i = 0; i < Nelec; ++i)
      ratios[i] = std::exp(Vat[i] - curAt);
  }

  inline LogValueType evaluateGL(const ParticleSet& P,
                                 ParticleSet::ParticleGradient_t& G,
                                 ParticleSet::ParticleLaplacian_t& L,
                                 bool fromscratch = false)
  {
    if (fromscratch)
      recompute(P);

    for (size_t iat = 0; iat < Nelec; ++iat)
      G[iat] += Grad[iat];
    for (size_t iat = 0; iat < Nelec; ++iat)
      L[iat] -= Lap[iat];
    return LogValue = -simd::accumulate_n(Vat.data(), Nelec, valT());
  }

  /** compute gradient and lap
   * @return lap
   */
  inline valT accumulateGL(const valT* restrict du, const valT* restrict d2u, const DisplRow& displ, posT& grad) const
  {
    valT lap(0);
    constexpr valT lapfac = OHMMS_DIM - RealType(1);
    //#pragma omp simd reduction(+:lap)
    for (int jat = 0; jat < Nions; ++jat)
      lap += d2u[jat] + lapfac * du[jat];
    for (int idim = 0; idim < OHMMS_DIM; ++idim)
    {
      const valT* restrict dX = displ.data(idim);
      valT s                  = valT();
      //#pragma omp simd reduction(+:s)
      for (int jat = 0; jat < Nions; ++jat)
        s += du[jat] * dX[jat];
      grad[idim] = s;
    }
    return lap;
  }

  /** compute U, dU and d2U 
   * @param P quantum particleset
   * @param iat the moving particle
   * @param dist starting address of the distances of the ions wrt the iat-th particle
   */
  inline void computeU3(const ParticleSet& P, int iat, const DistRow& dist)
  {
    if (NumGroups > 0)
    { //ions are grouped
      constexpr valT czero(0);
      std::fill_n(U.data(), Nions, czero);
      std::fill_n(dU.data(), Nions, czero);
      std::fill_n(d2U.data(), Nions, czero);

      for (int jg = 0; jg < NumGroups; ++jg)
      {
        if (F[jg] == nullptr)
          continue;
        F[jg]->evaluateVGL(-1, Ions.first(jg), Ions.last(jg), dist.data(), U.data(), dU.data(), d2U.data(),
                           DistCompressed.data(), DistIndice.data());
      }
    }
    else
    {
      for (int c = 0; c < Nions; ++c)
      {
        int gid = Ions.GroupID[c];
        if (F[gid] != nullptr)
        {
          U[c] = F[gid]->evaluate(dist[c], dU[c], d2U[c]);
          dU[c] /= dist[c];
        }
      }
    }
  }

  /** compute the gradient during particle-by-particle update
   * @param P quantum particleset
   * @param iat particle index
   */
  GradType evalGrad(ParticleSet& P, int iat) { return GradType(Grad[iat]); }

  /** compute the gradient during particle-by-particle update
   * @param P quantum particleset
   * @param iat particle index
   *
   * Using getTempDists(). curAt, curGrad and curLap are computed.
   */
  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
  {
    UpdateMode = ORB_PBYP_PARTIAL;

    computeU3(P, iat, P.getDistTable(myTableID).getTempDists());
    curLap = accumulateGL(dU.data(), d2U.data(), P.getDistTable(myTableID).getTempDispls(), curGrad);
    curAt  = simd::accumulate_n(U.data(), Nions, valT());
    grad_iat += curGrad;
    return std::exp(static_cast<PsiValueType>(Vat[iat] - curAt));
  }

  /** Rejected move. Nothing to do */
  inline void restore(int iat) {}

  /** Accpted move. Update Vat[iat],Grad[iat] and Lap[iat] */
  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false)
  {
    if (UpdateMode == ORB_PBYP_RATIO)
    {
      computeU3(P, iat, P.getDistTable(myTableID).getTempDists());
      curLap = accumulateGL(dU.data(), d2U.data(), P.getDistTable(myTableID).getTempDispls(), curGrad);
    }

    LogValue += Vat[iat] - curAt;
    Vat[iat]  = curAt;
    Grad[iat] = curGrad;
    Lap[iat]  = curLap;
  }


  inline void registerData(ParticleSet& P, WFBufferType& buf)
  {
    if (Bytes_in_WFBuffer == 0)
    {
      Bytes_in_WFBuffer = buf.current();
      buf.add(Vat.begin(), Vat.end());
      buf.add(Grad.begin(), Grad.end());
      buf.add(Lap.begin(), Lap.end());
      Bytes_in_WFBuffer = buf.current() - Bytes_in_WFBuffer;
      // free local space
      Vat.free();
      Grad.free();
      Lap.free();
    }
    else
    {
      buf.forward(Bytes_in_WFBuffer);
    }
  }

  inline LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false)
  {
    evaluateGL(P, P.G, P.L, false);
    buf.forward(Bytes_in_WFBuffer);
    return LogValue;
  }

  inline void copyFromBuffer(ParticleSet& P, WFBufferType& buf)
  {
    Vat.attachReference(buf.lendReference<valT>(Nelec), Nelec);
    Grad.attachReference(buf.lendReference<posT>(Nelec), Nelec);
    Lap.attachReference(buf.lendReference<valT>(Nelec), Nelec);
  }

  WaveFunctionComponentPtr makeClone(ParticleSet& tqp) const
  {
    J1OrbitalSoA<FT>* j1copy = new J1OrbitalSoA<FT>(myName, Ions, tqp);
    j1copy->Optimizable      = Optimizable;
    for (size_t i = 0, n = F.size(); i < n; ++i)
    {
      if (F[i] != nullptr)
        j1copy->addFunc(i, new FT(*F[i]));
    }
    if (dPsi)
    {
      j1copy->dPsi = dPsi->makeClone(tqp);
    }
    return j1copy;
  }

  /**@{ WaveFunctionComponent virtual functions that are not essential for the development */
  void reportStatus(std::ostream& os)
  {
    for (size_t i = 0, n = F.size(); i < n; ++i)
    {
      if (F[i] != nullptr)
        F[i]->myVars.print(os);
    }
  }

  void checkInVariables(opt_variables_type& active)
  {
    myVars.clear();
    for (size_t i = 0, n = F.size(); i < n; ++i)
    {
      if (F[i] != nullptr)
      {
        F[i]->checkInVariables(active);
        F[i]->checkInVariables(myVars);
      }
    }
  }
  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
    Optimizable = myVars.is_optimizable();
    for (size_t i = 0, n = F.size(); i < n; ++i)
      if (F[i] != nullptr)
        F[i]->checkOutVariables(active);
    if (dPsi)
      dPsi->checkOutVariables(active);
  }

  void resetParameters(const opt_variables_type& active)
  {
    if (!Optimizable)
      return;
    for (size_t i = 0, n = F.size(); i < n; ++i)
      if (F[i] != nullptr)
        F[i]->resetParameters(active);

    for (int i = 0; i < myVars.size(); ++i)
    {
      int ii = myVars.Index[i];
      if (ii >= 0)
        myVars[i] = active[ii];
    }
    if (dPsi)
      dPsi->resetParameters(active);
  }
  /**@} */

  inline GradType evalGradSource(ParticleSet& P, ParticleSet& source, int isrc)
  {
    GradType g_return(0.0);
    const DistanceTableData& d_ie(P.getDistTable(myTableID));
    for (int iat = 0; iat < Nelec; ++iat)
    {
      const auto& dist  = d_ie.getDistRow(iat);
      const auto& displ = d_ie.getDisplRow(iat);
      int gid           = source.GroupID[isrc];
      RealType r        = dist[isrc];
      RealType rinv     = 1.0 / r;
      PosType dr        = displ[isrc];

      if (F[gid] != nullptr)
      {
        U[isrc] = F[gid]->evaluate(dist[isrc], dU[isrc], d2U[isrc], d3U[isrc]);
        g_return -= dU[isrc] * rinv * dr;
      }
    }
    return g_return;
  }

  inline GradType evalGradSource(ParticleSet& P,
                                 ParticleSet& source,
                                 int isrc,
                                 TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM>& grad_grad,
                                 TinyVector<ParticleSet::ParticleLaplacian_t, OHMMS_DIM>& lapl_grad)
  {
    GradType g_return(0.0);
    const DistanceTableData& d_ie(P.getDistTable(myTableID));
    for (int iat = 0; iat < Nelec; ++iat)
    {
      const auto& dist  = d_ie.getDistRow(iat);
      const auto& displ = d_ie.getDisplRow(iat);
      int gid           = source.GroupID[isrc];
      RealType r        = dist[isrc];
      RealType rinv     = 1.0 / r;
      PosType dr        = displ[isrc];

      if (F[gid] != nullptr)
      {
        U[isrc] = F[gid]->evaluate(dist[isrc], dU[isrc], d2U[isrc], d3U[isrc]);
      }
      else
      {
        APP_ABORT("J1OrbitalSoa::evaluateGradSource:  F[gid]==nullptr")
      }

      g_return -= dU[isrc] * rinv * dr;

      //The following terms depend only on the radial component r.  Thus,
      //we compute them and mix with position vectors to acquire the full
      //cartesian vector objects.
      valT grad_component = (d2U[isrc] - dU[isrc] * rinv);
      valT lapl_component = d3U[isrc] + 2 * rinv * grad_component;

      for (int idim = 0; idim < OHMMS_DIM; idim++)
      {
        grad_grad[idim][iat] += dr[idim] * dr * rinv * rinv * grad_component;
        grad_grad[idim][iat][idim] += rinv * dU[isrc];

        lapl_grad[idim][iat] -= lapl_component * rinv * dr[idim];
      }
    }
    return g_return;
  }
};


} // namespace qmcplusplus
#endif
