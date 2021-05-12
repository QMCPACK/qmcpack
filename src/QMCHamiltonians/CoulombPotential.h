//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_COULOMBPOTENTIAL_H
#define QMCPLUSPLUS_COULOMBPOTENTIAL_H
#include "ParticleSet.h"
#include "WalkerSetRef.h"
#include "DistanceTableData.h"
#include "MCWalkerConfiguration.h"
#include "QMCHamiltonians/ForceBase.h"
#include "QMCHamiltonians/OperatorBase.h"
#include "QMCHamiltonians/OperatorBase.h"
#include <numeric>

namespace qmcplusplus
{
using WP = WalkerProperties::Indexes;

/** CoulombPotential
 * @tparam T type of the elementary data
 *
 * Hamiltonian operator for the Coulomb interaction for both AA and AB type for open systems.
 */
template<typename T>
struct CoulombPotential : public OperatorBase, public ForceBase
{
  ///source particle set
  ParticleSet& Pa;
  ///target particle set
  ParticleSet& Pb;
  ///distance table index
  const int myTableIndex;
  ///true if the table is AA
  const bool is_AA;
  ///true, if CoulombAA for quantum particleset
  bool is_active;
  ///number of centers
  int nCenters;
#if !defined(REMOVE_TRACEMANAGER)
  ///single particle trace samples
  Array<TraceReal, 1>* Va_sample;
  Array<TraceReal, 1>* Vb_sample;
#endif

  /// Flag for whether to compute forces or not
  bool ComputeForces;

  /** constructor for AA
   * @param s source particleset
   * @param active if true, new Value is computed whenver evaluate is used.
   * @param computeForces if true, computes forces between inactive species
   */
  inline CoulombPotential(ParticleSet& s, bool active, bool computeForces, bool copy = false)
      : ForceBase(s, s),
        Pa(s),
        Pb(s),
        myTableIndex(s.addTable(s)),
        is_AA(true),
        is_active(active),
        ComputeForces(computeForces)
  {
    set_energy_domain(potential);
    two_body_quantum_domain(s, s);
    nCenters = s.getTotalNum();
    prefix   = "F_AA";

    if (!is_active) //precompute the value
    {
      if (!copy)
        s.update();
      Value = evaluateAA(s.getDistTable(myTableIndex), s.Z.first_address());
      if (ComputeForces)
        evaluateAAForces(s.getDistTable(myTableIndex), s.Z.first_address());
    }
  }

  /** constructor for AB
   * @param s source particleset
   * @param t target particleset
   * @param active if true, new Value is computed whenver evaluate is used.
   * @param ComputeForces is not implemented for AB
   */
  inline CoulombPotential(ParticleSet& s, ParticleSet& t, bool active, bool copy = false)
      : ForceBase(s, t),
        Pa(s),
        Pb(t),
        myTableIndex(t.addTable(s)),
        is_AA(false),
        is_active(active),
        ComputeForces(false)
  {
    set_energy_domain(potential);
    two_body_quantum_domain(s, t);
    nCenters = s.getTotalNum();
  }

#if !defined(REMOVE_TRACEMANAGER)
  virtual void contribute_particle_quantities() { request.contribute_array(myName); }

  virtual void checkout_particle_quantities(TraceManager& tm)
  {
    streaming_particles = request.streaming_array(myName);
    if (streaming_particles)
    {
      Va_sample = tm.checkout_real<1>(myName, Pa);
      if (!is_AA)
      {
        Vb_sample = tm.checkout_real<1>(myName, Pb);
      }
      else if (!is_active)
        evaluate_spAA(Pa.getDistTable(myTableIndex), Pa.Z.first_address());
    }
  }

  virtual void delete_particle_quantities()
  {
    if (streaming_particles)
    {
      delete Va_sample;
      if (!is_AA)
        delete Vb_sample;
    }
  }
#endif

  inline void addObservables(PropertySetType& plist, BufferType& collectables)
  {
    addValue(plist);
    if (ComputeForces)
      addObservablesF(plist);
  }

  /** evaluate AA-type interactions */
  inline T evaluateAA(const DistanceTableData& d, const ParticleScalar_t* restrict Z)
  {
    T res = 0.0;
#if !defined(REMOVE_TRACEMANAGER)
    if (streaming_particles)
      res = evaluate_spAA(d, Z);
    else
#endif
      for (size_t iat = 1; iat < nCenters; ++iat)
      {
        const auto& dist = d.getDistRow(iat);
        T q              = Z[iat];
        for (size_t j = 0; j < iat; ++j)
          res += q * Z[j] / dist[j];
      }
    return res;
  }


  /** evaluate AA-type forces */
  inline void evaluateAAForces(const DistanceTableData& d, const ParticleScalar_t* restrict Z)
  {
    forces = 0.0;
    for (size_t iat = 1; iat < nCenters; ++iat)
    {
      const auto& dist  = d.getDistRow(iat);
      const auto& displ = d.getDisplRow(iat);
      T q               = Z[iat];
      for (size_t j = 0; j < iat; ++j)
      {
        forces[iat] += -q * Z[j] * displ[j] / (dist[j] * dist[j] * dist[j]);
        forces[j] -= -q * Z[j] * displ[j] / (dist[j] * dist[j] * dist[j]);
      }
    }
  }


  /** JNKIM: Need to check the precision */
  inline T evaluateAB(const DistanceTableData& d,
                      const ParticleScalar_t* restrict Za,
                      const ParticleScalar_t* restrict Zb)
  {
    constexpr T czero(0);
    T res = czero;
#if !defined(REMOVE_TRACEMANAGER)
    if (streaming_particles)
      res = evaluate_spAB(d, Za, Zb);
    else
#endif
    {
      const size_t nTargets = d.targets();
      for (size_t b = 0; b < nTargets; ++b)
      {
        const auto& dist = d.getDistRow(b);
        T e              = czero;
        for (size_t a = 0; a < nCenters; ++a)
          e += Za[a] / dist[a];
        res += e * Zb[b];
      }
    }
    return res;
  }


#if !defined(REMOVE_TRACEMANAGER)
  /** evaluate AA-type interactions */
  inline T evaluate_spAA(const DistanceTableData& d, const ParticleScalar_t* restrict Z)
  {
    T res = 0.0;
    T pairpot;
    Array<RealType, 1>& Va_samp = *Va_sample;
    Va_samp                     = 0.0;
    for (size_t iat = 1; iat < nCenters; ++iat)
    {
      const auto& dist = d.getDistRow(iat);
      T q              = Z[iat];
      for (size_t j = 0; j < iat; ++j)
      {
        pairpot = 0.5 * q * Z[j] / dist[j];
        Va_samp(iat) += pairpot;
        Va_samp(j) += pairpot;
        res += pairpot;
      }
    }
    res *= 2.0;
#if defined(TRACE_CHECK)
    auto sptmp          = streaming_particles;
    streaming_particles = false;
    T Vnow              = res;
    T Vsum              = Va_samp.sum();
    T Vorig             = evaluateAA(d, Z);
    if (std::abs(Vorig - Vnow) > TraceManager::trace_tol)
    {
      app_log() << "versiontest: CoulombPotential::evaluateAA()" << std::endl;
      app_log() << "versiontest:   orig:" << Vorig << std::endl;
      app_log() << "versiontest:    mod:" << Vnow << std::endl;
      APP_ABORT("Trace check failed");
    }
    if (std::abs(Vsum - Vnow) > TraceManager::trace_tol)
    {
      app_log() << "accumtest: CoulombPotential::evaluateAA()" << std::endl;
      app_log() << "accumtest:   tot:" << Vnow << std::endl;
      app_log() << "accumtest:   sum:" << Vsum << std::endl;
      APP_ABORT("Trace check failed");
    }
    streaming_particles = sptmp;
#endif
    return res;
  }


  inline T evaluate_spAB(const DistanceTableData& d,
                         const ParticleScalar_t* restrict Za,
                         const ParticleScalar_t* restrict Zb)
  {
    T res = 0.0;
    T pairpot;
    Array<RealType, 1>& Va_samp = *Va_sample;
    Array<RealType, 1>& Vb_samp = *Vb_sample;
    Va_samp                     = 0.0;
    Vb_samp                     = 0.0;
    const size_t nTargets       = d.targets();
    for (size_t b = 0; b < nTargets; ++b)
    {
      const auto& dist = d.getDistRow(b);
      T z              = 0.5 * Zb[b];
      for (size_t a = 0; a < nCenters; ++a)
      {
        pairpot = z * Za[a] / dist[a];
        Va_samp(a) += pairpot;
        Vb_samp(b) += pairpot;
        res += pairpot;
      }
    }
    res *= 2.0;

#if defined(TRACE_CHECK)
    auto sptmp          = streaming_particles;
    streaming_particles = false;
    T Vnow              = res;
    T Vasum             = Va_samp.sum();
    T Vbsum             = Vb_samp.sum();
    T Vsum              = Vasum + Vbsum;
    T Vorig             = evaluateAB(d, Za, Zb);
    if (std::abs(Vorig - Vnow) > TraceManager::trace_tol)
    {
      app_log() << "versiontest: CoulombPotential::evaluateAB()" << std::endl;
      app_log() << "versiontest:   orig:" << Vorig << std::endl;
      app_log() << "versiontest:    mod:" << Vnow << std::endl;
      APP_ABORT("Trace check failed");
    }
    if (std::abs(Vsum - Vnow) > TraceManager::trace_tol)
    {
      app_log() << "accumtest: CoulombPotential::evaluateAB()" << std::endl;
      app_log() << "accumtest:   tot:" << Vnow << std::endl;
      app_log() << "accumtest:   sum:" << Vsum << std::endl;
      APP_ABORT("Trace check failed");
    }
    if (std::abs(Vasum - Vbsum) > TraceManager::trace_tol)
    {
      app_log() << "sharetest: CoulombPotential::evaluateAB()" << std::endl;
      app_log() << "sharetest:   a share:" << Vasum << std::endl;
      app_log() << "sharetest:   b share:" << Vbsum << std::endl;
      APP_ABORT("Trace check failed");
    }
    streaming_particles = sptmp;
#endif
    return res;
  }
#endif


  void resetTargetParticleSet(ParticleSet& P)
  {
    //myTableIndex is the same
  }

  ~CoulombPotential() {}

  void update_source(ParticleSet& s)
  {
    if (is_AA)
    {
      Value = evaluateAA(s.getDistTable(myTableIndex), s.Z.first_address());
    }
  }

  inline Return_t evaluate(ParticleSet& P)
  {
    if (is_active)
    {
      if (is_AA)
        Value = evaluateAA(P.getDistTable(myTableIndex), P.Z.first_address());
      else
        Value = evaluateAB(P.getDistTable(myTableIndex), Pa.Z.first_address(), P.Z.first_address());
    }
    return Value;
  }

  inline Return_t evaluateWithIonDerivs(ParticleSet& P,
                                        ParticleSet& ions,
                                        TrialWaveFunction& psi,
                                        ParticleSet::ParticlePos_t& hf_terms,
                                        ParticleSet::ParticlePos_t& pulay_terms)
  {
    if (is_active)
      Value = evaluate(P); // No forces for the active
    else
      hf_terms -= forces; // No Pulay here
    return Value;
  }

  bool put(xmlNodePtr cur) { return true; }

  bool get(std::ostream& os) const
  {
    if (myTableIndex)
      os << "CoulombAB source=" << Pa.getName() << std::endl;
    else
      os << "CoulombAA source/target " << Pa.getName() << std::endl;
    return true;
  }

  void setObservables(PropertySetType& plist)
  {
    OperatorBase::setObservables(plist);
    if (ComputeForces)
      setObservablesF(plist);
  }

  void setParticlePropertyList(PropertySetType& plist, int offset)
  {
    OperatorBase::setParticlePropertyList(plist, offset);
    if (ComputeForces)
      setParticleSetF(plist, offset);
  }

  OperatorBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    if (is_AA)
    {
      if (is_active)
        return new CoulombPotential(qp, true, ComputeForces);
      else
        // Ye Luo April 16th, 2015
        // avoid recomputing ion-ion DistanceTable when reusing ParticleSet
        return new CoulombPotential(Pa, false, ComputeForces, true);
    }
    else
      return new CoulombPotential(Pa, qp, true);
  }

  void addEnergy(MCWalkerConfiguration& W, std::vector<RealType>& LocalEnergy)
  {
    auto& walkers = W.WalkerList;
    if (is_active)
    {
      for (int iw = 0; iw < W.getActiveWalkers(); iw++)
      {
        W.loadWalker(*walkers[iw], false);
        W.update();
        Value = evaluate(W);
        walkers[iw]->getPropertyBase()[WP::NUMPROPERTIES + myIndex] = Value;
        LocalEnergy[iw] += Value;
      }
    }
    else
      // assuminig the same results for all the walkers when the set is not active
      for (int iw = 0; iw < LocalEnergy.size(); iw++)
      {
        walkers[iw]->getPropertyBase()[WP::NUMPROPERTIES + myIndex] = Value;
        LocalEnergy[iw] += Value;
      }
  }
};

} // namespace qmcplusplus
#endif
