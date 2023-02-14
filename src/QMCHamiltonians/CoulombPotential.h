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
#include "DistanceTable.h"
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
        myTableIndex(s.addTable(s, DTModes::NEED_FULL_TABLE_ON_HOST_AFTER_DONEPBYP)),
        is_AA(true),
        is_active(active),
        ComputeForces(computeForces)
  {
    setEnergyDomain(POTENTIAL);
    twoBodyQuantumDomain(s, s);
    nCenters = s.getTotalNum();
    prefix_   = "F_AA";

    if (!is_active) //precompute the value
    {
      if (!copy)
        s.update();
      value_ = evaluateAA(s.getDistTableAA(myTableIndex), s.Z.first_address());
      if (ComputeForces)
        evaluateAAForces(s.getDistTableAA(myTableIndex), s.Z.first_address());
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
    setEnergyDomain(POTENTIAL);
    twoBodyQuantumDomain(s, t);
    nCenters = s.getTotalNum();
  }

  std::string getClassName() const override { return "CoulombPotential"; }

#if !defined(REMOVE_TRACEMANAGER)
  void contributeParticleQuantities() override { request_.contribute_array(name_); }

  void checkoutParticleQuantities(TraceManager& tm) override
  {
    streaming_particles_ = request_.streaming_array(name_);
    if (streaming_particles_)
    {
      Va_sample = tm.checkout_real<1>(name_, Pa);
      if (!is_AA)
      {
        Vb_sample = tm.checkout_real<1>(name_, Pb);
      }
      else if (!is_active)
        evaluate_spAA(Pa.getDistTableAA(myTableIndex), Pa.Z.first_address());
    }
  }

  void deleteParticleQuantities() override
  {
    if (streaming_particles_)
    {
      delete Va_sample;
      if (!is_AA)
        delete Vb_sample;
    }
  }
#endif

  inline void addObservables(PropertySetType& plist, BufferType& collectables) override
  {
    addValue(plist);
    if (ComputeForces)
      addObservablesF(plist);
  }

  /** evaluate AA-type interactions */
  inline T evaluateAA(const DistanceTableAA& d, const ParticleScalar* restrict Z)
  {
    T res = 0.0;
#if !defined(REMOVE_TRACEMANAGER)
    if (streaming_particles_)
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
  inline void evaluateAAForces(const DistanceTableAA& d, const ParticleScalar* restrict Z)
  {
    forces_ = 0.0;
    for (size_t iat = 1; iat < nCenters; ++iat)
    {
      const auto& dist  = d.getDistRow(iat);
      const auto& displ = d.getDisplRow(iat);
      T q               = Z[iat];
      for (size_t j = 0; j < iat; ++j)
      {
        forces_[iat] += -q * Z[j] * displ[j] / (dist[j] * dist[j] * dist[j]);
        forces_[j] -= -q * Z[j] * displ[j] / (dist[j] * dist[j] * dist[j]);
      }
    }
  }


  /** JNKIM: Need to check the precision */
  inline T evaluateAB(const DistanceTableAB& d, const ParticleScalar* restrict Za, const ParticleScalar* restrict Zb)
  {
    constexpr T czero(0);
    T res = czero;
#if !defined(REMOVE_TRACEMANAGER)
    if (streaming_particles_)
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
  inline T evaluate_spAA(const DistanceTableAA& d, const ParticleScalar* restrict Z)
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
    auto sptmp           = streaming_particles_;
    streaming_particles_ = false;
    T Vnow               = res;
    T Vsum               = Va_samp.sum();
    T Vorig              = evaluateAA(d, Z);
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
    streaming_particles_ = sptmp;
#endif
    return res;
  }


  inline T evaluate_spAB(const DistanceTableAB& d, const ParticleScalar* restrict Za, const ParticleScalar* restrict Zb)
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
    auto sptmp           = streaming_particles_;
    streaming_particles_ = false;
    T Vnow               = res;
    T Vasum              = Va_samp.sum();
    T Vbsum              = Vb_samp.sum();
    T Vsum               = Vasum + Vbsum;
    T Vorig              = evaluateAB(d, Za, Zb);
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
    streaming_particles_ = sptmp;
#endif
    return res;
  }
#endif


  void resetTargetParticleSet(ParticleSet& P) override
  {
    //myTableIndex is the same
  }

  ~CoulombPotential() override {}

  void updateSource(ParticleSet& s) override
  {
    if (is_AA)
    {
      value_ = evaluateAA(s.getDistTableAA(myTableIndex), s.Z.first_address());
    }
  }

  inline Return_t evaluate(ParticleSet& P) override
  {
    if (is_active)
    {
      if (is_AA)
        value_ = evaluateAA(P.getDistTableAA(myTableIndex), P.Z.first_address());
      else
        value_ = evaluateAB(P.getDistTableAB(myTableIndex), Pa.Z.first_address(), P.Z.first_address());
    }
    return value_;
  }

  inline Return_t evaluateWithIonDerivs(ParticleSet& P,
                                        ParticleSet& ions,
                                        TrialWaveFunction& psi,
                                        ParticleSet::ParticlePos& hf_terms,
                                        ParticleSet::ParticlePos& pulay_terms) override
  {
    if (is_active)
      value_ = evaluate(P); // No forces for the active
    else
      hf_terms -= forces_;   // No Pulay here
    return value_;
  }

  bool put(xmlNodePtr cur) override { return true; }

  bool get(std::ostream& os) const override
  {
    if (myTableIndex)
      os << "CoulombAB source=" << Pa.getName() << std::endl;
    else
      os << "CoulombAA source/target " << Pa.getName() << std::endl;
    return true;
  }

  void setObservables(PropertySetType& plist) override
  {
    OperatorBase::setObservables(plist);
    if (ComputeForces)
      setObservablesF(plist);
  }

  void setParticlePropertyList(PropertySetType& plist, int offset) override
  {
    OperatorBase::setParticlePropertyList(plist, offset);
    if (ComputeForces)
      setParticleSetF(plist, offset);
  }


  std::unique_ptr<OperatorBase> makeClone(ParticleSet& qp, TrialWaveFunction& psi) override
  {
    if (is_AA)
    {
      if (is_active)
        return std::make_unique<CoulombPotential>(qp, true, ComputeForces);
      else
        // Ye Luo April 16th, 2015
        // avoid recomputing ion-ion DistanceTable when reusing ParticleSet
        return std::make_unique<CoulombPotential>(Pa, false, ComputeForces, true);
    }
    else
      return std::make_unique<CoulombPotential>(Pa, qp, true);
  }
};

} // namespace qmcplusplus
#endif
