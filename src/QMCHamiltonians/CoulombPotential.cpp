//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "CoulombPotential.h"
#include <numeric>
#include <Resource.h>

namespace qmcplusplus
{
using WP       = WalkerProperties::Indexes;
using Return_t = CoulombPotential::Return_t;

struct CoulombPotential::CoulombPotentialMultiWalkerResource : public Resource
{
  CoulombPotentialMultiWalkerResource() : Resource("CoulombPotential") {}

  std::unique_ptr<Resource> makeClone() const override
  {
    return std::make_unique<CoulombPotentialMultiWalkerResource>(*this);
  }

  /// a crowds worth of per particle local ecp potential values
  Vector<RealType> va_samples;
  Vector<RealType> vb_samples;
};

/** Hamiltonian operator for simple particle particle coulomb interaction.
 *
 *  This is more of a problem than you might thing because A can be either electron or ion for AA.
 *  A is ion and B is elec for AB interactions
 *
 *  Because ion (static particles) and elec (dynamic or target particles) aren't on equal footing.
 *  ion's are almost entirely state, elec are passed as arguements.
 *
 *  So this class tries to generically support any particle particle coulomb interaction but consistency isn't
 *  really possible.
 */
CoulombPotential::CoulombPotential(ParticleSet& s, bool active, bool computeForces, bool copy)
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
  prefix_  = "F_AA";

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
CoulombPotential::CoulombPotential(ParticleSet& s, ParticleSet& t, bool active, bool copy)
    : ForceBase(s, t), Pa(s), Pb(t), myTableIndex(t.addTable(s)), is_AA(false), is_active(active), ComputeForces(false)
{
  setEnergyDomain(POTENTIAL);
  twoBodyQuantumDomain(s, t);
  nCenters = s.getTotalNum();
}

std::string CoulombPotential::getClassName() const { return "CoulombPotential"; }

#if !defined(REMOVE_TRACEMANAGER)
void CoulombPotential::contributeParticleQuantities() { request_.contribute_array(name_); }

void CoulombPotential::checkoutParticleQuantities(TraceManager& tm)
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

void CoulombPotential::deleteParticleQuantities()
{
  if (streaming_particles_)
  {
    delete Va_sample;
    if (!is_AA)
      delete Vb_sample;
  }
}
#endif

void CoulombPotential::addObservables(PropertySetType& plist, BufferType& collectables)
{
  addValue(plist);
  if (ComputeForces)
    addObservablesF(plist);
}

/** evaluate AA-type interactions */
Return_t CoulombPotential::evaluateAA(const DistanceTableAA& d, const ParticleScalar* restrict Z)
{
  Return_t res = 0.0;
#if !defined(REMOVE_TRACEMANAGER)
  if (streaming_particles_)
    res = evaluate_spAA(d, Z);
  else
#endif
    for (size_t iat = 1; iat < nCenters; ++iat)
    {
      const auto& dist = d.getDistRow(iat);
      Return_t q       = Z[iat];
      for (size_t j = 0; j < iat; ++j)
        res += q * Z[j] / dist[j];
    }
  return res;
}


/** evaluate AA-type forces */
void CoulombPotential::evaluateAAForces(const DistanceTableAA& d, const ParticleScalar* restrict Z)
{
  forces_ = 0.0;
  for (size_t iat = 1; iat < nCenters; ++iat)
  {
    const auto& dist  = d.getDistRow(iat);
    const auto& displ = d.getDisplRow(iat);
    Return_t q        = Z[iat];
    for (size_t j = 0; j < iat; ++j)
    {
      forces_[iat] += -q * Z[j] * displ[j] / (dist[j] * dist[j] * dist[j]);
      forces_[j] -= -q * Z[j] * displ[j] / (dist[j] * dist[j] * dist[j]);
    }
  }
}


/** JNKIM: Need to check the precision */
Return_t CoulombPotential::evaluateAB(const DistanceTableAB& d,
                                      const ParticleScalar* restrict Za,
                                      const ParticleScalar* restrict Zb)
{
  constexpr Return_t czero(0);
  Return_t res = czero;
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
      Return_t e       = czero;
      for (size_t a = 0; a < nCenters; ++a)
        e += Za[a] / dist[a];
      res += e * Zb[b];
    }
  }
  return res;
}

Return_t CoulombPotential::evaluate_spAA(const DistanceTableAA& d,
                                         const ParticleScalar* restrict Z,
                                         Vector<RealType>& va_sample,
                                         const std::vector<ListenerVector<RealType>>& listeners)
{
  Return_t res = 0.0;
  Return_t pairpot;
  int num_centers = d.centers();
  std::fill(va_sample.begin(), va_sample.end(), 0.0);

  for (size_t iat = 1; iat < num_centers; ++iat)
  {
    const auto& dist = d.getDistRow(iat);
    Return_t q       = Z[iat];
    for (size_t j = 0; j < iat; ++j)
    {
      pairpot = 0.5 * q * Z[j] / dist[j];
      va_sample[iat] += pairpot;
      va_sample[j] += pairpot;
      res += pairpot;
    }
  }
  res *= 2.0;
  return res;
}

Return_t CoulombPotential::evaluate_spAB(const DistanceTableAB& d,
                                         const ParticleScalar* restrict Za,
                                         const ParticleScalar* restrict Zb,
                                         Vector<RealType>& va_sample,
                                         Vector<RealType>& vb_sample,
                                         const std::vector<ListenerVector<RealType>>& listeners,
                                         const std::vector<ListenerVector<RealType>>& ion_listeners)
{
  Return_t res = 0.0;
  Return_t pairpot;
  std::fill(va_sample.begin(), va_sample.end(), 0.0);
  std::fill(vb_sample.begin(), vb_sample.end(), 0.0);
  const size_t num_targets = d.targets();
  const int num_sources    = d.sources();
  for (size_t b = 0; b < num_targets; ++b)
  {
    const auto& dist = d.getDistRow(b);
    Return_t z       = 0.5 * Zb[b];
    for (size_t a = 0; a < num_sources; ++a)
    {
      pairpot = z * Za[a] / dist[a];
      vb_sample[b] += pairpot;
      va_sample[a] += pairpot;
      res += pairpot;
    }
  }
  res *= 2.0;
  return res;
}

#if !defined(REMOVE_TRACEMANAGER)
/** evaluate AA-type interactions */
Return_t CoulombPotential::evaluate_spAA(const DistanceTableAA& d, const ParticleScalar* restrict Z)
{
  Return_t res = 0.0;
  Return_t pairpot;
  Array<RealType, 1>& Va_samp = *Va_sample;
  Va_samp                     = 0.0;
  for (size_t iat = 1; iat < nCenters; ++iat)
  {
    const auto& dist = d.getDistRow(iat);
    Return_t q       = Z[iat];
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
  Return_t Vnow        = res;
  Return_t Vsum        = Va_samp.sum();
  Return_t Vorig       = evaluateAA(d, Z);
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


Return_t CoulombPotential::evaluate_spAB(const DistanceTableAB& d,
                                         const ParticleScalar* restrict Za,
                                         const ParticleScalar* restrict Zb)
{
  Return_t res = 0.0;
  Return_t pairpot;
  Array<RealType, 1>& Va_samp = *Va_sample;
  Array<RealType, 1>& Vb_samp = *Vb_sample;
  Va_samp                     = 0.0;
  Vb_samp                     = 0.0;
  const size_t nTargets       = d.targets();
  for (size_t b = 0; b < nTargets; ++b)
  {
    const auto& dist = d.getDistRow(b);
    Return_t z       = 0.5 * Zb[b];
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
  Return_t Vnow        = res;
  Return_t Vasum       = Va_samp.sum();
  Return_t Vbsum       = Vb_samp.sum();
  Return_t Vsum        = Vasum + Vbsum;
  Return_t Vorig       = evaluateAB(d, Za, Zb);
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

void CoulombPotential::updateSource(ParticleSet& s)
{
  if (is_AA)
  {
    value_ = evaluateAA(s.getDistTableAA(myTableIndex), s.Z.first_address());
  }
}

Return_t CoulombPotential::evaluate(ParticleSet& P)
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

void CoulombPotential::mw_evaluatePerParticle(const RefVectorWithLeader<OperatorBase>& o_list,
                                              const RefVectorWithLeader<ParticleSet>& p_list,
                                              const std::vector<ListenerVector<RealType>>& listeners,
                                              const std::vector<ListenerVector<RealType>>& ion_listeners) const
{
  auto& o_leader = o_list.getCastedLeader<CoulombPotential>();
  auto& p_leader = p_list.getLeader();
  assert(this == &o_list.getLeader());
  auto& mw_res                = o_leader.mw_res_handle_.getResource();
  Vector<RealType>& va_sample = mw_res.va_samples;
  Vector<RealType>& vb_sample = mw_res.vb_samples;
  auto is_AA                  = o_leader.is_AA;
  auto is_active              = o_leader.is_active;
  auto myTableIndex           = o_leader.myTableIndex;
  auto nCenters               = o_leader.nCenters;

  auto evaluate_walker = [&va_sample, &vb_sample, myTableIndex, is_AA,
                          is_active](int walker_index, const ParticleSet& pset, const ParticleSet& pset_ions,
                                     const std::vector<ListenerVector<RealType>>& listeners,
                                     const std::vector<ListenerVector<RealType>>& ion_listeners) -> Return_t {
    Return_t value = 0;
    if (is_AA)
      if (is_active)
        value = evaluate_spAA(pset.getDistTableAA(myTableIndex), pset.Z.first_address(), va_sample, listeners);
      else
        value = evaluate_spAA(pset.getDistTableAA(myTableIndex), pset.Z.first_address(), va_sample, ion_listeners);
    else
      value = evaluate_spAB(pset.getDistTableAB(myTableIndex), pset_ions.Z.first_address(), pset.Z.first_address(),
                            va_sample, vb_sample, listeners, ion_listeners);
    return value;
  };

  auto name = o_leader.name_;
  for (int iw = 0; iw < o_list.size(); ++iw)
  {
    auto& coulomb_pot = o_list.getCastedElement<CoulombPotential>(iw);
    if (coulomb_pot.is_active)
    {
      if (is_AA)
        evaluate_walker(iw, p_list[iw], coulomb_pot.Pa, listeners, ion_listeners);
      else
        evaluate_walker(iw, p_list[iw], coulomb_pot.Pa, listeners, ion_listeners);
    }
    if (coulomb_pot.is_active)
    {
      if (is_AA)
        for (const ListenerVector<RealType>& listener : listeners)
          listener.report(iw, name, va_sample);
      else
        for (const ListenerVector<RealType>& listener : listeners)
          listener.report(iw, name, vb_sample);
      if (!is_AA)
        for (const ListenerVector<RealType>& ion_listener : ion_listeners)
          ion_listener.report(iw, name, va_sample);
    }
    else // its my belief that the only case here is ion AA.
      for (const ListenerVector<RealType>& ion_listener : ion_listeners)
        ion_listener.report(iw, name, va_sample);
  }
}

void CoulombPotential::evaluateIonDerivs(ParticleSet& P,
                                         ParticleSet& ions,
                                         TrialWaveFunction& psi,
                                         ParticleSet::ParticlePos& hf_terms,
                                         ParticleSet::ParticlePos& pulay_terms)
{
  if (!is_active)
    hf_terms -= forces_;
  // No Pulay here
}

bool CoulombPotential::put(xmlNodePtr cur) { return true; }

bool CoulombPotential::get(std::ostream& os) const
{
  if (myTableIndex)
    os << "CoulombAB source=" << Pa.getName() << std::endl;
  else
    os << "CoulombAA source/target " << Pa.getName() << std::endl;
  return true;
}

void CoulombPotential::setObservables(PropertySetType& plist)
{
  OperatorBase::setObservables(plist);
  if (ComputeForces)
    setObservablesF(plist);
}

void CoulombPotential::setParticlePropertyList(PropertySetType& plist, int offset)
{
  OperatorBase::setParticlePropertyList(plist, offset);
  if (ComputeForces)
    setParticleSetF(plist, offset);
}


std::unique_ptr<OperatorBase> CoulombPotential::makeClone(ParticleSet& qp)
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

} // namespace qmcplusplus
