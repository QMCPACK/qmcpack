//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "NonLocalECPotential.h"
#include <DistanceTable.h>
#include <IteratorUtility.h>
#include <ResourceCollection.h>
#include "NonLocalECPComponent.h"
#include "NLPPJob.h"

namespace qmcplusplus
{
struct NonLocalECPotentialMultiWalkerResource : public Resource
{
  NonLocalECPotentialMultiWalkerResource() : Resource("NonLocalECPotential"), collection("NLPPcollection") {}

  Resource* makeClone() const override { return new NonLocalECPotentialMultiWalkerResource(*this); }

  ResourceCollection collection;
};

void NonLocalECPotential::resetTargetParticleSet(ParticleSet& P) {}

/** constructor
 *\param ions the positions of the ions
 *\param els the positions of the electrons
 *\param psi trial wavefunction
 */
NonLocalECPotential::NonLocalECPotential(ParticleSet& ions,
                                         ParticleSet& els,
                                         TrialWaveFunction& psi,
                                         bool computeForces,
                                         bool enable_DLA)
    : ForceBase(ions, els),
      myRNG(nullptr),
      IonConfig(ions),
      Psi(psi),
      ComputeForces(computeForces),
      use_DLA(enable_DLA),
      Peln(els),
      ElecNeighborIons(els),
      IonNeighborElecs(ions),
      UseTMove(TMOVE_OFF)
{
  setEnergyDomain(POTENTIAL);
  twoBodyQuantumDomain(ions, els);
  myTableIndex = els.addTable(ions);
  NumIons      = ions.getTotalNum();
  //els.resizeSphere(NumIons);
  PP.resize(NumIons, nullptr);
  prefix = "FNL";
  PPset.resize(IonConfig.getSpeciesSet().getTotalNum());
  PulayTerm.resize(NumIons);
  update_mode_.set(NONLOCAL, 1);
  nlpp_jobs.resize(els.groups());
  for (size_t ig = 0; ig < els.groups(); ig++)
  {
    // this should be enough in most calculations assuming that every electron cannot be in more than two pseudo regions.
    nlpp_jobs[ig].reserve(2 * els.groupsize(ig));
  }
}

NonLocalECPotential::~NonLocalECPotential() = default;

#if !defined(REMOVE_TRACEMANAGER)
void NonLocalECPotential::contributeParticleQuantities() { request_.contribute_array(name_); }

void NonLocalECPotential::checkoutParticleQuantities(TraceManager& tm)
{
  streaming_particles_ = request_.streaming_array(name_);
  if (streaming_particles_)
  {
    Ve_sample = tm.checkout_real<1>(name_, Peln);
    Vi_sample = tm.checkout_real<1>(name_, IonConfig);
  }
}

void NonLocalECPotential::deleteParticleQuantities()
{
  if (streaming_particles_)
  {
    delete Ve_sample;
    delete Vi_sample;
  }
}
#endif

NonLocalECPotential::Return_t NonLocalECPotential::evaluate(ParticleSet& P)
{
  evaluateImpl(P, false);
  return value_;
}

NonLocalECPotential::Return_t NonLocalECPotential::evaluateDeterministic(ParticleSet& P)
{
  evaluateImpl(P, false, true);
  return value_;
}

void NonLocalECPotential::mw_evaluate(const RefVectorWithLeader<OperatorBase>& o_list,
                                      const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                      const RefVectorWithLeader<ParticleSet>& p_list) const
{
  mw_evaluateImpl(o_list, wf_list, p_list, false);
}

NonLocalECPotential::Return_t NonLocalECPotential::evaluateWithToperator(ParticleSet& P)
{
  if (UseTMove == TMOVE_V0 || UseTMove == TMOVE_V3)
    evaluateImpl(P, true);
  else
    evaluateImpl(P, false);
  return value_;
}

void NonLocalECPotential::mw_evaluateWithToperator(const RefVectorWithLeader<OperatorBase>& o_list,
                                                   const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                                   const RefVectorWithLeader<ParticleSet>& p_list) const
{
  if (UseTMove == TMOVE_V0 || UseTMove == TMOVE_V3)
    mw_evaluateImpl(o_list, wf_list, p_list, true);
  else
    mw_evaluateImpl(o_list, wf_list, p_list, false);
}

void NonLocalECPotential::evaluateImpl(ParticleSet& P, bool Tmove, bool keepGrid)
{
  if (Tmove)
    tmove_xy_.clear();

  value_ = 0.0;
#if !defined(REMOVE_TRACEMANAGER)
  auto& Ve_samp = *Ve_sample;
  auto& Vi_samp = *Vi_sample;
  if (streaming_particles_)
  {
    Ve_samp = 0.0;
    Vi_samp = 0.0;
  }
#endif
  for (int ipp = 0; ipp < PPset.size(); ipp++)
    if (PPset[ipp])
      if (!keepGrid)
        PPset[ipp]->randomize_grid(*myRNG);
  //loop over all the ions
  const auto& myTable = P.getDistTableAB(myTableIndex);
  // clear all the electron and ion neighbor lists
  for (int iat = 0; iat < NumIons; iat++)
    IonNeighborElecs.getNeighborList(iat).clear();
  for (int jel = 0; jel < P.getTotalNum(); jel++)
    ElecNeighborIons.getNeighborList(jel).clear();

  if (ComputeForces)
  {
    forces = 0;
    for (int ig = 0; ig < P.groups(); ++ig) //loop over species
    {
      Psi.prepareGroup(P, ig);
      for (int jel = P.first(ig); jel < P.last(ig); ++jel)
      {
        const auto& dist               = myTable.getDistRow(jel);
        const auto& displ              = myTable.getDisplRow(jel);
        std::vector<int>& NeighborIons = ElecNeighborIons.getNeighborList(jel);
        for (int iat = 0; iat < NumIons; iat++)
          if (PP[iat] != nullptr && dist[iat] < PP[iat]->getRmax())
          {
            RealType pairpot = PP[iat]->evaluateOneWithForces(P, iat, Psi, jel, dist[iat], -displ[iat], forces[iat]);
            if (Tmove)
              PP[iat]->contributeTxy(jel, tmove_xy_);
            value_ += pairpot;
            NeighborIons.push_back(iat);
            IonNeighborElecs.getNeighborList(iat).push_back(jel);
          }
      }
    }
  }
  else
  {
    for (int ig = 0; ig < P.groups(); ++ig) //loop over species
    {
      Psi.prepareGroup(P, ig);
      for (int jel = P.first(ig); jel < P.last(ig); ++jel)
      {
        const auto& dist               = myTable.getDistRow(jel);
        const auto& displ              = myTable.getDisplRow(jel);
        std::vector<int>& NeighborIons = ElecNeighborIons.getNeighborList(jel);
        for (int iat = 0; iat < NumIons; iat++)
          if (PP[iat] != nullptr && dist[iat] < PP[iat]->getRmax())
          {
            RealType pairpot = PP[iat]->evaluateOne(P, iat, Psi, jel, dist[iat], -displ[iat], use_DLA);
            if (Tmove)
              PP[iat]->contributeTxy(jel, tmove_xy_);
            value_ += pairpot;
            NeighborIons.push_back(iat);
            IonNeighborElecs.getNeighborList(iat).push_back(jel);
            if (streaming_particles_)
            {
              Ve_samp(jel) += 0.5 * pairpot;
              Vi_samp(iat) += 0.5 * pairpot;
            }
          }
      }
    }
  }

#if !defined(TRACE_CHECK)
  if (streaming_particles_)
  {
    Return_t Vnow  = value_;
    RealType Visum = Vi_sample->sum();
    RealType Vesum = Ve_sample->sum();
    RealType Vsum  = Vesum + Visum;
    if (std::abs(Vsum - Vnow) > TraceManager::trace_tol)
    {
      app_log() << "accumtest: NonLocalECPotential::evaluate()" << std::endl;
      app_log() << "accumtest:   tot:" << Vnow << std::endl;
      app_log() << "accumtest:   sum:" << Vsum << std::endl;
      APP_ABORT("Trace check failed");
    }
    if (std::abs(Vesum - Visum) > TraceManager::trace_tol)
    {
      app_log() << "sharetest: NonLocalECPotential::evaluate()" << std::endl;
      app_log() << "sharetest:   e share:" << Vesum << std::endl;
      app_log() << "sharetest:   i share:" << Visum << std::endl;
      APP_ABORT("Trace check failed");
    }
  }
#endif
}

void NonLocalECPotential::mw_evaluateImpl(const RefVectorWithLeader<OperatorBase>& o_list,
                                          const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                          const RefVectorWithLeader<ParticleSet>& p_list,
                                          bool Tmove)
{
  auto& O_leader           = o_list.getCastedLeader<NonLocalECPotential>();
  ParticleSet& pset_leader = p_list.getLeader();
  const size_t nw          = o_list.size();

  if (O_leader.ComputeForces)
    APP_ABORT("NonLocalECPotential::mw_evaluateImpl(): Forces not imlpemented\n");

  for (size_t iw = 0; iw < nw; iw++)
  {
    auto& O = o_list.getCastedElement<NonLocalECPotential>(iw);
    const ParticleSet& P(p_list[iw]);

    if (Tmove)
      O.tmove_xy_.clear();

    for (int ipp = 0; ipp < O.PPset.size(); ipp++)
      if (O.PPset[ipp])
        O.PPset[ipp]->randomize_grid(*O.myRNG);

    //loop over all the ions
    const auto& myTable = P.getDistTableAB(O.myTableIndex);
    // clear all the electron and ion neighbor lists
    for (int iat = 0; iat < O.NumIons; iat++)
      O.IonNeighborElecs.getNeighborList(iat).clear();
    for (int jel = 0; jel < P.getTotalNum(); jel++)
      O.ElecNeighborIons.getNeighborList(jel).clear();

    for (int ig = 0; ig < P.groups(); ++ig) //loop over species
    {
      auto& joblist = O.nlpp_jobs[ig];
      joblist.clear();

      for (int jel = P.first(ig); jel < P.last(ig); ++jel)
      {
        const auto& dist               = myTable.getDistRow(jel);
        const auto& displ              = myTable.getDisplRow(jel);
        std::vector<int>& NeighborIons = O.ElecNeighborIons.getNeighborList(jel);
        for (int iat = 0; iat < O.NumIons; iat++)
          if (O.PP[iat] != nullptr && dist[iat] < O.PP[iat]->getRmax())
          {
            NeighborIons.push_back(iat);
            O.IonNeighborElecs.getNeighborList(iat).push_back(jel);
            joblist.emplace_back(iat, jel, P.R[jel], dist[iat], -displ[iat]);
          }
      }
    }

    O.value_ = 0.0;
  }

  // make this class unit tests friendly without the need of setup resources.
  if (!O_leader.mw_res_)
  {
    app_warning() << "NonLocalECPotential: This message should not be seen in production (performance bug) runs "
                     "but only unit tests (expected)."
                  << std::endl;
    O_leader.mw_res_ = std::make_unique<NonLocalECPotentialMultiWalkerResource>();
    for (int ig = 0; ig < O_leader.PPset.size(); ++ig)
      if (O_leader.PPset[ig]->getVP())
      {
        O_leader.PPset[ig]->getVP()->createResource(O_leader.mw_res_->collection);
        break;
      }
  }

  auto pp_component = std::find_if(O_leader.PPset.begin(), O_leader.PPset.end(), [](auto& ptr) { return bool(ptr); });
  assert(pp_component != std::end(O_leader.PPset));

  RefVector<NonLocalECPotential> ecp_potential_list;
  RefVectorWithLeader<NonLocalECPComponent> ecp_component_list(**pp_component);
  RefVectorWithLeader<ParticleSet> pset_list(pset_leader);
  RefVectorWithLeader<TrialWaveFunction> psi_list(O_leader.Psi);
  // we are moving away from internally stored Psi, double check before Psi gets finally removed.
  assert(&O_leader.Psi == &wf_list.getLeader());
  for (size_t iw = 0; iw < nw; iw++)
    assert(&o_list.getCastedElement<NonLocalECPotential>(iw).Psi == &wf_list[iw]);

  RefVector<const NLPPJob<RealType>> batch_list;
  std::vector<RealType> pairpots(nw);

  ecp_potential_list.reserve(nw);
  ecp_component_list.reserve(nw);
  pset_list.reserve(nw);
  psi_list.reserve(nw);
  batch_list.reserve(nw);

  for (int ig = 0; ig < pset_leader.groups(); ++ig) //loop over species
  {
    TrialWaveFunction::mw_prepareGroup(wf_list, p_list, ig);

    // find the max number of jobs of all the walkers
    size_t max_num_jobs = 0;
    for (size_t iw = 0; iw < nw; iw++)
    {
      const auto& O = o_list.getCastedElement<NonLocalECPotential>(iw);
      max_num_jobs  = std::max(max_num_jobs, O.nlpp_jobs[ig].size());
    }

    for (size_t jobid = 0; jobid < max_num_jobs; jobid++)
    {
      ecp_potential_list.clear();
      ecp_component_list.clear();
      pset_list.clear();
      psi_list.clear();
      batch_list.clear();
      for (size_t iw = 0; iw < nw; iw++)
      {
        auto& O = o_list.getCastedElement<NonLocalECPotential>(iw);
        if (jobid < O.nlpp_jobs[ig].size())
        {
          const auto& job = O.nlpp_jobs[ig][jobid];
          ecp_potential_list.push_back(O);
          ecp_component_list.push_back(*O.PP[job.ion_id]);
          pset_list.push_back(p_list[iw]);
          psi_list.push_back(wf_list[iw]);
          batch_list.push_back(job);
        }
      }

      NonLocalECPComponent::mw_evaluateOne(ecp_component_list, pset_list, psi_list, batch_list, pairpots,
                                           O_leader.mw_res_->collection, O_leader.use_DLA);

      for (size_t j = 0; j < ecp_potential_list.size(); j++)
      {
        if (false)
        { // code usefully for debugging
          RealType check_value =
              ecp_component_list[j].evaluateOne(pset_list[j], batch_list[j].get().ion_id, psi_list[j],
                                                batch_list[j].get().electron_id, batch_list[j].get().ion_elec_dist,
                                                batch_list[j].get().ion_elec_displ, O_leader.use_DLA);
          if (std::abs(check_value - pairpots[j]) > 1e-5)
            std::cout << "check " << check_value << " wrong " << pairpots[j] << " diff "
                      << std::abs(check_value - pairpots[j]) << std::endl;
        }
        ecp_potential_list[j].get().value_ += pairpots[j];
        if (Tmove)
          ecp_component_list[j].contributeTxy(batch_list[j].get().electron_id, ecp_potential_list[j].get().tmove_xy_);
      }
    }
  }
}


void NonLocalECPotential::evalIonDerivsImpl(ParticleSet& P,
                                            ParticleSet& ions,
                                            TrialWaveFunction& psi,
                                            ParticleSet::ParticlePos& hf_terms,
                                            ParticleSet::ParticlePos& pulay_terms,
                                            bool keepGrid)
{
  //We're going to ignore psi and use the internal Psi.
  //
  //Dummy vector for now.  Tmoves not implemented
  bool Tmove = false;

  forces    = 0;
  PulayTerm = 0;

  value_ = 0.0;
  if (!keepGrid)
  {
    for (int ipp = 0; ipp < PPset.size(); ipp++)
      if (PPset[ipp])
        PPset[ipp]->randomize_grid(*myRNG);
  }
  //loop over all the ions
  const auto& myTable = P.getDistTableAB(myTableIndex);
  // clear all the electron and ion neighbor lists
  for (int iat = 0; iat < NumIons; iat++)
    IonNeighborElecs.getNeighborList(iat).clear();
  for (int jel = 0; jel < P.getTotalNum(); jel++)
    ElecNeighborIons.getNeighborList(jel).clear();

  for (int ig = 0; ig < P.groups(); ++ig) //loop over species
  {
    Psi.prepareGroup(P, ig);
    for (int jel = P.first(ig); jel < P.last(ig); ++jel)
    {
      const auto& dist               = myTable.getDistRow(jel);
      const auto& displ              = myTable.getDisplRow(jel);
      std::vector<int>& NeighborIons = ElecNeighborIons.getNeighborList(jel);
      for (int iat = 0; iat < NumIons; iat++)
        if (PP[iat] != nullptr && dist[iat] < PP[iat]->getRmax())
        {
          value_ +=
              PP[iat]->evaluateOneWithForces(P, ions, iat, Psi, jel, dist[iat], -displ[iat], forces[iat], PulayTerm);
          if (Tmove)
            PP[iat]->contributeTxy(jel, tmove_xy_);
          NeighborIons.push_back(iat);
          IonNeighborElecs.getNeighborList(iat).push_back(jel);
        }
    }
  }

  hf_terms -= forces;
  pulay_terms -= PulayTerm;
}

NonLocalECPotential::Return_t NonLocalECPotential::evaluateWithIonDerivs(ParticleSet& P,
                                                                         ParticleSet& ions,
                                                                         TrialWaveFunction& psi,
                                                                         ParticleSet::ParticlePos& hf_terms,
                                                                         ParticleSet::ParticlePos& pulay_terms)
{
  evalIonDerivsImpl(P, ions, psi, hf_terms, pulay_terms);
  return value_;
}

NonLocalECPotential::Return_t NonLocalECPotential::evaluateWithIonDerivsDeterministic(
    ParticleSet& P,
    ParticleSet& ions,
    TrialWaveFunction& psi,
    ParticleSet::ParticlePos& hf_terms,
    ParticleSet::ParticlePos& pulay_terms)
{
  evalIonDerivsImpl(P, ions, psi, hf_terms, pulay_terms, true);
  return value_;
}

void NonLocalECPotential::computeOneElectronTxy(ParticleSet& P, const int ref_elec)
{
  tmove_xy_.clear();
  const auto& myTable                  = P.getDistTableAB(myTableIndex);
  const std::vector<int>& NeighborIons = ElecNeighborIons.getNeighborList(ref_elec);

  const auto& dist  = myTable.getDistRow(ref_elec);
  const auto& displ = myTable.getDisplRow(ref_elec);
  for (int atom_index = 0; atom_index < NeighborIons.size(); atom_index++)
  {
    const int iat = NeighborIons[atom_index];
    PP[iat]->evaluateOne(P, iat, Psi, ref_elec, dist[iat], -displ[iat], use_DLA);
    PP[iat]->contributeTxy(ref_elec, tmove_xy_);
  }
}

void NonLocalECPotential::evaluateOneBodyOpMatrix(ParticleSet& P,
                                                  const TWFFastDerivWrapper& psi,
                                                  std::vector<ValueMatrix>& B)
{
  bool keepGrid = true;
  for (int ipp = 0; ipp < PPset.size(); ipp++)
    if (PPset[ipp])
      if (!keepGrid)
        PPset[ipp]->randomize_grid(*myRNG);
  //loop over all the ions
  const auto& myTable = P.getDistTableAB(myTableIndex);
  // clear all the electron and ion neighbor lists
  for (int iat = 0; iat < NumIons; iat++)
    IonNeighborElecs.getNeighborList(iat).clear();
  for (int jel = 0; jel < P.getTotalNum(); jel++)
    ElecNeighborIons.getNeighborList(jel).clear();

  for (int ig = 0; ig < P.groups(); ++ig) //loop over species
  {
    for (int jel = P.first(ig); jel < P.last(ig); ++jel)
    {
      const auto& dist               = myTable.getDistRow(jel);
      const auto& displ              = myTable.getDisplRow(jel);
      auto& NeighborIons = ElecNeighborIons.getNeighborList(jel);
      for (int iat = 0; iat < NumIons; iat++)
        if (PP[iat] != nullptr && dist[iat] < PP[iat]->getRmax())
        {
          PP[iat]->evaluateOneBodyOpMatrixContribution(P, iat, psi, jel, dist[iat], -displ[iat], B);
          NeighborIons.push_back(iat);
          IonNeighborElecs.getNeighborList(iat).push_back(jel);
        }
    }
  }
}

void NonLocalECPotential::evaluateOneBodyOpMatrixForceDeriv(ParticleSet& P,
                                                            const ParticleSet& source,
                                                            const TWFFastDerivWrapper& psi,
                                                            const int iat_source,
                                                            std::vector<std::vector<ValueMatrix>>& Bforce)
{
  bool keepGrid = true;
  for (int ipp = 0; ipp < PPset.size(); ipp++)
    if (PPset[ipp])
      if (!keepGrid)
        PPset[ipp]->randomize_grid(*myRNG);
  //loop over all the ions
  const auto& myTable = P.getDistTableAB(myTableIndex);
  // clear all the electron and ion neighbor lists
  for (int iat = 0; iat < NumIons; iat++)
    IonNeighborElecs.getNeighborList(iat).clear();
  for (int jel = 0; jel < P.getTotalNum(); jel++)
    ElecNeighborIons.getNeighborList(jel).clear();

  for (int ig = 0; ig < P.groups(); ++ig) //loop over species
  {
    for (int jel = P.first(ig); jel < P.last(ig); ++jel)
    {
      const auto& dist               = myTable.getDistRow(jel);
      const auto& displ              = myTable.getDisplRow(jel);
      auto& NeighborIons = ElecNeighborIons.getNeighborList(jel);
      for (int iat = 0; iat < NumIons; iat++)
        if (PP[iat] != nullptr && dist[iat] < PP[iat]->getRmax())
        {
          PP[iat]->evaluateOneBodyOpMatrixdRContribution(P, source, iat, iat_source, psi, jel, dist[iat], -displ[iat],
                                                         Bforce);
          NeighborIons.push_back(iat);
          IonNeighborElecs.getNeighborList(iat).push_back(jel);
        }
    }
  }
}
int NonLocalECPotential::makeNonLocalMovesPbyP(ParticleSet& P)
{
  int NonLocalMoveAccepted = 0;
  RandomGenerator& RandomGen(*myRNG);
  if (UseTMove == TMOVE_V0)
  {
    const NonLocalData* oneTMove = nonLocalOps.selectMove(RandomGen(), tmove_xy_);
    //make a non-local move
    if (oneTMove)
    {
      int iat = oneTMove->PID;
      Psi.prepareGroup(P, P.getGroupID(iat));
      if (P.makeMoveAndCheck(iat, oneTMove->Delta))
      {
        GradType grad_iat;
        Psi.calcRatioGrad(P, iat, grad_iat);
        Psi.acceptMove(P, iat, true);
        P.acceptMove(iat);
        NonLocalMoveAccepted++;
      }
    }
  }
  else if (UseTMove == TMOVE_V1)
  {
    GradType grad_iat;
    //make a non-local move per particle
    for (int ig = 0; ig < P.groups(); ++ig) //loop over species
    {
      Psi.prepareGroup(P, ig);
      for (int iat = P.first(ig); iat < P.last(ig); ++iat)
      {
        computeOneElectronTxy(P, iat);
        const NonLocalData* oneTMove = nonLocalOps.selectMove(RandomGen(), tmove_xy_);
        if (oneTMove)
        {
          if (P.makeMoveAndCheck(iat, oneTMove->Delta))
          {
            Psi.calcRatioGrad(P, iat, grad_iat);
            Psi.acceptMove(P, iat, true);
            P.acceptMove(iat);
            NonLocalMoveAccepted++;
          }
        }
      }
    }
  }
  else if (UseTMove == TMOVE_V3)
  {
    elecTMAffected.assign(P.getTotalNum(), false);
    nonLocalOps.groupByElectron(P.getTotalNum(), tmove_xy_);
    GradType grad_iat;
    //make a non-local move per particle
    for (int ig = 0; ig < P.groups(); ++ig) //loop over species
    {
      Psi.prepareGroup(P, ig);
      for (int iat = P.first(ig); iat < P.last(ig); ++iat)
      {
        const NonLocalData* oneTMove;
        if (elecTMAffected[iat])
        {
          // recompute Txy for the given electron effected by Tmoves
          computeOneElectronTxy(P, iat);
          oneTMove = nonLocalOps.selectMove(RandomGen(), tmove_xy_);
        }
        else
          oneTMove = nonLocalOps.selectMove(RandomGen(), iat);
        if (oneTMove)
        {
          if (P.makeMoveAndCheck(iat, oneTMove->Delta))
          {
            Psi.calcRatioGrad(P, iat, grad_iat);
            Psi.acceptMove(P, iat, true);
            // mark all affected electrons
            markAffectedElecs(P.getDistTableAB(myTableIndex), iat);
            P.acceptMove(iat);
            NonLocalMoveAccepted++;
          }
        }
      }
    }
  }

  if (NonLocalMoveAccepted > 0)
  {
    Psi.completeUpdates();
    // this step also updates electron positions on the device.
    P.donePbyP(true);
  }

  return NonLocalMoveAccepted;
}

void NonLocalECPotential::markAffectedElecs(const DistanceTableAB& myTable, int iel)
{
  std::vector<int>& NeighborIons = ElecNeighborIons.getNeighborList(iel);
  for (int iat = 0; iat < NumIons; iat++)
  {
    if (PP[iat] == nullptr)
      continue;
    RealType old_distance = 0.0;
    RealType new_distance = 0.0;
    old_distance          = myTable.getDistRow(iel)[iat];
    new_distance          = myTable.getTempDists()[iat];
    bool moved            = false;
    // move out
    if (old_distance < PP[iat]->getRmax() && new_distance >= PP[iat]->getRmax())
    {
      moved                           = true;
      std::vector<int>& NeighborElecs = IonNeighborElecs.getNeighborList(iat);
      auto iter_at                    = std::find(NeighborIons.begin(), NeighborIons.end(), iat);
      auto iter_el                    = std::find(NeighborElecs.begin(), NeighborElecs.end(), iel);
      *iter_at                        = NeighborIons.back();
      *iter_el                        = NeighborElecs.back();
      NeighborIons.pop_back();
      NeighborElecs.pop_back();
      elecTMAffected[iel] = true;
    }
    // move in
    if (old_distance >= PP[iat]->getRmax() && new_distance < PP[iat]->getRmax())
    {
      moved                           = true;
      std::vector<int>& NeighborElecs = IonNeighborElecs.getNeighborList(iat);
      NeighborElecs.push_back(iel);
      NeighborIons.push_back(iat);
    }
    // move around
    if (moved || (old_distance < PP[iat]->getRmax() && new_distance < PP[iat]->getRmax()))
    {
      std::vector<int>& NeighborElecs = IonNeighborElecs.getNeighborList(iat);
      for (int jel = 0; jel < NeighborElecs.size(); ++jel)
        elecTMAffected[NeighborElecs[jel]] = true;
    }
  }
}

void NonLocalECPotential::addComponent(int groupID, std::unique_ptr<NonLocalECPComponent>&& ppot)
{
  for (int iat = 0; iat < PP.size(); iat++)
    if (IonConfig.GroupID[iat] == groupID)
      PP[iat] = ppot.get();
  PPset[groupID] = std::move(ppot);
}

void NonLocalECPotential::createResource(ResourceCollection& collection) const
{
  auto new_res = std::make_unique<NonLocalECPotentialMultiWalkerResource>();
  for (int ig = 0; ig < PPset.size(); ++ig)
    if (PPset[ig]->getVP())
    {
      PPset[ig]->getVP()->createResource(new_res->collection);
      break;
    }
  auto resource_index = collection.addResource(std::move(new_res));
}

void NonLocalECPotential::acquireResource(ResourceCollection& collection,
                                          const RefVectorWithLeader<OperatorBase>& o_list) const
{
  auto& O_leader = o_list.getCastedLeader<NonLocalECPotential>();
  auto res_ptr   = dynamic_cast<NonLocalECPotentialMultiWalkerResource*>(collection.lendResource().release());
  if (!res_ptr)
    throw std::runtime_error("NonLocalECPotential::acquireResource dynamic_cast failed");
  O_leader.mw_res_.reset(res_ptr);
}

void NonLocalECPotential::releaseResource(ResourceCollection& collection,
                                          const RefVectorWithLeader<OperatorBase>& o_list) const
{
  auto& O_leader = o_list.getCastedLeader<NonLocalECPotential>();
  collection.takebackResource(std::move(O_leader.mw_res_));
}

std::unique_ptr<OperatorBase> NonLocalECPotential::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  std::unique_ptr<NonLocalECPotential> myclone =
      std::make_unique<NonLocalECPotential>(IonConfig, qp, psi, ComputeForces, use_DLA);
  for (int ig = 0; ig < PPset.size(); ++ig)
    if (PPset[ig])
      myclone->addComponent(ig, std::unique_ptr<NonLocalECPComponent>(PPset[ig]->makeClone(qp)));
  return myclone;
}


void NonLocalECPotential::addObservables(PropertySetType& plist, BufferType& collectables)
{
  OperatorBase::addValue(plist);
  if (ComputeForces)
  {
    if (FirstForceIndex < 0)
      FirstForceIndex = plist.size();
    for (int iat = 0; iat < Nnuc; iat++)
    {
      for (int x = 0; x < OHMMS_DIM; x++)
      {
        std::ostringstream obsName1, obsName2;
        obsName1 << "FNL"
                 << "_" << iat << "_" << x;
        plist.add(obsName1.str());
        //        obsName2 << "FNL_Pulay" << "_" << iat << "_" << x;
        //        plist.add(obsName2.str());
      }
    }
  }
}

void NonLocalECPotential::registerObservables(std::vector<ObservableHelper>& h5list, hid_t gid) const
{
  OperatorBase::registerObservables(h5list, gid);
  if (ComputeForces)
  {
    std::vector<int> ndim(2);
    ndim[0] = Nnuc;
    ndim[1] = OHMMS_DIM;
    h5list.emplace_back("FNL");
    auto& h5o1 = h5list.back();
    h5o1.set_dimensions(ndim, FirstForceIndex);
    h5o1.open(gid);
  }
}

void NonLocalECPotential::setObservables(QMCTraits::PropertySetType& plist)
{
  OperatorBase::setObservables(plist);
  if (ComputeForces)
  {
    int index = FirstForceIndex;
    for (int iat = 0; iat < Nnuc; iat++)
    {
      for (int x = 0; x < OHMMS_DIM; x++)
      {
        plist[index++] = forces[iat][x];
        //    plist[index++] = PulayTerm[iat][x];
      }
    }
  }
}


void NonLocalECPotential::setParticlePropertyList(QMCTraits::PropertySetType& plist, int offset)
{
  OperatorBase::setParticlePropertyList(plist, offset);
  if (ComputeForces)
  {
    int index = FirstForceIndex + offset;
    for (int iat = 0; iat < Nnuc; iat++)
    {
      for (int x = 0; x < OHMMS_DIM; x++)
      {
        plist[index++] = forces[iat][x];
        //        plist[index++] = PulayTerm[iat][x];
      }
    }
  }
}


} // namespace qmcplusplus
