//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//
// File created by: Cody A. Melton, cmelton@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "Particle/DistanceTable.h"
#include "SOECPotential.h"
#include "Utilities/IteratorUtility.h"
#include "NLPPJob.h"
#include <ResourceCollection.h>

namespace qmcplusplus
{

struct SOECPotential::SOECPotentialMultiWalkerResource : public Resource
{
  SOECPotentialMultiWalkerResource() : Resource("SOECPotential") {}

  Resource* makeClone() const override;

  ResourceCollection collection{"SOPPcollection"};
};

/** constructor
 *\param ionic positions
 *\param els electronic poitions
 *\param psi Trial wave function
*/
SOECPotential::SOECPotential(ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi)
    : myRNG_(nullptr), IonConfig_(ions), Psi_(psi), Peln_(els), ElecNeighborIons_(els), IonNeighborElecs_(ions)
{
  setEnergyDomain(POTENTIAL);
  twoBodyQuantumDomain(ions, els);
  myTableIndex_ = els.addTable(ions);
  NumIons_      = ions.getTotalNum();
  PP_.resize(NumIons_, nullptr);
  PPset_.resize(IonConfig_.getSpeciesSet().getTotalNum());
  sopp_jobs_.resize(els.groups());
  for (size_t ig = 0; ig < els.groups(); ig++)
    sopp_jobs_[ig].reserve(2 * els.groupsize(ig));
}

SOECPotential::~SOECPotential() = default;

void SOECPotential::resetTargetParticleSet(ParticleSet& P) {}

SOECPotential::Return_t SOECPotential::evaluate(ParticleSet& P)
{
  value_ = 0.0;
  for (int ipp = 0; ipp < PPset_.size(); ipp++)
    if (PPset_[ipp])
      PPset_[ipp]->rotateQuadratureGrid(generateRandomRotationMatrix(*myRNG_));

  const auto& myTable = P.getDistTableAB(myTableIndex_);
  for (int iat = 0; iat < NumIons_; iat++)
    IonNeighborElecs_.getNeighborList(iat).clear();
  for (int jel = 0; jel < P.getTotalNum(); jel++)
    ElecNeighborIons_.getNeighborList(jel).clear();

  for (int jel = 0; jel < P.getTotalNum(); jel++)
  {
    const auto& dist               = myTable.getDistRow(jel);
    const auto& displ              = myTable.getDisplRow(jel);
    std::vector<int>& NeighborIons = ElecNeighborIons_.getNeighborList(jel);
    for (int iat = 0; iat < NumIons_; iat++)
      if (PP_[iat] != nullptr && dist[iat] < PP_[iat]->getRmax())
      {
        RealType pairpot = PP_[iat]->evaluateOne(P, iat, Psi_, jel, dist[iat], -displ[iat]);
        value_ += pairpot;
        NeighborIons.push_back(iat);
        IonNeighborElecs_.getNeighborList(iat).push_back(jel);
      }
  }
  return value_;
}

SOECPotential::Return_t SOECPotential::evaluateValueAndDerivatives(ParticleSet& P,
                                                                   const opt_variables_type& optvars,
                                                                   const Vector<ValueType>& dlogpsi,
                                                                   Vector<ValueType>& dhpsioverpsi)
{
  value_ = 0.0;
  for (int ipp = 0; ipp < PPset_.size(); ipp++)
    if (PPset_[ipp])
      PPset_[ipp]->rotateQuadratureGrid(generateRandomRotationMatrix(*myRNG_));

  const auto& myTable = P.getDistTableAB(myTableIndex_);
  for (int jel = 0; jel < P.getTotalNum(); jel++)
  {
    const auto& dist  = myTable.getDistRow(jel);
    const auto& displ = myTable.getDisplRow(jel);
    for (int iat = 0; iat < NumIons_; iat++)
      if (PP_[iat] != nullptr && dist[iat] < PP_[iat]->getRmax())
        value_ += PP_[iat]->evaluateValueAndDerivatives(P, iat, Psi_, jel, dist[iat], -displ[iat], optvars, dlogpsi,
                                                        dhpsioverpsi);
  }
  return value_;
}

void SOECPotential::mw_evaluate(const RefVectorWithLeader<OperatorBase>& o_list,
                                const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                const RefVectorWithLeader<ParticleSet>& p_list) const
{
  auto& O_leader           = o_list.getCastedLeader<SOECPotential>();
  ParticleSet& pset_leader = p_list.getLeader();
  const size_t nw          = o_list.size();

  for (size_t iw = 0; iw < nw; iw++)
  {
    auto& O = o_list.getCastedElement<SOECPotential>(iw);
    const ParticleSet& P(p_list[iw]);

    for (size_t ipp = 0; ipp < O.PPset_.size(); ipp++)
      if (O.PPset_[ipp])
        O.PPset_[ipp]->rotateQuadratureGrid(generateRandomRotationMatrix(*O.myRNG_));

    //loop over all the ions
    const auto& myTable = P.getDistTableAB(O.myTableIndex_);
    //clear elec and ion neighbor lists
    for (size_t iat = 0; iat < O.NumIons_; iat++)
      O.IonNeighborElecs_.getNeighborList(iat).clear();
    for (size_t jel = 0; jel < P.getTotalNum(); jel++)
      O.ElecNeighborIons_.getNeighborList(jel).clear();

    for (size_t ig = 0; ig < P.groups(); ig++) // loop over species
    {
      auto& joblist = O.sopp_jobs_[ig];
      joblist.clear();

      for (size_t jel = P.first(ig); jel < P.last(ig); jel++)
      {
        const auto& dist               = myTable.getDistRow(jel);
        const auto& displ              = myTable.getDisplRow(jel);
        std::vector<int>& NeighborIons = O.ElecNeighborIons_.getNeighborList(jel);
        for (size_t iat = 0; iat < O.NumIons_; iat++)
          if (O.PP_[iat] != nullptr && dist[iat] < O.PP_[iat]->getRmax())
          {
            NeighborIons.push_back(iat);
            O.IonNeighborElecs_.getNeighborList(iat).push_back(jel);
            joblist.emplace_back(iat, jel, dist[iat], -displ[iat]);
          }
      }
    }
    O.value_ = 0.0;
  }

  auto pp_component = std::find_if(O_leader.PPset_.begin(), O_leader.PPset_.end(), [](auto& ptr) { return bool(ptr); });
  assert(pp_component != std::end(O_leader.PPset_));

  RefVector<SOECPotential> soecp_potential_list;
  RefVectorWithLeader<SOECPComponent> soecp_component_list(**pp_component);
  RefVectorWithLeader<ParticleSet> pset_list(pset_leader);
  RefVectorWithLeader<TrialWaveFunction> psi_list(O_leader.Psi_);
  assert(&O_leader.Psi == &wf_list.getLeader());
  for (size_t iw = 0; iw < nw; iw++)
    assert(&o_list.getCastedElement<SOECPotential>(iw).Psi == &wf_list[iw]);

  RefVector<const NLPPJob<RealType>> batch_list;
  std::vector<RealType> pairpots(nw);

  soecp_potential_list.reserve(nw);
  soecp_component_list.reserve(nw);
  pset_list.reserve(nw);
  psi_list.reserve(nw);
  batch_list.reserve(nw);

  for (size_t ig = 0; ig < pset_leader.groups(); ig++)
  {
    TrialWaveFunction::mw_prepareGroup(wf_list, p_list, ig);

    size_t max_num_jobs = 0;
    for (size_t iw = 0; iw < nw; iw++)
    {
      const auto& O = o_list.getCastedElement<SOECPotential>(iw);
      max_num_jobs  = std::max(max_num_jobs, O.sopp_jobs_[ig].size());
    }

    for (size_t jobid = 0; jobid < max_num_jobs; jobid++)
    {
      soecp_potential_list.clear();
      soecp_component_list.clear();
      pset_list.clear();
      psi_list.clear();
      batch_list.clear();
      for (size_t iw = 0; iw < nw; iw++)
      {
        auto& O = o_list.getCastedElement<SOECPotential>(iw);
        if (jobid < O.sopp_jobs_[ig].size())
        {
          const auto& job = O.sopp_jobs_[ig][jobid];
          soecp_potential_list.push_back(O);
          soecp_component_list.push_back(*O.PP_[job.ion_id]);
          pset_list.push_back(p_list[iw]);
          psi_list.push_back(wf_list[iw]);
          batch_list.push_back(job);
        }
      }

      SOECPComponent::mw_evaluateOne(soecp_component_list, pset_list, psi_list, batch_list, pairpots, O_leader.mw_res_->collection);

      for (size_t j = 0; j < soecp_potential_list.size(); j++)
        soecp_potential_list[j].get().value_ += pairpots[j];
    }
  }
}

std::unique_ptr<OperatorBase> SOECPotential::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  std::unique_ptr<SOECPotential> myclone = std::make_unique<SOECPotential>(IonConfig_, qp, psi);
  for (int ig = 0; ig < PPset_.size(); ++ig)
    if (PPset_[ig])
      myclone->addComponent(ig, std::unique_ptr<SOECPComponent>(PPset_[ig]->makeClone(qp)));
  return myclone;
}

void SOECPotential::addComponent(int groupID, std::unique_ptr<SOECPComponent>&& ppot)
{
  for (int iat = 0; iat < PP_.size(); iat++)
    if (IonConfig_.GroupID[iat] == groupID)
      PP_[iat] = ppot.get();
  PPset_[groupID] = std::move(ppot);
}

void SOECPotential::createResource(ResourceCollection& collection) const
{
  auto new_res = std::make_unique<SOECPotentialMultiWalkerResource>();
  for (int ig = 0; ig < PPset_.size(); ig++)
    if (PPset_[ig]->getVP())
    {
      PPset_[ig]->getVP()->createResource(new_res->collection);
      break;
    }
  auto resource_index = collection.addResource(std::move(new_res));
}

void SOECPotential::acquireResource(ResourceCollection& collection,
                                    const RefVectorWithLeader<OperatorBase>& o_list) const
{
  auto& O_leader = o_list.getCastedLeader<SOECPotential>();
  auto res_ptr   = dynamic_cast<SOECPotentialMultiWalkerResource*>(collection.lendResource().release());
  if (!res_ptr)
    throw std::runtime_error("SOECPotential::acquireResource dynamic_cast failed");
  O_leader.mw_res_.reset(res_ptr);
}

void SOECPotential::releaseResource(ResourceCollection& collection,
                                    const RefVectorWithLeader<OperatorBase>& o_list) const
{
  auto& O_leader = o_list.getCastedLeader<SOECPotential>();
  collection.takebackResource(std::move(O_leader.mw_res_));
}

Resource* SOECPotential::SOECPotentialMultiWalkerResource::makeClone() const
{
  return new SOECPotentialMultiWalkerResource(*this);
}

} // namespace qmcplusplus
