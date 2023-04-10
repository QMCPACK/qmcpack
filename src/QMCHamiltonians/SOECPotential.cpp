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
#include <optional>

namespace qmcplusplus
{

struct SOECPotential::SOECPotentialMultiWalkerResource : public Resource
{
  SOECPotentialMultiWalkerResource() : Resource("SOECPotential") {}

  std::unique_ptr<Resource> makeClone() const override
  {
    return std::make_unique<SOECPotentialMultiWalkerResource>(*this);
  }


  ResourceCollection collection{"SOPPcollection"};

  //crowds worth of per particle soecp values
  Matrix<Real> ve_samples;
  Matrix<Real> vi_samples;
};

/** constructor
 *\param ionic positions
 *\param els electronic poitions
 *\param psi Trial wave function
*/
SOECPotential::SOECPotential(ParticleSet& ions, ParticleSet& els, TrialWaveFunction& psi)
    : my_rng_(nullptr), ion_config_(ions), psi_(psi), peln_(els), elec_neighbor_ions_(els), ion_neighbor_elecs_(ions)
{
  setEnergyDomain(POTENTIAL);
  twoBodyQuantumDomain(ions, els);
  my_table_index_ = els.addTable(ions);
  num_ions_       = ions.getTotalNum();
  pp_.resize(num_ions_, nullptr);
  ppset_.resize(ion_config_.getSpeciesSet().getTotalNum());
  sopp_jobs_.resize(els.groups());
  for (size_t ig = 0; ig < els.groups(); ig++)
    sopp_jobs_[ig].reserve(2 * els.groupsize(ig));
}

SOECPotential::~SOECPotential() = default;

void SOECPotential::resetTargetParticleSet(ParticleSet& P) {}

SOECPotential::Return_t SOECPotential::evaluate(ParticleSet& P)
{
  evaluateImpl(P, false);
  return value_;
}

SOECPotential::Return_t SOECPotential::evaluateDeterministic(ParticleSet& P)
{
  evaluateImpl(P, true);
  return value_;
}

void SOECPotential::evaluateImpl(ParticleSet& P, bool keep_grid)
{
  value_ = 0.0;
  if (!keep_grid)
    for (int ipp = 0; ipp < ppset_.size(); ipp++)
      if (ppset_[ipp])
        ppset_[ipp]->rotateQuadratureGrid(generateRandomRotationMatrix(*my_rng_));

  const auto& ble = P.getDistTableAB(my_table_index_);
  for (int iat = 0; iat < num_ions_; iat++)
    ion_neighbor_elecs_.getNeighborList(iat).clear();
  for (int jel = 0; jel < P.getTotalNum(); jel++)
    elec_neighbor_ions_.getNeighborList(jel).clear();

  for (int jel = 0; jel < P.getTotalNum(); jel++)
  {
    const auto& dist               = ble.getDistRow(jel);
    const auto& displ              = ble.getDisplRow(jel);
    std::vector<int>& NeighborIons = elec_neighbor_ions_.getNeighborList(jel);
    for (int iat = 0; iat < num_ions_; iat++)
      if (pp_[iat] != nullptr && dist[iat] < pp_[iat]->getRmax())
      {
        RealType pairpot = pp_[iat]->evaluateOne(P, iat, psi_, jel, dist[iat], -displ[iat]);
        value_ += pairpot;
        NeighborIons.push_back(iat);
        ion_neighbor_elecs_.getNeighborList(iat).push_back(jel);
      }
  }
}

SOECPotential::Return_t SOECPotential::evaluateValueAndDerivatives(ParticleSet& P,
                                                                   const opt_variables_type& optvars,
                                                                   const Vector<ValueType>& dlogpsi,
                                                                   Vector<ValueType>& dhpsioverpsi)
{
  value_ = 0.0;
  for (int ipp = 0; ipp < ppset_.size(); ipp++)
    if (ppset_[ipp])
      ppset_[ipp]->rotateQuadratureGrid(generateRandomRotationMatrix(*my_rng_));

  const auto& ble = P.getDistTableAB(my_table_index_);
  for (int jel = 0; jel < P.getTotalNum(); jel++)
  {
    const auto& dist  = ble.getDistRow(jel);
    const auto& displ = ble.getDisplRow(jel);
    for (int iat = 0; iat < num_ions_; iat++)
      if (pp_[iat] != nullptr && dist[iat] < pp_[iat]->getRmax())
        value_ += pp_[iat]->evaluateValueAndDerivatives(P, iat, psi_, jel, dist[iat], -displ[iat], optvars, dlogpsi,
                                                        dhpsioverpsi);
  }
  return value_;
}

void SOECPotential::mw_evaluate(const RefVectorWithLeader<OperatorBase>& o_list,
                                const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                const RefVectorWithLeader<ParticleSet>& p_list) const
{
  mw_evaluateImpl(o_list, wf_list, p_list, std::nullopt);
}

void SOECPotential::mw_evaluatePerParticle(const RefVectorWithLeader<OperatorBase>& o_list,
                                           const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                           const RefVectorWithLeader<ParticleSet>& p_list,
                                           const std::vector<ListenerVector<Real>>& listeners,
                                           const std::vector<ListenerVector<Real>>& listeners_ions) const
{
  std::optional<ListenerOption<Real>> l_opt(std::in_place, listeners, listeners_ions);
  mw_evaluateImpl(o_list, wf_list, p_list, l_opt);
}

void SOECPotential::mw_evaluateImpl(const RefVectorWithLeader<OperatorBase>& o_list,
                                    const RefVectorWithLeader<TrialWaveFunction>& wf_list,
                                    const RefVectorWithLeader<ParticleSet>& p_list,
                                    const std::optional<ListenerOption<Real>> listeners,
                                    bool keep_grid)
{
  auto& O_leader           = o_list.getCastedLeader<SOECPotential>();
  ParticleSet& pset_leader = p_list.getLeader();
  const size_t nw          = o_list.size();

  for (size_t iw = 0; iw < nw; iw++)
  {
    auto& O = o_list.getCastedElement<SOECPotential>(iw);
    const ParticleSet& P(p_list[iw]);

    if (!keep_grid)
      for (size_t ipp = 0; ipp < O.ppset_.size(); ipp++)
        if (O.ppset_[ipp])
          O.ppset_[ipp]->rotateQuadratureGrid(generateRandomRotationMatrix(*O.my_rng_));

    //loop over all the ions
    const auto& ble = P.getDistTableAB(O.my_table_index_);
    //clear elec and ion neighbor lists
    for (size_t iat = 0; iat < O.num_ions_; iat++)
      O.ion_neighbor_elecs_.getNeighborList(iat).clear();
    for (size_t jel = 0; jel < P.getTotalNum(); jel++)
      O.elec_neighbor_ions_.getNeighborList(jel).clear();

    for (size_t ig = 0; ig < P.groups(); ig++) // loop over species
    {
      auto& joblist = O.sopp_jobs_[ig];
      joblist.clear();

      for (size_t jel = P.first(ig); jel < P.last(ig); jel++)
      {
        const auto& dist               = ble.getDistRow(jel);
        const auto& displ              = ble.getDisplRow(jel);
        std::vector<int>& NeighborIons = O.elec_neighbor_ions_.getNeighborList(jel);
        for (size_t iat = 0; iat < O.num_ions_; iat++)
          if (O.pp_[iat] != nullptr && dist[iat] < O.pp_[iat]->getRmax())
          {
            NeighborIons.push_back(iat);
            O.ion_neighbor_elecs_.getNeighborList(iat).push_back(jel);
            joblist.emplace_back(iat, jel, dist[iat], -displ[iat]);
          }
      }
    }
    O.value_ = 0.0;
  }

  if (listeners)
  {
    auto& ve_samples = O_leader.mw_res_handle_.getResource().ve_samples;
    auto& vi_samples = O_leader.mw_res_handle_.getResource().vi_samples;
    ve_samples.resize(nw, pset_leader.getTotalNum());
    vi_samples.resize(nw, O_leader.ion_config_.getTotalNum());
  }

  auto pp_component = std::find_if(O_leader.ppset_.begin(), O_leader.ppset_.end(), [](auto& ptr) { return bool(ptr); });
  assert(pp_component != std::end(O_leader.ppset_));

  RefVector<SOECPotential> soecp_potential_list;
  RefVectorWithLeader<SOECPComponent> soecp_component_list(**pp_component);
  RefVectorWithLeader<ParticleSet> pset_list(pset_leader);
  RefVectorWithLeader<TrialWaveFunction> psi_list(O_leader.psi_);
  assert(&O_leader.psi_ == &wf_list.getLeader());
  for (size_t iw = 0; iw < nw; iw++)
    assert(&o_list.getCastedElement<SOECPotential>(iw).psi_ == &wf_list[iw]);

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
          soecp_component_list.push_back(*O.pp_[job.ion_id]);
          pset_list.push_back(p_list[iw]);
          psi_list.push_back(wf_list[iw]);
          batch_list.push_back(job);
        }
      }

      SOECPComponent::mw_evaluateOne(soecp_component_list, pset_list, psi_list, batch_list, pairpots,
                                     O_leader.mw_res_handle_.getResource().collection);

      for (size_t j = 0; j < soecp_potential_list.size(); j++)
      {
        soecp_potential_list[j].get().value_ += pairpots[j];

        if (listeners)
        {
          auto& ve_samples = O_leader.mw_res_handle_.getResource().ve_samples;
          auto& vi_samples = O_leader.mw_res_handle_.getResource().vi_samples;
          int iw           = j;
          ve_samples(iw, batch_list[j].get().electron_id) += pairpots[j];
          vi_samples(iw, batch_list[j].get().ion_id) += pairpots[j];
        }
      }
    }
  }

  if (listeners)
  {
    auto& ve_samples  = O_leader.mw_res_handle_.getResource().ve_samples;
    auto& vi_samples  = O_leader.mw_res_handle_.getResource().vi_samples;
    int num_electrons = pset_leader.getTotalNum();
    for (int iw = 0; iw < nw; iw++)
    {
      Vector<Real> ve_sample(ve_samples.begin(iw), num_electrons);
      Vector<Real> vi_sample(vi_samples.begin(iw), O_leader.num_ions_);
      for (const ListenerVector<Real>& listener : listeners->electron_values)
        listener.report(iw, O_leader.getName(), ve_sample);
      for (const ListenerVector<Real>& listener : listeners->ion_values)
        listener.report(iw, O_leader.getName(), vi_sample);
    }
    ve_samples = 0.0;
    vi_samples = 0.0;
  }
}

std::unique_ptr<OperatorBase> SOECPotential::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  std::unique_ptr<SOECPotential> myclone = std::make_unique<SOECPotential>(ion_config_, qp, psi);
  for (int ig = 0; ig < ppset_.size(); ++ig)
    if (ppset_[ig])
      myclone->addComponent(ig, std::unique_ptr<SOECPComponent>(ppset_[ig]->makeClone(qp)));
  return myclone;
}

void SOECPotential::addComponent(int groupID, std::unique_ptr<SOECPComponent>&& ppot)
{
  for (int iat = 0; iat < pp_.size(); iat++)
    if (ion_config_.GroupID[iat] == groupID)
      pp_[iat] = ppot.get();
  ppset_[groupID] = std::move(ppot);
}

void SOECPotential::createResource(ResourceCollection& collection) const
{
  auto new_res = std::make_unique<SOECPotentialMultiWalkerResource>();
  for (int ig = 0; ig < ppset_.size(); ig++)
    if (ppset_[ig]->getVP())
    {
      ppset_[ig]->getVP()->createResource(new_res->collection);
      break;
    }
  auto resource_index = collection.addResource(std::move(new_res));
}

void SOECPotential::acquireResource(ResourceCollection& collection,
                                    const RefVectorWithLeader<OperatorBase>& o_list) const
{
  auto& O_leader          = o_list.getCastedLeader<SOECPotential>();
  O_leader.mw_res_handle_ = collection.lendResource<SOECPotentialMultiWalkerResource>();
}

void SOECPotential::releaseResource(ResourceCollection& collection,
                                    const RefVectorWithLeader<OperatorBase>& o_list) const
{
  auto& O_leader = o_list.getCastedLeader<SOECPotential>();
  collection.takebackResource(O_leader.mw_res_handle_);
}

} // namespace qmcplusplus
