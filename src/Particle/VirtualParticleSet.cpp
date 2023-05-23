//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file VirtualParticleSet.cpp
 * A proxy class to the quantum ParticleSet
 */

#include "Configuration.h"
#include "VirtualParticleSet.h"
#include "Particle/DistanceTable.h"
#include "Particle/createDistanceTable.h"
#include "QMCHamiltonians/NLPPJob.h"
#include "ResourceCollection.h"

namespace qmcplusplus
{

struct VPMultiWalkerMem : public Resource
{
  /// multi walker reference particle
  Vector<int, OffloadPinnedAllocator<int>> mw_refPctls;

  VPMultiWalkerMem() : Resource("VPMultiWalkerMem") {}

  VPMultiWalkerMem(const VPMultiWalkerMem&) : VPMultiWalkerMem() {}

  std::unique_ptr<Resource> makeClone() const override { return std::make_unique<VPMultiWalkerMem>(*this); }
};

VirtualParticleSet::VirtualParticleSet(const ParticleSet& p, int nptcl, size_t dt_count_limit)
    : ParticleSet(p.getSimulationCell())
{
  setName("virtual");

  //initialize local data structure
  setSpinor(p.isSpinor());
  TotalNum = nptcl;
  R.resize(nptcl);
  if (isSpinor())
    spins.resize(nptcl);
  coordinates_->resize(nptcl);

  //create distancetables
  assert(dt_count_limit <= p.getNumDistTables());
  if (dt_count_limit == 0)
    dt_count_limit = p.getNumDistTables();
  for (int i = 0; i < dt_count_limit; ++i)
    if (p.getDistTable(i).getModes() & DTModes::NEED_VP_FULL_TABLE_ON_HOST)
      addTable(p.getDistTable(i).get_origin());
    else
      addTable(p.getDistTable(i).get_origin(), DTModes::MW_EVALUATE_RESULT_NO_TRANSFER_TO_HOST);
}

VirtualParticleSet::~VirtualParticleSet() = default;

Vector<int, OffloadPinnedAllocator<int>>& VirtualParticleSet::getMultiWalkerRefPctls()
{
  return mw_mem_handle_.getResource().mw_refPctls;
}

const Vector<int, OffloadPinnedAllocator<int>>& VirtualParticleSet::getMultiWalkerRefPctls() const
{
  return mw_mem_handle_.getResource().mw_refPctls;
}

void VirtualParticleSet::createResource(ResourceCollection& collection) const
{
  collection.addResource(std::make_unique<VPMultiWalkerMem>());
  ParticleSet::createResource(collection);
}

void VirtualParticleSet::acquireResource(ResourceCollection& collection,
                                         const RefVectorWithLeader<VirtualParticleSet>& vp_list)
{
  auto& vp_leader          = vp_list.getLeader();
  vp_leader.mw_mem_handle_ = collection.lendResource<VPMultiWalkerMem>();

  auto p_list = RefVectorWithLeaderParticleSet(vp_list);
  ParticleSet::acquireResource(collection, p_list);
}

void VirtualParticleSet::releaseResource(ResourceCollection& collection,
                                         const RefVectorWithLeader<VirtualParticleSet>& vp_list)
{
  collection.takebackResource(vp_list.getLeader().mw_mem_handle_);
  auto p_list = RefVectorWithLeaderParticleSet(vp_list);
  ParticleSet::releaseResource(collection, p_list);
}

/// move virtual particles to new postions and update distance tables
void VirtualParticleSet::makeMoves(const ParticleSet& refp,
                                   int jel,
                                   const std::vector<PosType>& deltaV,
                                   bool sphere,
                                   int iat)
{
  if (sphere && iat < 0)
    throw std::runtime_error(
        "VirtualParticleSet::makeMoves is invoked incorrectly, the flag sphere=true requires iat specified!");
  onSphere      = sphere;
  refPS         = refp;
  refPtcl       = jel;
  refSourcePtcl = iat;
  assert(R.size() == deltaV.size());
  for (size_t ivp = 0; ivp < R.size(); ivp++)
    R[ivp] = refp.R[jel] + deltaV[ivp];
  if (refp.isSpinor())
    for (size_t ivp = 0; ivp < R.size(); ivp++)
      spins[ivp] = refp.spins[jel]; //no spin deltas in this API
  update();
}

/// move virtual particles to new postions and update distance tables
void VirtualParticleSet::makeMovesWithSpin(const ParticleSet& refp,
                                           int jel,
                                           const std::vector<PosType>& deltaV,
                                           const std::vector<RealType>& deltaS,
                                           bool sphere,
                                           int iat)
{
  assert(refp.isSpinor());
  if (sphere && iat < 0)
    throw std::runtime_error(
        "VirtualParticleSet::makeMovesWithSpin is invoked incorrectly, the flag sphere=true requires iat specified!");
  onSphere      = sphere;
  refPS         = refp;
  refPtcl       = jel;
  refSourcePtcl = iat;
  assert(R.size() == deltaV.size());
  assert(spins.size() == deltaS.size());
  for (size_t ivp = 0; ivp < R.size(); ivp++)
  {
    R[ivp]     = refp.R[jel] + deltaV[ivp];
    spins[ivp] = refp.spins[jel] + deltaS[ivp];
  }
  update();
}

void VirtualParticleSet::mw_makeMoves(const RefVectorWithLeader<VirtualParticleSet>& vp_list,
                                      const RefVectorWithLeader<ParticleSet>& refp_list,
                                      const RefVector<const std::vector<PosType>>& deltaV_list,
                                      const RefVector<const NLPPJob<RealType>>& joblist,
                                      bool sphere)
{
  auto& vp_leader    = vp_list.getLeader();
  vp_leader.onSphere = sphere;
  vp_leader.refPS    = refp_list.getLeader();

  const size_t nVPs = countVPs(vp_list);
  auto& mw_refPctls = vp_leader.getMultiWalkerRefPctls();
  mw_refPctls.resize(nVPs);

  RefVectorWithLeader<ParticleSet> p_list(vp_leader);
  p_list.reserve(vp_list.size());

  size_t ivp = 0;
  for (int iw = 0; iw < vp_list.size(); iw++)
  {
    VirtualParticleSet& vp(vp_list[iw]);
    const std::vector<PosType>& deltaV(deltaV_list[iw]);
    const NLPPJob<RealType>& job(joblist[iw]);

    vp.onSphere      = sphere;
    vp.refPS         = refp_list[iw];
    vp.refPtcl       = job.electron_id;
    vp.refSourcePtcl = job.ion_id;
    assert(vp.R.size() == deltaV.size());
    for (size_t k = 0; k < vp.R.size(); k++, ivp++)
    {
      vp.R[k] = refp_list[iw].R[vp.refPtcl] + deltaV[k];
      if (vp_leader.isSpinor())
        vp.spins[k] = refp_list[iw].spins[vp.refPtcl]; //no spin deltas in this API
      mw_refPctls[ivp] = vp.refPtcl;
    }
    p_list.push_back(vp);
  }
  assert(ivp == nVPs);

  mw_refPctls.updateTo();
  ParticleSet::mw_update(p_list);
}

void VirtualParticleSet::mw_makeMovesWithSpin(const RefVectorWithLeader<VirtualParticleSet>& vp_list,
                                              const RefVectorWithLeader<ParticleSet>& refp_list,
                                              const RefVector<const std::vector<PosType>>& deltaV_list,
                                              const RefVector<const std::vector<RealType>>& deltaS_list,
                                              const RefVector<const NLPPJob<RealType>>& joblist,
                                              bool sphere)
{
  auto& vp_leader = vp_list.getLeader();
  if (!vp_leader.isSpinor())
    throw std::runtime_error(
        "VirtualParticleSet::mw_makeMovesWithSpin should not be called if particle sets aren't spionor types");
  vp_leader.onSphere = sphere;
  vp_leader.refPS    = refp_list.getLeader();

  const size_t nVPs = countVPs(vp_list);
  auto& mw_refPctls = vp_leader.getMultiWalkerRefPctls();
  mw_refPctls.resize(nVPs);

  RefVectorWithLeader<ParticleSet> p_list(vp_leader);
  p_list.reserve(vp_list.size());

  size_t ivp = 0;
  for (int iw = 0; iw < vp_list.size(); iw++)
  {
    VirtualParticleSet& vp(vp_list[iw]);
    const std::vector<PosType>& deltaV(deltaV_list[iw]);
    const std::vector<RealType>& deltaS(deltaS_list[iw]);
    const NLPPJob<RealType>& job(joblist[iw]);

    vp.onSphere      = sphere;
    vp.refPS         = refp_list[iw];
    vp.refPtcl       = job.electron_id;
    vp.refSourcePtcl = job.ion_id;
    assert(vp.R.size() == deltaV.size());
    assert(vp.spins.size() == deltaS.size());
    assert(vp.R.size() == vp.spins.size());
    for (size_t k = 0; k < vp.R.size(); k++, ivp++)
    {
      vp.R[k]          = refp_list[iw].R[vp.refPtcl] + deltaV[k];
      vp.spins[k]      = refp_list[iw].spins[vp.refPtcl] + deltaS[k];
      mw_refPctls[ivp] = vp.refPtcl;
    }
    p_list.push_back(vp);
  }
  assert(ivp == nVPs);

  mw_refPctls.updateTo();
  ParticleSet::mw_update(p_list);
}

} // namespace qmcplusplus
