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

  Resource* makeClone() const override { return new VPMultiWalkerMem(*this); }
};

VirtualParticleSet::VirtualParticleSet(const ParticleSet& p, int nptcl) : ParticleSet(p.getSimulationCell()), refPS(p)
{
  setName("virtual");

  //initialize local data structure
  TotalNum = nptcl;
  R.resize(nptcl);
  coordinates_->resize(nptcl);

  //create distancetables
  for (int i = 0; i < refPS.getNumDistTables(); ++i)
    if (refPS.getDistTable(i).getModes() & DTModes::NEED_TEMP_DATA_ON_HOST)
      addTable(refPS.getDistTable(i).get_origin());
    else
      addTable(refPS.getDistTable(i).get_origin(), DTModes::MW_EVALUATE_RESULT_NO_TRANSFER_TO_HOST);
}

VirtualParticleSet::~VirtualParticleSet() = default;

Vector<int, OffloadPinnedAllocator<int>>& VirtualParticleSet::getMultiWalkerRefPctls()
{
  assert(mw_mem_ != nullptr);
  return mw_mem_->mw_refPctls;
}

const Vector<int, OffloadPinnedAllocator<int>>& VirtualParticleSet::getMultiWalkerRefPctls() const
{
  assert(mw_mem_ != nullptr);
  return mw_mem_->mw_refPctls;
}

void VirtualParticleSet::createResource(ResourceCollection& collection) const
{
  collection.addResource(std::make_unique<VPMultiWalkerMem>());

  ParticleSet::createResource(collection);
}

void VirtualParticleSet::acquireResource(ResourceCollection& collection,
                                         const RefVectorWithLeader<VirtualParticleSet>& vp_list)
{
  auto& vp_leader = vp_list.getLeader();
  auto res_ptr    = dynamic_cast<VPMultiWalkerMem*>(collection.lendResource().release());
  if (!res_ptr)
    throw std::runtime_error("VirtualParticleSet::acquireResource dynamic_cast failed");
  vp_leader.mw_mem_.reset(res_ptr);

  auto p_list = RefVectorWithLeaderParticleSet(vp_list);
  ParticleSet::acquireResource(collection, p_list);
}

void VirtualParticleSet::releaseResource(ResourceCollection& collection,
                                         const RefVectorWithLeader<VirtualParticleSet>& vp_list)
{
  collection.takebackResource(std::move(vp_list.getLeader().mw_mem_));
  auto p_list = RefVectorWithLeaderParticleSet(vp_list);
  ParticleSet::releaseResource(collection, p_list);
}

/// move virtual particles to new postions and update distance tables
void VirtualParticleSet::makeMoves(int jel,
                                   const PosType& ref_pos,
                                   const std::vector<PosType>& deltaV,
                                   bool sphere,
                                   int iat)
{
  if (sphere && iat < 0)
    throw std::runtime_error(
        "VirtualParticleSet::makeMoves is invoked incorrectly, the flag sphere=true requires iat specified!");
  onSphere      = sphere;
  refPtcl       = jel;
  refSourcePtcl = iat;
  assert(R.size() == deltaV.size());
  for (size_t ivp = 0; ivp < R.size(); ivp++)
    R[ivp] = ref_pos + deltaV[ivp];
  update();
}

void VirtualParticleSet::mw_makeMoves(const RefVectorWithLeader<VirtualParticleSet>& vp_list,
                                      const RefVector<const std::vector<PosType>>& deltaV_list,
                                      const RefVector<const NLPPJob<RealType>>& joblist,
                                      bool sphere)
{
  auto& vp_leader    = vp_list.getLeader();
  vp_leader.onSphere = sphere;

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
    vp.refPtcl       = job.electron_id;
    vp.refSourcePtcl = job.ion_id;
    assert(vp.R.size() == deltaV.size());
    for (size_t k = 0; k < vp.R.size(); k++, ivp++)
    {
      vp.R[k]          = job.elec_pos + deltaV[k];
      mw_refPctls[ivp] = vp.refPtcl;
    }
    p_list.push_back(vp);
  }
  assert(ivp == nVPs);

  mw_refPctls.updateTo();
  ParticleSet::mw_update(p_list);
}

} // namespace qmcplusplus
