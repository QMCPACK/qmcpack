//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
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
#include "Particle/VirtualParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "Particle/createDistanceTable.h"
#include "QMCHamiltonians/NLPPJob.h"

namespace qmcplusplus
{
VirtualParticleSet::VirtualParticleSet(const ParticleSet& p, int nptcl) : refPS(p)
{
  setName("virtual");

  //initialize local data structure
  Lattice  = p.Lattice;
  TotalNum = nptcl;
  R.resize(nptcl);
  coordinates_->resize(nptcl);

  //create distancetables
  for (int i = 0; i < refPS.getNumDistTables(); ++i)
    addTable(refPS.getDistTable(i).origin());
}

/// move virtual particles to new postions and update distance tables
void VirtualParticleSet::makeMoves(int jel,
                                   const PosType& ref_pos,
                                   const std::vector<PosType>& deltaV,
                                   bool sphere,
                                   int iat)
{
  if (sphere && iat < 0)
    APP_ABORT("VirtualParticleSet::makeMoves is invoked incorrectly, the flag sphere=true requires iat specified!");
  onSphere      = sphere;
  refPtcl       = jel;
  refSourcePtcl = iat;
  assert(R.size() == deltaV.size());
  for (size_t ivp = 0; ivp < R.size(); ivp++)
    R[ivp] = ref_pos + deltaV[ivp];
  update();
}

void VirtualParticleSet::flex_makeMoves(const RefVector<VirtualParticleSet>& vp_list,
                                        const RefVector<const std::vector<PosType>>& deltaV_list,
                                        const RefVector<const NLPPJob<RealType>>& joblist,
                                        bool sphere)
{
  if (vp_list.size() > 1)
  {
    RefVector<ParticleSet> p_list;
    p_list.reserve(vp_list.size());

    for (int iw = 0; iw < vp_list.size(); iw++)
    {
      VirtualParticleSet& vp(vp_list[iw]);
      const std::vector<PosType>& deltaV(deltaV_list[iw]);
      const NLPPJob<RealType>& job(joblist[iw]);

      vp.onSphere      = sphere;
      vp.refPtcl       = job.electron_id;
      vp.refSourcePtcl = job.ion_id;
      assert(vp.R.size() == deltaV.size());
      for (size_t ivp = 0; ivp < vp.R.size(); ivp++)
        vp.R[ivp] = job.elec_pos + deltaV[ivp];
      p_list.push_back(vp);
    }

    ParticleSet::flex_update(p_list);
  }
  else if (vp_list.size() == 1)
    vp_list[0].get().makeMoves(joblist[0].get().electron_id, joblist[0].get().electron_id, deltaV_list[0].get(), sphere,
                               joblist[0].get().ion_id);
}

} // namespace qmcplusplus
