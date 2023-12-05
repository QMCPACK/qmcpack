//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file HamiltonianPool.cpp
 * @brief Implements HamiltonianPool operators.
 */
#include "HamiltonianPool.h"
#include "QMCWaveFunctions/WaveFunctionPool.h"
#include "Particle/ParticleSetPool.h"
#include "OhmmsData/AttributeSet.h"
#include "Concurrency/OpenMP.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus
{
HamiltonianPool::HamiltonianPool(ParticleSetPool& pset_pool,
                                 WaveFunctionPool& psi_pool,
                                 Communicate* c,
                                 const char* aname)
    : MPIObjectBase(c), curH(0), ptcl_pool_(pset_pool), psi_pool_(psi_pool), curDoc(0)
{
  ClassName = "HamiltonianPool";
  myName    = aname;
}

HamiltonianPool::~HamiltonianPool()
{
  PoolType::iterator it(myPool.begin());
  while (it != myPool.end())
  {
    delete (*it).second;
    ++it;
  }
}

bool HamiltonianPool::put(xmlNodePtr cur)
{
  ReportEngine PRE("HamiltonianPool", "put");
  std::string id("h0"), target("e"), role("extra");
  OhmmsAttributeSet hAttrib;
  hAttrib.add(id, "id");
  hAttrib.add(id, "name");
  hAttrib.add(role, "role");
  hAttrib.add(target, "target");
  hAttrib.put(cur);
  ParticleSet* qp = ptcl_pool_.getParticleSet(target);
  if (qp == 0)
  {
    //never a good thing
    PRE.error("No target particle " + target + " exists.");
    return false;
  }
  bool set2Primary = false;
  //first Hamiltonian is set to the primary Hamiltonian
  if (myPool.empty() || role == "primary")
    set2Primary = true;
  HamiltonianFactory* curH = 0;
  PoolType::iterator hit(myPool.find(id));
  if (hit == myPool.end())
  {
    curH = new HamiltonianFactory(id, *qp, ptcl_pool_.getPool(), psi_pool_.getPool(), myComm);
    curH->setName(id);
    myPool[id] = curH;
  }
  else
    curH = (*hit).second;
  bool success = curH->put(cur);
  if (set2Primary)
    primaryH = curH->getH();
  return success;
}

bool HamiltonianPool::get(std::ostream& os) const
{
  for(auto& [name, factory] : myPool)
  {
    os << "  Hamiltonian " << name << std::endl;
    factory->getH()->get(os);
  }
  os << std::endl;
  return true;
}

void HamiltonianPool::reset() {}
} // namespace qmcplusplus
