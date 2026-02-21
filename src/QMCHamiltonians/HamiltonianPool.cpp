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
#include "HamiltonianFactory.h"
#include "QMCWaveFunctions/WaveFunctionPool.h"
#include "Particle/ParticleSetPool.h"
#include "OhmmsData/AttributeSet.h"
#include "Concurrency/OpenMP.h"
#include "Utilities/ProgressReportEngine.h"
#include "Message/UniformCommunicateError.h"

namespace qmcplusplus
{
HamiltonianPool::HamiltonianPool(ParticleSetPool& pset_pool,
                                 WaveFunctionPool& psi_pool,
                                 Communicate* c,
                                 const char* aname)
    : MPIObjectBase(c), ptcl_pool_(pset_pool), psi_pool_(psi_pool), curDoc(0)
{}

HamiltonianPool::~HamiltonianPool() = default;

bool HamiltonianPool::put(xmlNodePtr cur)
{
  ReportEngine PRE("HamiltonianPool", "put");
  std::string id("h0"), target("e"), role("extra"), psi_name;
  OhmmsAttributeSet hAttrib;
  hAttrib.add(id, "id");
  hAttrib.add(id, "name");
  hAttrib.add(role, "role");
  hAttrib.add(target, "target");
  hAttrib.add(psi_name, "wavefunction");
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

  if (myPool.find(id) != myPool.end())
    throw UniformCommunicateError(
        "Hamiltonian object named \"" + id +
        "\" already exists in the pool! Please set a different name using \"name\" attribute of \"hamiltonian\" node.");

  HamiltonianFactory ham_fac(id, *qp, ptcl_pool_.getPool(), psi_pool_.getWaveFunction(psi_name), myComm);
  ham_fac.put(cur);
  myPool.emplace(id, ham_fac.releaseHamiltonian());
  if (set2Primary)
    primaryH = myPool[id].get();
  return true;
}

bool HamiltonianPool::get(std::ostream& os) const
{
  for (auto& [name, ham] : myPool)
  {
    os << "  Hamiltonian " << name << std::endl;
    ham->get(os);
  }
  os << std::endl;
  return true;
}

void HamiltonianPool::reset() {}
} // namespace qmcplusplus
