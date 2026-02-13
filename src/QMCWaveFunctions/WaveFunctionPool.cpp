//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file WaveFunctionPool.cpp
 * @brief Implements WaveFunctionPool operators.
 */
#include <sstream>
#include "WaveFunctionPool.h"
#include "Message/UniformCommunicateError.h"
#include "Particle/ParticleSetPool.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
template class ObjectPool<TrialWaveFunction>;

WaveFunctionPool::WaveFunctionPool(const RuntimeOptions& runtime_options, ParticleSetPool& pset_pool, Communicate* c)
    : MPIObjectBase(c),
      ObjectPool("wavefunction"),
      runtime_options_(runtime_options),
      primary_psi_(nullptr),
      ptcl_pool_(pset_pool)
{}

WaveFunctionPool::~WaveFunctionPool() = default;

bool WaveFunctionPool::put(xmlNodePtr cur)
{
  std::string target("e"), role("extra");
  OhmmsAttributeSet pAttrib;
  pAttrib.add(target, "target");
  pAttrib.add(target, "ref");
  pAttrib.add(role, "role");
  pAttrib.put(cur);

  ParticleSet* qp = ptcl_pool_.getParticleSet(target);

  if (qp == nullptr)
    myComm->barrier_and_abort("target particle set named '" + target + "' not found");

  WaveFunctionFactory psiFactory(*qp, ptcl_pool_.getPool(), myComm);
  auto psi = psiFactory.buildTWF(cur, runtime_options_);
  addFactory(std::move(psi), empty() || role == "primary");
  return true;
}

void WaveFunctionPool::addFactory(std::unique_ptr<TrialWaveFunction> psi, bool primary)
{
  if (contains(psi->getName()))
    throw UniformCommunicateError("wavefunction " + psi->getName() + " exists. Cannot be added to the pool.");

  app_log() << "  Adding " << psi->getName() << " TrialWaveFunction to the pool" << std::endl;

  if (primary)
    primary_psi_ = psi.get();
  add(psi->getName(), std::move(psi));
}

TrialWaveFunction* WaveFunctionPool::getWaveFunction(const std::string& pname)
{
  try
  {
    return getObject(pname).get();
  }
  catch (const std::exception& re)
  {
    throw UniformCommunicateError(re.what());
  }
}

xmlNodePtr WaveFunctionPool::getWaveFunctionNode(const std::string& id)
{
  return getWaveFunction(id)->getNode();
}
} // namespace qmcplusplus
