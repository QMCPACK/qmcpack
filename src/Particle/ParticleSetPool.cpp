//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file ParticleSetPool.cpp
 * @brief Implements ParticleSetPool operators.
 */
#include "ParticleSetPool.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "ParticleIO/XMLParticleIO.h"
#include "ParticleIO/LatticeIO.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Particle/InitMolecularSystem.h"
#include "LongRange/LRCoulombSingleton.h"
#include <Message/UniformCommunicateError.h>
#include <PlatformSelector.hpp>

namespace qmcplusplus
{
ParticleSetPool::ParticleSetPool(Communicate* c, const char* aname)
    : MPIObjectBase(c), simulation_cell_(std::make_unique<SimulationCell>())
{
  ClassName = "ParticleSetPool";
  myName    = aname;
}

ParticleSetPool::ParticleSetPool(ParticleSetPool&& other) noexcept
    : MPIObjectBase(other.myComm), simulation_cell_(std::move(other.simulation_cell_)), myPool(std::move(other.myPool))
{
  ClassName = other.ClassName;
  myName    = other.myName;
}

ParticleSetPool::~ParticleSetPool() = default;

ParticleSet* ParticleSetPool::getParticleSet(const std::string& pname)
{
  if (auto pit = myPool.find(pname); pit == myPool.end())
    return nullptr;
  else
    return pit->second.get();
}

MCWalkerConfiguration* ParticleSetPool::getWalkerSet(const std::string& pname)
{
  ParticleSet* mc = 0;
  if (myPool.size() == 1)
    mc = myPool.begin()->second.get();
  else
    mc = getParticleSet(pname);
  if (mc == 0)
  {
    throw std::runtime_error("ParticleSePool::getWalkerSet missing " + pname);
  }
  return dynamic_cast<MCWalkerConfiguration*>(mc);
}

void ParticleSetPool::addParticleSet(std::unique_ptr<ParticleSet>&& p)
{
  const auto pit(myPool.find(p->getName()));
  if (pit == myPool.end())
  {
    auto& pname = p->getName();
    LOGMSG("  Adding " << pname << " ParticleSet to the pool")
    if (&p->getSimulationCell() != simulation_cell_.get())
      throw std::runtime_error("bug mandate");
    myPool.emplace(pname, std::move(p));
  }
  else
    throw std::runtime_error(p->getName() + " exists. Cannot be added again.");
}

bool ParticleSetPool::readSimulationCellXML(xmlNodePtr cur)
{
  ReportEngine PRE("ParticleSetPool", "putLattice");

  LatticeParser a(simulation_cell_->lattice_);
  bool lattice_defined = a.put(cur);
  if (lattice_defined)
  {
    app_log() << "  Overwriting global supercell " << std::endl;
    simulation_cell_->resetLRBox();
    if (outputManager.isHighActive())
      simulation_cell_->lattice_.print(app_log(), 2);
    else
      simulation_cell_->lattice_.print(app_summary(), 1);
  }
  return lattice_defined;
}

/** process an xml element
 * @param cur current xmlNodePtr
 * @return true, if successful.
 *
 * Creating MCWalkerConfiguration for all the ParticleSet
 * objects.
 */
bool ParticleSetPool::put(xmlNodePtr cur)
{
  ReportEngine PRE("ParticleSetPool", "put");
  std::string id("e");
  std::string role("none");
  std::string randomR("no");
  std::string randomsrc;
  std::string useGPU;
  std::string spinor;
  OhmmsAttributeSet pAttrib;
  pAttrib.add(id, "id");
  pAttrib.add(id, "name");
  pAttrib.add(role, "role");
  pAttrib.add(randomR, "random");
  pAttrib.add(randomsrc, "randomsrc");
  pAttrib.add(randomsrc, "random_source");
  pAttrib.add(spinor, "spinor", {"no", "yes"});
  pAttrib.add(useGPU, "gpu", CPUOMPTargetSelector::candidate_values);
  pAttrib.put(cur);
  //backward compatibility
  if (id == "e" && role == "none")
    role = "MC";
  ParticleSet* pTemp = getParticleSet(id);
  if (pTemp == 0)
  {
    const bool use_offload = CPUOMPTargetSelector::selectPlatform(useGPU) == PlatformKind::OMPTARGET;
    app_summary() << std::endl;
    app_summary() << " Particle Set" << std::endl;
    app_summary() << " ------------" << std::endl;
    app_summary() << "  Name: " << id << "   Offload : " << (use_offload ? "yes" : "no") << std::endl;
    app_summary() << std::endl;

    // select OpenMP offload implementation in ParticleSet.
    if (use_offload)
      pTemp = new MCWalkerConfiguration(*simulation_cell_, DynamicCoordinateKind::DC_POS_OFFLOAD);
    else
      pTemp = new MCWalkerConfiguration(*simulation_cell_, DynamicCoordinateKind::DC_POS);

    myPool.emplace(id, pTemp);

    try
    {
      XMLParticleParser pread(*pTemp);
      pread.readXML(cur);
    }
    catch (const UniformCommunicateError& ue)
    {
      myComm->barrier_and_abort(ue.what());
    }

    //if random_source is given, create a node <init target="" soruce=""/>
    if (randomR == "yes" && !randomsrc.empty())
    {
      xmlNodePtr anode = xmlNewNode(NULL, (const xmlChar*)"init");
      xmlNewProp(anode, (const xmlChar*)"source", (const xmlChar*)randomsrc.c_str());
      xmlNewProp(anode, (const xmlChar*)"target", (const xmlChar*)id.c_str());
      randomize_nodes.push_back(anode);
    }
    pTemp->setName(id);
    pTemp->setSpinor(spinor == "yes");
    app_summary() << "  Particle set size: " << pTemp->getTotalNum() << "   Groups : " << pTemp->groups() << std::endl;
    app_summary() << std::endl;
    return true;
  }
  else
  {
    app_warning() << "Particle set " << id << " is already created. Ignoring this section." << std::endl;
  }
  app_summary() << std::endl;
  return true;
}

void ParticleSetPool::randomize()
{
  app_log() << "ParticleSetPool::randomize " << randomize_nodes.size() << " ParticleSet"
            << (randomize_nodes.size() == 1 ? "" : "s") << "." << std::endl;
  bool success = true;
  for (int i = 0; i < randomize_nodes.size(); ++i)
  {
    InitMolecularSystem moinit(*this);
    success &= moinit.put(randomize_nodes[i]);
    xmlFreeNode(randomize_nodes[i]);
  }
  randomize_nodes.clear();
  if (!success)
    throw std::runtime_error("ParticleSePool::randomize failed to randomize some Particlesets!");
}

bool ParticleSetPool::get(std::ostream& os) const
{
  os << "ParticleSetPool has: " << std::endl << std::endl;
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.precision(14);
  PoolType::const_iterator it(myPool.begin()), it_end(myPool.end());
  while (it != it_end)
  {
    (*it).second->get(os);
    ++it;
  }
  return true;
}

void ParticleSetPool::output_particleset_info(Libxml2Document& doc, xmlNodePtr root)
{
  xmlNodePtr particles_info = doc.addChild(root, "particles");
  PoolType::const_iterator it(myPool.begin()), it_end(myPool.end());
  while (it != it_end)
  {
    xmlNodePtr particle = doc.addChild(particles_info, "particle");
    doc.addChild(particle, "name", (*it).second->getName());
    doc.addChild(particle, "size", (*it).second->getTotalNum());
    ++it;
  }
}

/** reset is used to initialize and evaluate the distance tables
 */
void ParticleSetPool::reset()
{
  for (const auto& [key, pset] : myPool)
    pset->update();
}

} // namespace qmcplusplus
