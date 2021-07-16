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
#include "ParticleIO/ParticleLayoutIO.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Particle/InitMolecularSystem.h"
#include "LongRange/LRCoulombSingleton.h"

namespace qmcplusplus
{
ParticleSetPool::ParticleSetPool(Communicate* c, const char* aname) : MPIObjectBase(c), TileMatrix(0)
{
  TileMatrix.diagonal(1);
  ClassName = "ParticleSetPool";
  myName    = aname;
}

ParticleSetPool::ParticleSetPool(ParticleSetPool&& other) noexcept
    : MPIObjectBase(other.myComm),
      SimulationCell(std::move(other.SimulationCell)),
      TileMatrix(other.TileMatrix),
      myPool(std::move(other.myPool))
{
  ClassName = other.ClassName;
  myName    = other.myName;
}

ParticleSetPool::~ParticleSetPool()
{
  PoolType::const_iterator it(myPool.begin()), it_end(myPool.end());
  while (it != it_end)
  {
    delete (*it).second;
    it++;
  }
}

ParticleSet* ParticleSetPool::getParticleSet(const std::string& pname)
{
  std::map<std::string, ParticleSet*>::iterator pit(myPool.find(pname));
  if (pit == myPool.end())
  {
    return 0;
  }
  else
  {
    return (*pit).second;
  }
}

MCWalkerConfiguration* ParticleSetPool::getWalkerSet(const std::string& pname)
{
  ParticleSet* mc = 0;
  if (myPool.size() == 1)
    mc = (*myPool.begin()).second;
  else
    mc = getParticleSet(pname);
  if (mc == 0)
  {
    APP_ABORT("ParticleSePool::getWalkerSet missing " + pname);
  }
  return dynamic_cast<MCWalkerConfiguration*>(mc);
}

void ParticleSetPool::addParticleSet(std::unique_ptr<ParticleSet>&& p)
{
  PoolType::iterator pit(myPool.find(p->getName()));
  if (pit == myPool.end())
  {
    auto& pname = p->getName();
    LOGMSG("  Adding " << pname << " ParticleSet to the pool")
    myPool[pname] = p.release();
  }
  else
  {
    WARNMSG("  " << p->getName() << " exists. Ignore addition")
  }
}

bool ParticleSetPool::putTileMatrix(xmlNodePtr cur)
{
  TileMatrix = 0;
  TileMatrix.diagonal(1);
  OhmmsAttributeSet pAttrib;
  pAttrib.add(TileMatrix, "tilematrix");
  pAttrib.put(cur);
  return true;
}

bool ParticleSetPool::putLattice(xmlNodePtr cur)
{
  ReportEngine PRE("ParticleSetPool", "putLattice");
  bool printcell = false;
  if (!SimulationCell)
  {
    app_debug() << "  Creating global supercell " << std::endl;
    SimulationCell = std::make_unique<ParticleSet::ParticleLayout_t>();
    printcell      = true;
  }
  else
  {
    app_log() << "  Overwriting global supercell " << std::endl;
  }
  LatticeParser a(*SimulationCell);
  bool lattice_defined = a.put(cur);
  if (printcell && lattice_defined)
  {
    if (outputManager.isHighActive())
    {
      SimulationCell->print(app_log(), 2);
    }
    else
    {
      SimulationCell->print(app_summary(), 1);
    }
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
  //const ParticleSet::ParticleLayout_t* sc=DistanceTable::getSimulationCell();
  //ParticleSet::ParticleLayout_t* sc=0;
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
#if defined(ENABLE_OFFLOAD)
  pAttrib.add(useGPU, "gpu", {"yes", "no"});
#endif
  pAttrib.put(cur);
  //backward compatibility
  if (id == "e" && role == "none")
    role = "MC";
  ParticleSet* pTemp = getParticleSet(id);
  if (pTemp == 0)
  {
    app_summary() << std::endl;
    app_summary() << " Particle Set" << std::endl;
    app_summary() << " ------------" << std::endl;
    app_summary() << "  Name: " << id << "   Offload : " << useGPU << std::endl;
    app_summary() << std::endl;

    // select OpenMP offload implementation in ParticleSet.
    if (useGPU == "yes")
      pTemp = new MCWalkerConfiguration(DynamicCoordinateKind::DC_POS_OFFLOAD);
    else
      pTemp = new MCWalkerConfiguration(DynamicCoordinateKind::DC_POS);
    //if(role == "MC")
    //  pTemp = new MCWalkerConfiguration;
    //else
    //  pTemp = new ParticleSet;
    if (SimulationCell)
    {
      app_log() << "  Initializing the lattice by the global supercell" << std::endl;
      pTemp->Lattice = *SimulationCell;
    }
    myPool[id] = pTemp;
    XMLParticleParser pread(*pTemp, TileMatrix);
    bool success = pread.put(cur);
    //if random_source is given, create a node <init target="" soruce=""/>
    if (randomR == "yes" && !randomsrc.empty())
    {
      xmlNodePtr anode = xmlNewNode(NULL, (const xmlChar*)"init");
      xmlNewProp(anode, (const xmlChar*)"source", (const xmlChar*)randomsrc.c_str());
      xmlNewProp(anode, (const xmlChar*)"target", (const xmlChar*)id.c_str());
      randomize_nodes.push_back(anode);
    }
    pTemp->setName(id);
    pTemp->is_spinor_ = spinor == "yes";
    app_summary() << "  Particle set size: " << pTemp->getTotalNum() << std::endl;
    app_summary() << std::endl;
    return success;
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
    APP_ABORT("ParticleSePool::randomize failed to randomize some Particlesets!");
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
  PoolType::iterator it(myPool.begin()), it_end(myPool.end());
  while (it != it_end)
  {
    ParticleSet* pt((*it).second);
    pt->update();
    ++it;
  }
}

} // namespace qmcplusplus
