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
#if OHMMS_DIM == 3
#include "ParticleIO/ESHDFParticleParser.h"
#endif
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

ParticleSetPool::ParticleSetPool(ParticleSetPool&& other)
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
    pTemp->IsSpinor = spinor == "yes";
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

/** Create particlesets from ES-HDF file
 */
ParticleSet* ParticleSetPool::createESParticleSet(xmlNodePtr cur, const std::string& target, ParticleSet* qp)
{
  //TinyVector<int,OHMMS_DIM> tilefactor;
  Tensor<int, OHMMS_DIM> eshdf_tilematrix(0);
  eshdf_tilematrix.diagonal(1);
  double lr_cut = 15;
  std::string h5name;
  std::string source("i");
  std::string bc("p p p");
  std::string spotype("0");
  std::string lr_handler("opt_breakup");
  OhmmsAttributeSet attribs;
  attribs.add(h5name, "href");
  attribs.add(eshdf_tilematrix, "tilematrix");
  attribs.add(source, "source");
  attribs.add(bc, "bconds");
  attribs.add(lr_cut, "LR_dim_cutoff");
  attribs.add(lr_handler, "LR_handler");
  attribs.add(spotype, "type");
  attribs.put(cur);

  if (h5name.empty())
    return qp;

#if OHMMS_DIM == 3
  ParticleSet* ions = getParticleSet(source);
  if (ions == 0)
  {
    ions = new MCWalkerConfiguration;
    ions->setName(source);
    //set the boundary condition
    ions->Lattice.LR_dim_cutoff = lr_cut;
    std::istringstream is(bc);
    char c;
    int idim = 0;
    while (!is.eof() && idim < OHMMS_DIM)
    {
      if (is >> c)
        ions->Lattice.BoxBConds[idim++] = (c == 'p');
    }

    tolower(lr_handler);
    if (lr_handler == "ewald")
    {
      LRCoulombSingleton::this_lr_type = LRCoulombSingleton::EWALD;
    }
    else if (lr_handler == "opt_breakup")
    {
      LRCoulombSingleton::this_lr_type = LRCoulombSingleton::ESLER;
    }
    else if (lr_handler == "opt_breakup_original")
    {
      LRCoulombSingleton::this_lr_type = LRCoulombSingleton::NATOLI;
    }
    else
    {
      APP_ABORT("Long range breakup handler not recognized\n");
    }

    //initialize ions from hdf5
    hid_t h5 = -1;
    if (myComm->rank() == 0)
    {
      //Rather than turn off all H5errors, we're going to
      //temporarily disable it.
      //
      //old_func:  function pointer to current function that
      //           displays when H5 encounters error.
      H5E_auto2_t old_func;
      //old_client_data:  null pointer to associated error stream.
      void* old_client_data;
      //Grab the current handler info.
      H5Eget_auto2(H5E_DEFAULT, &old_func, &old_client_data);
      //Now kill error notifications.
      H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
      h5 = H5Fopen(h5name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      //and restore to defaults.
      H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
      if (h5 < 0)
      {
        app_error() << "Could not open HDF5 file \"" << h5name
                    << "\" in ParticleSetPool::createESParticleSet.  Aborting.\n"
                    << "(Please ensure that your path is correct, the file exists, and that "
                    << "you have read permissions.)\n";
        APP_ABORT("ParticleSetPool::createESParticleSet");
      }
    }
    ESHDFIonsParser ap(*ions, h5, myComm);
    ap.put(cur);
    ap.expand(eshdf_tilematrix);
    if (h5 > -1)
      H5Fclose(h5);
    //failed to initialize the ions
    if (ions->getTotalNum() == 0)
      return 0;
    typedef ParticleSet::SingleParticleIndex_t SingleParticleIndex_t;
    std::vector<SingleParticleIndex_t> grid(OHMMS_DIM, SingleParticleIndex_t(1));
    ions->Lattice.reset();

    myPool[source] = ions;
  }

  if (SimulationCell == 0)
  {
    SimulationCell = std::make_unique<ParticleSet::ParticleLayout_t>(ions->Lattice);
  }

  if (qp == 0)
  {
    //create the electrons
    qp = new MCWalkerConfiguration;
    qp->setName(target);
    qp->Lattice = ions->Lattice;

    app_log() << "  Simulation cell radius = " << qp->Lattice.SimulationCellRadius << std::endl;
    app_log() << "  Wigner-Seitz cell radius = " << qp->Lattice.WignerSeitzRadius << std::endl;
    SimulationCell->print(app_log());

    // Goback to the // and OhmmsXPathObject handles it internally
    OhmmsXPathObject det("//determinant", cur);

    if (det.size() > 2)
      APP_ABORT("Only two electron groups are supported.");

    std::vector<int> num_spin(det.size(), 0);
    for (int i = 0; i < det.size(); ++i)
    {
      OhmmsAttributeSet a;
      a.add(num_spin[i], "size");
      a.put(det[i]);
    }

    {
      //create species
      SpeciesSet& species = qp->getSpeciesSet();
      //add up and down
      species.addSpecies("u");
      if (num_spin.size() > 1)
        species.addSpecies("d");
      int chid = species.addAttribute("charge");
      for (int i = 0; i < num_spin.size(); ++i)
        species(chid, i) = -1.0;
      int mid = species.addAttribute("membersize");
      for (int i = 0; i < num_spin.size(); ++i)
        species(mid, i) = num_spin[i];
      mid = species.addAttribute("mass");
      for (int i = 0; i < num_spin.size(); ++i)
        species(mid, i) = 1.0;
      qp->create(num_spin);
    }
    //name it with the target
    qp->setName(target);

    //if(qp->getTotalNum() == 0 || ions->getTotalNum() == 0)
    //{
    //  delete qp;
    //  delete ions;
    //  APP_ABORT("ParticleSetPool failed to create particlesets for the electron structure calculation");
    //  return 0;
    //}
    //for PPP, use uniform random
    if (qp->Lattice.SuperCellEnum == SUPERCELL_BULK)
    {
      makeUniformRandom(qp->R);
      qp->R.setUnit(PosUnit::Lattice);
      qp->convert2Cart(qp->R);
    }
    else
    {
      //assign non-trivial positions for the quanmtum particles
      InitMolecularSystem mole(*this);
      mole.initMolecule(ions, qp);
      qp->R.setUnit(PosUnit::Cartesian);
    }
    //for(int i=0; i<qp->getTotalNum(); ++i)
    //  std::cout << qp->GroupID[i] << " " << qp->R[i] << std::endl;

#if !defined(QMC_CUDA)
    makeUniformRandom(qp->spins);
    qp->spins *= 2 * M_PI;
#endif

    if (qp->Lattice.SuperCellEnum)
      qp->createSK();
    qp->resetGroups();
    myPool[target] = qp;
  }

#endif
  return qp;
}
} // namespace qmcplusplus
