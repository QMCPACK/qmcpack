//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "SPOSetBuilderFactory.h"
#include "QMCWaveFunctions/SPOSetScanner.h"
#include "QMCWaveFunctions/ElectronGas/ElectronGasOrbitalBuilder.h"
#include "QMCWaveFunctions/HarmonicOscillator/SHOSetBuilder.h"
#if OHMMS_DIM == 3
#include "QMCWaveFunctions/LCAO/LCAOrbitalBuilder.h"

#if defined(QMC_COMPLEX)
#include "QMCWaveFunctions/EinsplineSpinorSetBuilder.h"
#include "QMCWaveFunctions/LCAO/LCAOSpinorBuilder.h"
#endif

#if defined(HAVE_EINSPLINE)
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#endif
#endif
#include "QMCWaveFunctions/CompositeSPOSet.h"
#include "Utilities/ProgressReportEngine.h"
#include "Utilities/IteratorUtility.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/MPIObjectBase.h"


namespace qmcplusplus
{

SPOSet* SPOSetBuilderFactory::getSPOSet(const std::string& name) const
{
  int nfound  = 0;
  SPOSet* spo = nullptr;
  for (auto it = spo_builders.begin(); it != spo_builders.end(); ++it)
  {
    std::vector<SPOSet*>& sposets = it->second->sposets;
    for (int i = 0; i < sposets.size(); ++i)
    {
      SPOSet* sposet = sposets[i];
      if (sposet->getName() == name)
      {
        spo = sposet;
        nfound++;
      }
    }
  }
  if (nfound > 1)
  {
    write_spo_builders();
    throw std::runtime_error("getSPOSet: requested sposet " + name + " is not unique!");
  }
  // keep this commented until legacy input styles are moved.
  // In legacy input styles, this look up may fail and need to build SPOSet on the fly.
  //else if (nfound == 0)
  //  throw std::runtime_error("getSPOSet: requested sposet " + name + " is not found!");
  return spo;
}


void SPOSetBuilderFactory::write_spo_builders(const std::string& pad) const
{
  std::string pad2 = pad + "  ";
  for (auto it = spo_builders.begin(); it != spo_builders.end(); ++it)
  {
    const std::string& type       = it->first;
    std::vector<SPOSet*>& sposets = it->second->sposets;
    app_log() << pad << "sposets for SPOSetBuilder of type " << type << std::endl;
    for (int i = 0; i < sposets.size(); ++i)
    {
      app_log() << pad2 << "sposet " << sposets[i]->getName() << std::endl;
    }
  }
}


/** constructor
 * \param els reference to the electrons
 * \param psi reference to the wavefunction
 * \param ions reference to the ions
 */
SPOSetBuilderFactory::SPOSetBuilderFactory(Communicate* comm, ParticleSet& els, PtclPoolType& psets)
    : MPIObjectBase(comm), last_builder(nullptr), targetPtcl(els), ptclPool(psets)
{
  ClassName = "SPOSetBuilderFactory";
}

SPOSetBuilderFactory::~SPOSetBuilderFactory() { DEBUG_MEMORY("SPOSetBuilderFactory::~SPOSetBuilderFactory"); }

SPOSetBuilder* SPOSetBuilderFactory::createSPOSetBuilder(xmlNodePtr rootNode)
{
  ReportEngine PRE(ClassName, "createSPOSetBuilder");
  std::string sourceOpt("ion0");
  std::string type("");
  std::string name("");
  std::string keyOpt("NMO");       //gaussian Molecular Orbital
  std::string transformOpt("yes"); //numerical Molecular Orbital
  std::string cuspC("no");         // cusp correction
  std::string cuspInfo("");        // file with precalculated cusp correction info
  std::string MOH5Ref("");         // Path to H5 file for MO calculations
  OhmmsAttributeSet aAttrib;
  aAttrib.add(sourceOpt, "source");
  aAttrib.add(cuspC, "cuspCorrection");
  aAttrib.add(type, "type");
  aAttrib.add(keyOpt, "keyword");
  aAttrib.add(keyOpt, "key");
  aAttrib.add(name, "name");
  aAttrib.add(transformOpt, "transform");
  aAttrib.add(cuspInfo, "cuspInfo");
  aAttrib.add(MOH5Ref, "href");

  std::cout << "IS SPINOR: " << targetPtcl.IsSpinor << std::endl;

  if (rootNode != NULL)
    aAttrib.put(rootNode);

  std::string type_in = type;
  tolower(type);

  //when name is missing, type becomes the input
  if (name.empty())
    name = type_in;

  SPOSetBuilder* bb = 0;

  //check if builder can be reused
  std::map<std::string, SPOSetBuilder*>::iterator bbit = spo_builders.find(name);
  if (bbit != spo_builders.end())
  {
    app_log() << "Reuse SPOSetBuilder \"" << name << "\" type " << type_in << std::endl;
    app_log().flush();
    bb                  = (*bbit).second;
    return last_builder = bb;
  }

  //assign last_builder
  bb = last_builder;

  if (type == "composite")
  {
    app_log() << "Composite SPO set with existing SPOSets." << std::endl;
    bb = new CompositeSPOSetBuilder(myComm, *this);
  }
  else if (type == "jellium" || type == "heg")
  {
    app_log() << "Electron gas SPO set" << std::endl;
    bb = new ElectronGasSPOBuilder(targetPtcl, myComm, rootNode);
  }
  else if (type == "sho")
  {
    app_log() << "Harmonic Oscillator SPO set" << std::endl;
    bb = new SHOSetBuilder(targetPtcl, myComm);
  }
#if OHMMS_DIM == 3
  else if (type.find("spline") < type.size())
  {
    if (targetPtcl.IsSpinor)
    {
#ifdef QMC_COMPLEX
      app_log() << "Einspline Spinor Set\n";
      bb = new EinsplineSpinorSetBuilder(targetPtcl, ptclPool, myComm, rootNode);
#else
      PRE.error("Use of einspline spinors requires QMC_COMPLEX=1.  Rebuild with this option");
#endif
    }
    else
    {
#if defined(HAVE_EINSPLINE)
      PRE << "EinsplineSetBuilder:  using libeinspline for B-spline orbitals.\n";
      bb = new EinsplineSetBuilder(targetPtcl, ptclPool, myComm, rootNode);
#else
      PRE.error("Einspline is missing for B-spline orbitals", true);
#endif
    }
  }
  else if (type == "molecularorbital" || type == "mo")
  {
      std::cout << "HERE " << std::endl;
    ParticleSet* ions = 0;
    //initialize with the source tag
    PtclPoolType::iterator pit(ptclPool.find(sourceOpt));
    if (pit == ptclPool.end())
      PRE.error("Missing basisset/@source.", true);
    else
      ions = (*pit).second;
    if (targetPtcl.IsSpinor)
#ifdef QMC_COMPLEX
      bb = new LCAOSpinorBuilder(targetPtcl, *ions, myComm, rootNode);
#else
      PRE.error("Use of lcao spinors requires QMC_COMPLEX=1.  Rebuild with this option");
#endif
    else
      bb = new LCAOrbitalBuilder(targetPtcl, *ions, myComm, rootNode);
  }
#endif //OHMMS_DIM==3
  PRE.flush();

  if (bb == 0)
    APP_ABORT_TRACE(__FILE__, __LINE__, "SPOSetBuilderFactory::createSPOSetBuilder\n  SPOSetBuilder creation failed.");

  if (bb == last_builder)
    app_log() << " Missing both \"@name\" and \"@type\". Use the last SPOSetBuilder." << std::endl;
  else
  {
    app_log() << "  Created SPOSet builder named '" << name << "' of type " << type << std::endl;
    spo_builders[name] = bb; //use name, if missing type is used
  }
  last_builder = bb;

  return bb;
}


SPOSet* SPOSetBuilderFactory::createSPOSet(xmlNodePtr cur)
{
  std::string sname("");
  std::string type("");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(sname, "name");
  aAttrib.add(type, "type");
  aAttrib.put(cur);

  SPOSetBuilder* bb = nullptr;
  if (type == "")
    bb = last_builder;
  else if (spo_builders.find(type) != spo_builders.end())
    bb = spo_builders[type];

  if (bb)
  {
    SPOSet* spo = bb->createSPOSet(cur);
    spo->setName(sname);
    return spo;
  }
  else
    throw std::runtime_error("Cannot find any SPOSetBuilder to build SPOSet!");
}

void SPOSetBuilderFactory::buildSPOSetCollection(xmlNodePtr cur)
{
  std::string collection_name;
  std::string collection_type;
  OhmmsAttributeSet attrib;
  attrib.add(collection_name, "name");
  attrib.add(collection_type, "type");
  attrib.put(cur);

  // use collection_type as collection_name if collection_name is not given
  if (collection_name.empty())
    collection_name = collection_type;

  app_summary() << std::endl;
  app_summary() << "   Single particle orbitals (SPO) collection" << std::endl;
  app_summary() << "   -----------------------------------------" << std::endl;
  app_summary() << "    Name: " << collection_name << "   Type input: " << collection_type << std::endl;
  app_summary() << std::endl;

  // create the SPOSet builder
  SPOSetBuilder* bb = createSPOSetBuilder(cur);

  // going through a list of sposet entries
  int nsposets = 0;
  processChildren(cur, [&](const std::string& cname, const xmlNodePtr element) {
    if (cname == "sposet")
    {
      SPOSet* spo = bb->createSPOSet(element);
      nsposets++;
    }
  });

  if (nsposets == 0)
    myComm->barrier_and_abort("SPOSetBuilderFactory::buildSPOSetCollection  no <sposet/> elements found");

  // going through a list of spo_scanner entries
  processChildren(cur, [&](const std::string& cname, const xmlNodePtr element) {
    if (cname == "spo_scanner")
      if (myComm->rank() == 0)
      {
        SPOSetScanner ascanner(bb->sposets, targetPtcl, ptclPool);
        ascanner.put(element);
      }
  });
}

std::string SPOSetBuilderFactory::basisset_tag = "basisset";

} // namespace qmcplusplus
