//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
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


#include "QMCWaveFunctions/SPOSetBuilderFactory.h"
#include "QMCWaveFunctions/ElectronGas/ElectronGasOrbitalBuilder.h"
#include "QMCWaveFunctions/HarmonicOscillator/SHOSetBuilder.h"
#if OHMMS_DIM == 3
#if defined(ENABLE_SOA)
#include "QMCWaveFunctions/lcao/LCAOrbitalBuilder.h"
#else
#if !defined(QMC_COMPLEX)
#include "QMCWaveFunctions/MolecularOrbitals/NGOBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/GTOBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/STOBuilder.h"
#include "QMCWaveFunctions/MolecularOrbitals/MolecularSPOBuilder.h"
#endif
#endif

#if !defined(QMC_COMPLEX)
#include "QMCWaveFunctions/RotatedSPOs.h"
#endif

#if defined(HAVE_EINSPLINE)
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/EinsplineSpinorSetBuilder.h"
#endif
#endif
#include "QMCWaveFunctions/CompositeSPOSet.h"
#include "Utilities/ProgressReportEngine.h"
#include "Utilities/IteratorUtility.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/MPIObjectBase.h"


namespace qmcplusplus
{
//initialization of the static data of SPOSetBuilderFactory
std::map<std::string, SPOSetBuilder*> SPOSetBuilderFactory::spo_builders;
SPOSetBuilder* SPOSetBuilderFactory::last_builder = 0;

void SPOSetBuilderFactory::clear()
{
  spo_builders.clear();
  last_builder = nullptr;
}

SPOSet* get_sposet(const std::string& name)
{
  int nfound  = 0;
  SPOSet* spo = 0;
  std::map<std::string, SPOSetBuilder*>::iterator it;
  for (it = SPOSetBuilderFactory::spo_builders.begin(); it != SPOSetBuilderFactory::spo_builders.end(); ++it)
  {
    std::vector<SPOSet*>& sposets = it->second->sposets;
    for (int i = 0; i < sposets.size(); ++i)
    {
      SPOSet* sposet = sposets[i];
      if (sposet->objectName == name)
      {
        spo = sposet;
        nfound++;
      }
    }
  }
  if (nfound > 1)
  {
    write_spo_builders();
    APP_ABORT_TRACE(__FILE__, __LINE__, "get_sposet: requested sposet " + name + " is not unique");
  }
  //else if(spo==NULL)
  //{
  //  write_spo_builders();
  //  APP_ABORT("get_sposet: requested sposet "+name+" does not exist");
  //}
  return spo;
}


void write_spo_builders(const std::string& pad)
{
  std::string pad2 = pad + "  ";
  std::map<std::string, SPOSetBuilder*>::iterator it;
  for (it = SPOSetBuilderFactory::spo_builders.begin(); it != SPOSetBuilderFactory::spo_builders.end(); ++it)
  {
    const std::string& type       = it->first;
    std::vector<SPOSet*>& sposets = it->second->sposets;
    app_log() << pad << "sposets for SPOSetBuilder of type " << type << std::endl;
    for (int i = 0; i < sposets.size(); ++i)
    {
      app_log() << pad2 << "sposet " << sposets[i]->objectName << std::endl;
    }
  }
}


/** constructor
 * \param els reference to the electrons
 * \param psi reference to the wavefunction
 * \param ions reference to the ions
 */
SPOSetBuilderFactory::SPOSetBuilderFactory(Communicate* comm, ParticleSet& els, PtclPoolType& psets)
    : MPIObjectBase(comm), targetPtcl(els), ptclPool(psets)
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
    bb = new CompositeSPOSetBuilder(myComm);
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
    name = type_in;
#if defined(HAVE_EINSPLINE)
    PRE << "EinsplineSetBuilder:  using libeinspline for B-spline orbitals.\n";
   // bb = new EinsplineSetBuilder(targetPtcl, ptclPool, myComm, rootNode);
    bb = new EinsplineSpinorSetBuilder(targetPtcl, ptclPool, myComm, rootNode);
#else
    PRE.error("Einspline is missing for B-spline orbitals", true);
#endif
  }
  else if (type == "molecularorbital" || type == "mo")
  {
    ParticleSet* ions = 0;
    //initialize with the source tag
    PtclPoolType::iterator pit(ptclPool.find(sourceOpt));
    if (pit == ptclPool.end())
      PRE.error("Missing basisset/@source.", true);
    else
      ions = (*pit).second;
#if defined(ENABLE_SOA)
    bb = new LCAOrbitalBuilder(targetPtcl, *ions, myComm, rootNode);
#else
#if defined(QMC_COMPLEX)
    PRE.error("Complex molecular orbitals not supported by AoS builds!", true);
#else
    if (transformOpt == "yes")
    {
      app_log() << "Using MolecularSPOBuilder<NGOBuilder>" << std::endl;
      bb = new MolecularSPOBuilder<NGOBuilder>(targetPtcl, *ions, myComm, cuspC == "yes", cuspInfo, MOH5Ref);
    }
    else
    {
      if (cuspC == "yes")
        app_log() << " ****** Cusp Correction algorithm is only implemented in combination with numerical radial "
                     "orbitals. Use transform=yes to enable this option. \n";
      if (keyOpt == "GTO")
        bb = new MolecularSPOBuilder<GTOBuilder>(targetPtcl, *ions, myComm);
      else if (keyOpt == "STO")
        bb = new MolecularSPOBuilder<STOBuilder>(targetPtcl, *ions, myComm);
    }
#endif //QMC_COMPLEX
#endif
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
  std::string bname("");
  std::string sname("");
  std::string type("");
  std::string rotation("no");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(bname, "basisset");
  aAttrib.add(sname, "name");
  aAttrib.add(type, "type");
  aAttrib.add(rotation, "optimize");
  //aAttrib.put(rcur);
  aAttrib.put(cur);

  //tolower(type);

  SPOSetBuilder* bb;
  if (bname == "")
    bname = type;
  if (type == "")
    bb = last_builder;
  else if (spo_builders.find(type) != spo_builders.end())
    bb = spo_builders[type];
  else
  {
    std::string cname("");
    xmlNodePtr tcur = cur->children;
    if (tcur != NULL)
      getNodeName(cname, tcur);
    if (cname == basisset_tag)
    {
      bb = createSPOSetBuilder(cur);
      bb->loadBasisSetFromXML(tcur);
    }
    else
      bb = createSPOSetBuilder(cur);
  }
  if (bb)
  {
    app_log() << "  Building SPOSet '" << sname << "' with '" << bname << "' basis set." << std::endl;
    SPOSet* spo     = bb->createSPOSet(cur);
    spo->objectName = sname;
    if (rotation == "yes")
    {
#ifdef QMC_COMPLEX
      app_error() << "Orbital optimization via rotation doesn't support complex wavefunction yet.\n";
      abort();
#else
      auto* rot_spo   = new RotatedSPOs(spo);
      xmlNodePtr tcur = cur->xmlChildrenNode;
      while (tcur != NULL)
      {
        std::string cname((const char*)(tcur->name));
        if (cname == "opt_vars")
        {
          rot_spo->params_supplied = true;
          putContent(rot_spo->params, tcur);
        }
        tcur = tcur->next;
      }
      spo = rot_spo;
      spo->objectName = sname;
#endif
    }
    return spo;
  }
  else
  {
    APP_ABORT("SPOSetBuilderFactory::createSPOSet Failed to create a SPOSet. SPOSetBuilder is empty.");
    return 0;
  }
}

void SPOSetBuilderFactory::build_sposet_collection(xmlNodePtr cur)
{
  xmlNodePtr parent = cur;
  std::string type("");
  OhmmsAttributeSet attrib;
  attrib.add(type, "type");
  attrib.put(cur);
  //tolower(type);

  app_log() << "building sposet collection of type " << type << std::endl;

  SPOSetBuilder* bb  = createSPOSetBuilder(cur);
  xmlNodePtr element = parent->children;
  int nsposets       = 0;
  while (element != NULL)
  {
    std::string cname((const char*)(element->name));
    if (cname == "sposet")
    {
      std::string name("");
      OhmmsAttributeSet attrib;
      attrib.add(name, "name");
      attrib.put(element);

      app_log() << "  Building SPOSet \"" << name << "\" with " << type << " SPOSetBuilder" << std::endl;
      SPOSet* spo     = bb->createSPOSet(element);
      spo->objectName = name;
      nsposets++;
    }
    element = element->next;
  }
  if (nsposets == 0)
    APP_ABORT("SPOSetBuilderFactory::build_sposet_collection  no <sposet/> elements found");
}

std::string SPOSetBuilderFactory::basisset_tag = "basisset";

} // namespace qmcplusplus
