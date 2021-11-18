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
  for (const auto& sposet_builder : sposet_builders_)
  {
    auto& sposets = sposet_builder->sposets;
    for (auto& sposet : sposets)
      if (sposet->getName() == name)
      {
        spo = sposet.get();
        nfound++;
      }
  }
  if (nfound > 1)
  {
    write_sposet_builders_();
    throw std::runtime_error("getSPOSet: requested sposet " + name + " is not unique!");
  }
  // keep this commented until legacy input styles are moved.
  // In legacy input styles, this look up may fail and need to build SPOSet on the fly.
  //else if (nfound == 0)
  //  throw std::runtime_error("getSPOSet: requested sposet " + name + " is not found!");
  return spo;
}


void SPOSetBuilderFactory::write_sposet_builders_(const std::string& pad) const
{
  std::string pad2 = pad + "  ";
  for (const auto& sposet_builder : sposet_builders_)
  {
    auto& sposets = sposet_builder->sposets;
    app_log() << pad << "sposets for SPOSetBuilder of type " << sposet_builder->getTypeName() << std::endl;
    for (int i = 0; i < sposets.size(); ++i)
      app_log() << pad2 << "sposet " << sposets[i]->getName() << std::endl;
  }
}

SPOSetBuilder& SPOSetBuilderFactory::getLastBuilder()
{
  if (sposet_builders_.empty())
    myComm->barrier_and_abort("SPOSetBuilderFactory::getLastBuilder BUG! No SPOSetBuilder has been created.");
  return *sposet_builders_.back();
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

SPOSetBuilder& SPOSetBuilderFactory::createSPOSetBuilder(xmlNodePtr rootNode)
{
  ReportEngine PRE(ClassName, "createSPOSetBuilder");
  std::string sourceOpt("ion0");
  std::string type("");
  std::string name("");
  OhmmsAttributeSet aAttrib;
  aAttrib.add(sourceOpt, "source");
  aAttrib.add(type, "type");
  aAttrib.add(name, "name");

  if (rootNode != NULL)
    aAttrib.put(rootNode);

  std::string type_in = type;
  tolower(type);

  //when name is missing, type becomes the input
  if (name.empty())
    name = type_in;

  std::unique_ptr<SPOSetBuilder> bb;

  if (type == "composite")
  {
    app_log() << "Composite SPO set with existing SPOSets." << std::endl;
    bb = std::make_unique<CompositeSPOSetBuilder>(myComm, *this);
  }
  else if (type == "jellium" || type == "heg")
  {
    app_log() << "Electron gas SPO set" << std::endl;
    bb = std::make_unique<ElectronGasSPOBuilder>(targetPtcl, myComm, rootNode);
  }
  else if (type == "sho")
  {
    app_log() << "Harmonic Oscillator SPO set" << std::endl;
    bb = std::make_unique<SHOSetBuilder>(targetPtcl, myComm);
  }
#if OHMMS_DIM == 3
  else if (type.find("spline") < type.size())
  {
    if (targetPtcl.is_spinor_)
    {
#ifdef QMC_COMPLEX
      app_log() << "Einspline Spinor Set\n";
      bb = std::make_unique<EinsplineSpinorSetBuilder>(targetPtcl, ptclPool, myComm, rootNode);
#else
      PRE.error("Use of einspline spinors requires QMC_COMPLEX=1.  Rebuild with this option");
#endif
    }
    else
    {
#if defined(HAVE_EINSPLINE)
      PRE << "EinsplineSetBuilder:  using libeinspline for B-spline orbitals.\n";
      bb = std::make_unique<EinsplineSetBuilder>(targetPtcl, ptclPool, myComm, rootNode);
#else
      PRE.error("Einspline is missing for B-spline orbitals", true);
#endif
    }
  }
  else if (type == "molecularorbital" || type == "mo")
  {
    ParticleSet* ions = nullptr;
    //initialize with the source tag
    auto pit(ptclPool.find(sourceOpt));
    if (pit == ptclPool.end())
      PRE.error("Missing basisset/@source.", true);
    else
      ions = (*pit).second;
    if (targetPtcl.is_spinor_)
#ifdef QMC_COMPLEX
      bb = std::make_unique<LCAOSpinorBuilder>(targetPtcl, *ions, myComm, rootNode);
#else
      PRE.error("Use of lcao spinors requires QMC_COMPLEX=1.  Rebuild with this option");
#endif
    else
      bb = std::make_unique<LCAOrbitalBuilder>(targetPtcl, *ions, myComm, rootNode);
  }
#endif //OHMMS_DIM==3
  PRE.flush();

  if (bb == 0)
    myComm->barrier_and_abort("SPOSetBuilderFactory::createSPOSetBuilder SPOSetBuilder creation failed.");

  app_log() << "  Created SPOSet builder named '" << name << "' of type " << type << std::endl;
  sposet_builders_.push_back(std::move(bb));

  return *sposet_builders_.back();
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
  auto& bb = createSPOSetBuilder(cur);

  // going through a list of sposet entries
  int nsposets = 0;
  processChildren(cur, [&](const std::string& cname, const xmlNodePtr element) {
    if (cname == "sposet")
    {
      SPOSet* spo = bb.createSPOSet(element);
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
        SPOSetScanner ascanner(bb.sposets, targetPtcl, ptclPool);
        ascanner.put(element);
      }
  });
}

std::string SPOSetBuilderFactory::basisset_tag = "basisset";

} // namespace qmcplusplus
