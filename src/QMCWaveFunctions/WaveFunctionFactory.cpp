//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/**@file WaveFunctionFactory.cpp
 *@brief Definition of a WaveFunctionFactory
 */
#include "WaveFunctionFactory.h"
#include "QMCWaveFunctions/Jastrow/JastrowBuilder.h"
#include "QMCWaveFunctions/Fermion/SlaterDetBuilder.h"
#include "QMCWaveFunctions/LatticeGaussianProductBuilder.h"
#include "QMCWaveFunctions/ExampleHeBuilder.h"

#if defined(QMC_COMPLEX)
#include "QMCWaveFunctions/ElectronGas/ElectronGasComplexOrbitalBuilder.h"
#else
#include "QMCWaveFunctions/ElectronGas/ElectronGasOrbitalBuilder.h"
#endif

#include "QMCWaveFunctions/PlaneWave/PWOrbitalBuilder.h"
#if OHMMS_DIM == 3 && !defined(QMC_COMPLEX)
#include "QMCWaveFunctions/AGPDeterminantBuilder.h"
#endif

#include "Utilities/ProgressReportEngine.h"
#include "Utilities/IteratorUtility.h"
#include "OhmmsData/AttributeSet.h"
namespace qmcplusplus
{
WaveFunctionFactory::WaveFunctionFactory(const std::string& psiName,
                                         ParticleSet& qp,
                                         PtclPoolType& pset,
                                         Communicate* c,
                                         bool tasking)
    : MPIObjectBase(c),
      targetPsi(std::make_unique<TrialWaveFunction>(psiName, tasking)),
      targetPtcl(qp),
      ptclPool(pset),
      myNode(NULL),
      sposet_builder_factory_(c, qp, pset)
{
  ClassName = "WaveFunctionFactory";
  myName    = psiName;
  targetPsi->setMassTerm(targetPtcl);
}

bool WaveFunctionFactory::build(xmlNodePtr cur, bool buildtree)
{
  app_summary() << std::endl;
  app_summary() << " Many-body wavefunction" << std::endl;
  app_summary() << " -------------------" << std::endl;
  app_summary() << "  Name: " << myName << "   Tasking: " << (targetPsi->use_tasking() ? "yes" : "no") << std::endl;
  app_summary() << std::endl;

  ReportEngine PRE(ClassName, "build");
  if (cur == NULL)
    return false;
  bool attach2Node = false;
  if (buildtree)
  {
    if (myNode == NULL)
    {
      myNode = xmlCopyNode(cur, 1);
    }
    else
    {
      attach2Node = true;
    }
  }
  cur          = cur->children;
  bool success = true;
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "sposet_builder" || cname == "sposet_collection")
      sposet_builder_factory_.buildSPOSetCollection(cur);
    else if (cname == WaveFunctionComponentBuilder::detset_tag)
    {
      success = addFermionTerm(cur);
      bool foundtwist(false);
      xmlNodePtr kcur = cur->children;
      while (kcur != NULL)
      {
        std::string kname((const char*)(kcur->name));
        if (kname == "h5tag")
        {
          std::string hdfName;
          OhmmsAttributeSet attribs;
          attribs.add(hdfName, "name");
          if (hdfName == "twistAngle")
          {
            std::vector<ParticleSet::RealType> tsts(3, 0);
            putContent(tsts, kcur);
            targetPsi->setTwist(tsts);
            foundtwist = true;
          }
        }
        kcur = kcur->next;
      }
      if (!foundtwist)
      {
        //default twist is [0 0 0]
        std::vector<ParticleSet::RealType> tsts(3, 0);
        targetPsi->setTwist(tsts);
      }
    }
    else if (cname == WaveFunctionComponentBuilder::jastrow_tag)
    {
      WaveFunctionComponentBuilder* jbuilder = new JastrowBuilder(myComm, targetPtcl, ptclPool);
      targetPsi->addComponent(jbuilder->buildComponent(cur));
      success = true;
      addNode(jbuilder, cur);
    }
    else if (cname == "fdlrwfn")
    {
      APP_ABORT("FDLR wave functions are not currently supported.");
    }
    else if (cname == WaveFunctionComponentBuilder::ionorb_tag)
    {
      LatticeGaussianProductBuilder* builder = new LatticeGaussianProductBuilder(myComm, targetPtcl, ptclPool);
      targetPsi->addComponent(builder->buildComponent(cur));
      success = true;
      addNode(builder, cur);
    }
    else if ((cname == "Molecular") || (cname == "molecular"))
    {
      APP_ABORT("  Removed Helium Molecular terms from qmcpack ");
      success = false;
    }
    else if (cname == "example_he")
    {
      WaveFunctionComponentBuilder* exampleHe_builder = new ExampleHeBuilder(myComm, targetPtcl, ptclPool);
      targetPsi->addComponent(exampleHe_builder->buildComponent(cur));
      success = true;
      addNode(exampleHe_builder, cur);
    }
#if !defined(QMC_COMPLEX) && OHMMS_DIM == 3
    else if (cname == "agp")
    {
      AGPDeterminantBuilder* agpbuilder = new AGPDeterminantBuilder(myComm, targetPtcl, ptclPool);
      targetPsi->addComponent(agpbuilder->buildComponent(cur));
      success = true;
      addNode(agpbuilder, cur);
    }
#endif
    if (attach2Node)
      xmlAddChild(myNode, xmlCopyNode(cur, 1));
    cur = cur->next;
  }
  //{
  //  ReportEngine PREA("TrialWaveFunction","print");
  //  targetPsi->VarList.print(app_log());
  //}
  // synch all parameters. You don't want some being different if same name.
  opt_variables_type dummy;
  targetPsi->checkInVariables(dummy);
  dummy.resetIndex();
  targetPsi->checkOutVariables(dummy);
  targetPsi->resetParameters(dummy);
  return success;
}


bool WaveFunctionFactory::addFermionTerm(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "addFermionTerm");
  std::string orbtype("MolecularOrbital");
  std::string nuclei("i");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(orbtype, "type");
  oAttrib.add(nuclei, "source");
  oAttrib.put(cur);
  WaveFunctionComponentBuilder* detbuilder = 0;
  if (orbtype == "electron-gas")
  {
#if defined(QMC_COMPLEX)
    detbuilder = new ElectronGasComplexOrbitalBuilder(myComm, targetPtcl);
#else
    detbuilder = new ElectronGasOrbitalBuilder(myComm, targetPtcl);
#endif
  }
  else if (orbtype == "PWBasis" || orbtype == "PW" || orbtype == "pw")
  {
    detbuilder = new PWOrbitalBuilder(myComm, targetPtcl, ptclPool);
  }
  else
    detbuilder = new SlaterDetBuilder(myComm, sposet_builder_factory_, targetPtcl, *targetPsi, ptclPool);
  targetPsi->addComponent(detbuilder->buildComponent(cur));
  addNode(detbuilder, cur);
  return true;
}


bool WaveFunctionFactory::addNode(WaveFunctionComponentBuilder* b, xmlNodePtr cur)
{
  psiBuilder.emplace_back(b);
  ///if(myNode != NULL) {
  ///  std::cout << ">>>> Adding " << (const char*)cur->name << std::endl;
  ///  xmlAddChild(myNode,xmlCopyNode(cur,1));
  ///}
  return true;
}

bool WaveFunctionFactory::put(xmlNodePtr cur) { return build(cur, true); }

} // namespace qmcplusplus
