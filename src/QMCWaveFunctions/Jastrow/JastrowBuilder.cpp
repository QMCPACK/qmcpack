//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "JastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/kSpaceJastrowBuilder.h"
#if OHMMS_DIM == 3
#include "QMCWaveFunctions/Jastrow/eeI_JastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/CountingJastrowBuilder.h"
#endif
#include "QMCWaveFunctions/Jastrow/RadialJastrowBuilder.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{
JastrowBuilder::JastrowBuilder(Communicate* comm, ParticleSet& p, const PtclPoolType& psets)
    : WaveFunctionComponentBuilder(comm, p), ptclPool(psets)
{
  resetOptions();
  ClassName = "JastrowBuilder";
}

void JastrowBuilder::resetOptions()
{
  JastrowType  = 0;
  nameOpt      = "0";
  typeOpt      = "Two";
  funcOpt      = "any";
  spinOpt      = "yes";
  transformOpt = "no";
  sourceOpt    = targetPtcl.getName();
}

std::unique_ptr<WaveFunctionComponent> JastrowBuilder::buildComponent(xmlNodePtr cur)
{
  myNode = cur;
  resetOptions();
  OhmmsAttributeSet oAttrib;
  oAttrib.add(typeOpt, "type");
  oAttrib.add(nameOpt, "name");
  oAttrib.add(funcOpt, "function");
  oAttrib.add(transformOpt, "transform");
  oAttrib.add(sourceOpt, "source");
  oAttrib.add(spinOpt, "spin");
  oAttrib.put(cur);
  if (nameOpt[0] == '0')
  {
    APP_ABORT("  JastrowBuilder::put does not have name!\n");
    return nullptr;
  }
  app_summary() << std::endl;
  app_summary() << "   Jastrow" << std::endl;
  app_summary() << "   -------" << std::endl;
  app_summary() << "    Name: " << nameOpt << "   Type: " << typeOpt << "   Function: " << funcOpt << std::endl;
  app_summary() << std::endl;

  if (typeOpt.find("One") < typeOpt.size())
    return buildOneBody(cur);
  else if (typeOpt.find("Two") < typeOpt.size())
    return buildTwoBody(cur);
  else if (typeOpt.find("eeI") < typeOpt.size())
    return build_eeI(cur);
  else if (typeOpt.find("kSpace") < typeOpt.size())
    return buildkSpace(cur);
  else if (typeOpt.find("Counting") < typeOpt.size())
    return buildCounting(cur);
  else
  {
    APP_ABORT("  JastrowBuilder::buildComponent unknown type!\n");
    return nullptr;
  }
}

std::unique_ptr<WaveFunctionComponent> JastrowBuilder::buildCounting(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "addCounting(xmlNodePtr)");
  std::unique_ptr<CountingJastrowBuilder> cjb;
  auto pa_it(ptclPool.find(sourceOpt));
  if (pa_it != ptclPool.end() && sourceOpt != targetPtcl.getName()) // source is not target
  {
    ParticleSet* sourcePtcl = (*pa_it).second;
    cjb                     = std::make_unique<CountingJastrowBuilder>(myComm, targetPtcl, *sourcePtcl);
  }
  else
    cjb = std::make_unique<CountingJastrowBuilder>(myComm, targetPtcl);
  return cjb->buildComponent(cur);
}

std::unique_ptr<WaveFunctionComponent> JastrowBuilder::buildkSpace(xmlNodePtr cur)
{
  app_log() << "  JastrowBuilder::buildkSpace(xmlNodePtr)" << std::endl;
  auto pa_it(ptclPool.find(sourceOpt));
  if (pa_it == ptclPool.end())
  {
    app_warning() << "  JastrowBuilder::buildkSpace failed. " << sourceOpt << " does not exist" << std::endl;
    return nullptr;
  }
  ParticleSet* sourcePtcl = (*pa_it).second;
  app_log() << "\n  Using kSpaceJastrowBuilder for reciprocal-space Jastrows" << std::endl;
  kSpaceJastrowBuilder sBuilder(myComm, targetPtcl, *sourcePtcl);
  return sBuilder.buildComponent(cur);
}

std::unique_ptr<WaveFunctionComponent> JastrowBuilder::buildOneBody(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "addOneBody(xmlNodePtr)");
  if (sourceOpt == targetPtcl.getName())
  {
    PRE.error("One-Body Jastrow Function needs a source different from " + targetPtcl.getName() +
              "\nExit JastrowBuilder::buildOneBody.\n");
    return nullptr;
  }
  auto pa_it(ptclPool.find(sourceOpt));
  if (pa_it == ptclPool.end())
  {
    PRE.error("JastrowBuilder::buildOneBody failed. " + sourceOpt + " does not exist.");
    return nullptr;
  }
  ParticleSet* sourcePtcl = (*pa_it).second;
  //use lowercase, to be handled by parser later
  RadialJastrowBuilder rb(myComm, targetPtcl, *sourcePtcl);
  return rb.buildComponent(cur);
}

std::unique_ptr<WaveFunctionComponent> JastrowBuilder::build_eeI(xmlNodePtr cur)
{
#if OHMMS_DIM == 3
  ReportEngine PRE(ClassName, "add_eeI(xmlNodePtr)");
  auto pit(ptclPool.find(sourceOpt));
  if (pit == ptclPool.end())
  {
    app_error() << "     JastrowBuilder::build_eeI requires a source attribute. " << sourceOpt << " is invalid "
                << std::endl;
    APP_ABORT("  JastrowBuilder::build_eeI");
    return nullptr;
  }
  ParticleSet& sourcePtcl = *((*pit).second);
  eeI_JastrowBuilder jb(myComm, targetPtcl, sourcePtcl);
  return jb.buildComponent(cur);
#else
  APP_ABORT("  eeI is not valid for OHMMS_DIM != 3 ");
  return nullptr;
#endif
}


std::unique_ptr<WaveFunctionComponent> JastrowBuilder::buildTwoBody(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "addTwoBody(xmlNodePtr)");
  RadialJastrowBuilder rb(myComm, targetPtcl);
  return rb.buildComponent(cur);
}

} // namespace qmcplusplus
