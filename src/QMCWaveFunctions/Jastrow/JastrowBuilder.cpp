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


#include "QMCWaveFunctions/Jastrow/JastrowBuilder.h"
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
JastrowBuilder::JastrowBuilder(ParticleSet& p, TrialWaveFunction& psi, PtclPoolType& psets)
    : WaveFunctionComponentBuilder(p, psi), ptclPool(psets)
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

bool JastrowBuilder::put(xmlNodePtr cur)
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
    app_warning() << "  JastrowBuilder::put does not have name " << std::endl;
    return false;
  }
  app_summary() << std::endl;
  app_summary() << "   Jastrow" << std::endl;
  app_summary() << "   -------" << std::endl;
  app_summary() << "    Name: " << nameOpt << "   Type: " << typeOpt << "   Function: " << funcOpt << std::endl;
  app_summary() << std::endl;

  if (typeOpt.find("One") < typeOpt.size())
    return addOneBody(cur);
  if (typeOpt.find("Two") < typeOpt.size())
    return addTwoBody(cur);
  if (typeOpt.find("eeI") < typeOpt.size())
    return add_eeI(cur);
  if (typeOpt.find("kSpace") < typeOpt.size())
    return addkSpace(cur);
  if (typeOpt.find("Counting") < typeOpt.size())
    return addCounting(cur);
  return false;
}

bool JastrowBuilder::addCounting(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "addCounting(xmlNodePtr)");
  CountingJastrowBuilder* cjb;
  std::map<std::string, ParticleSet*>::iterator pa_it(ptclPool.find(sourceOpt));
  if (pa_it != ptclPool.end() && sourceOpt != targetPtcl.getName()) // source is not target
  {
    ParticleSet* sourcePtcl = (*pa_it).second;
    cjb = new CountingJastrowBuilder(targetPtcl, targetPsi, *sourcePtcl);
  }
  else
    cjb = new CountingJastrowBuilder(targetPtcl, targetPsi);
  return cjb->put(cur);
}

bool JastrowBuilder::addkSpace(xmlNodePtr cur)
{
  app_log() << "  JastrowBuilder::addkSpace(xmlNodePtr)" << std::endl;
  std::map<std::string, ParticleSet*>::iterator pa_it(ptclPool.find(sourceOpt));
  if (pa_it == ptclPool.end())
  {
    app_warning() << "  JastrowBuilder::addkSpace failed. " << sourceOpt << " does not exist" << std::endl;
    return false;
  }
  ParticleSet* sourcePtcl = (*pa_it).second;
  app_log() << "\n  Using kSpaceJastrowBuilder for reciprocal-space Jastrows" << std::endl;
  WaveFunctionComponentBuilder* sBuilder = new kSpaceJastrowBuilder(targetPtcl, targetPsi, *sourcePtcl);
  Children.push_back(sBuilder);
  return sBuilder->put(cur);
}

bool JastrowBuilder::addOneBody(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "addOneBody(xmlNodePtr)");
  if (sourceOpt == targetPtcl.getName())
  {
    PRE.warning("One-Body Jastrow Function needs a source different from " + targetPtcl.getName() +
                "\nExit JastrowBuilder::addOneBody.\n");
    return false;
  }
  std::map<std::string, ParticleSet*>::iterator pa_it(ptclPool.find(sourceOpt));
  if (pa_it == ptclPool.end())
  {
    PRE.warning("JastrowBuilder::addOneBody failed. " + sourceOpt + " does not exist.");
    return false;
  }
  bool success            = false;
  ParticleSet* sourcePtcl = (*pa_it).second;
  //use lowercase, to be handled by parser later
  RadialJastrowBuilder rb(targetPtcl, targetPsi, *sourcePtcl);
  success = rb.put(cur);
  return success;
}

bool JastrowBuilder::add_eeI(xmlNodePtr cur)
{
#if OHMMS_DIM == 3
  ReportEngine PRE(ClassName, "add_eeI(xmlNodePtr)");
  PtclPoolType::iterator pit(ptclPool.find(sourceOpt));
  if (pit == ptclPool.end())
  {
    app_error() << "     JastrowBuilder::add_eeI requires a source attribute. " << sourceOpt << " is invalid "
                << std::endl;
    APP_ABORT("  JastrowBuilder::add_eeI");
    return false;
  }
  ParticleSet& sourcePtcl = *((*pit).second);
  eeI_JastrowBuilder jb(targetPtcl, targetPsi, sourcePtcl);
  return jb.put(cur);
#else
  APP_ABORT("  eeI is not valid for OHMMS_DIM != 3 ");
  return true;
#endif
}


bool JastrowBuilder::addTwoBody(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName, "addTwoBody(xmlNodePtr)");
  RadialJastrowBuilder rb(targetPtcl, targetPsi);
  return rb.put(cur);
}

} // namespace qmcplusplus
