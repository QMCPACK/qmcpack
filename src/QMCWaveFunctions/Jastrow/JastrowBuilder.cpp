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
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCWaveFunctions/Jastrow/JastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/PadeJastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/RPAJastrow.h"
#include "QMCWaveFunctions/Jastrow/BsplineJastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/kSpaceJastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/singleRPAJastrowBuilder.h"
#if OHMMS_DIM ==3
#include "QMCWaveFunctions/Jastrow/eeI_JastrowBuilder.h"
#endif
#include "QMCWaveFunctions/Jastrow/JABBuilder.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

JastrowBuilder::JastrowBuilder(ParticleSet& p, TrialWaveFunction& psi, PtclPoolType& psets):
  OrbitalBuilderBase(p,psi), ptclPool(psets)
{
  resetOptions();
  ClassName="JastrowBuilder";
}

void JastrowBuilder::resetOptions()
{
  JastrowType = 0;
  nameOpt="0";
  typeOpt="Two";
  funcOpt="any";
  spinOpt="yes";
  transformOpt="no";
  sourceOpt=targetPtcl.getName();
}

bool JastrowBuilder::put(xmlNodePtr cur)
{
  myNode=cur;
  resetOptions();
  OhmmsAttributeSet oAttrib;
  oAttrib.add(typeOpt,"type");
  oAttrib.add(nameOpt,"name");
  oAttrib.add(funcOpt,"function");
  oAttrib.add(transformOpt,"transform");
  oAttrib.add(sourceOpt,"source");
  oAttrib.add(spinOpt,"spin");
  oAttrib.put(cur);
  if(nameOpt[0] == '0')
  {
    app_warning() << "  JastrowBuilder::put does not have name "<< std::endl;
    return false;
  }
  if(typeOpt.find("One") < typeOpt.size())
    return addOneBody(cur);
  if(typeOpt.find("Two") < typeOpt.size())
    return addTwoBody(cur);
  if(typeOpt.find("eeI") < typeOpt.size())
    return add_eeI(cur);
  if(typeOpt.find("kSpace") < typeOpt.size())
    return addkSpace(cur);
  return false;
}

bool JastrowBuilder::addkSpace(xmlNodePtr cur)
{
  app_log() << "  JastrowBuilder::addkSpace(xmlNodePtr)" << std::endl;
  std::map<std::string,ParticleSet*>::iterator pa_it(ptclPool.find(sourceOpt));
  if(pa_it == ptclPool.end())
  {
    app_warning() << "  JastrowBuilder::addkSpace failed. "
                  << sourceOpt << " does not exist" << std::endl;
    return false;
  }
  ParticleSet* sourcePtcl= (*pa_it).second;
  app_log() << "\n  Using kSpaceJastrowBuilder for reciprocal-space Jastrows" << std::endl;
  OrbitalBuilderBase* sBuilder = new kSpaceJastrowBuilder (targetPtcl, targetPsi, *sourcePtcl);
  Children.push_back(sBuilder);
  return sBuilder->put(cur);
}

bool JastrowBuilder::addOneBody(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName,"addOneBody(xmlNodePtr)");
  if(sourceOpt == targetPtcl.getName())
  {
    PRE.warning("One-Body Jastrow Function needs a source different from "+targetPtcl.getName()
                +"\nExit JastrowBuilder::addOneBody.\n");
    return false;
  }
  std::map<std::string,ParticleSet*>::iterator pa_it(ptclPool.find(sourceOpt));
  if(pa_it == ptclPool.end())
  {
    PRE.warning("JastrowBuilder::addOneBody failed. "+sourceOpt+" does not exist.");
    return false;
  }
  bool success=false;
  ParticleSet* sourcePtcl= (*pa_it).second;
  //use lowercase, to be handled by parser later
  tolower(funcOpt);
  if (funcOpt == "bspline" )
  {
    app_log() << "\n  Using BsplineBuilder for one-body jastrow with B-spline functions" << std::endl;
    BsplineJastrowBuilder jb(targetPtcl,targetPsi,*sourcePtcl);
    success=jb.put(cur);
  }
  else
    if (funcOpt == "rpa" )
    {
#if OHMMS_DIM ==3
      app_log() << "\n  Using RPA for one-body jastrow" << std::endl;
      singleRPAJastrowBuilder jb(targetPtcl, targetPsi, *sourcePtcl);
      success= jb.put(cur);
#else
      APP_ABORT("RPA for one-body jastrow is only available for 3D");
#endif
    }
    else
    {
      app_log() << "\n  Using JABBuilder for one-body jastrow with analytic functions" << std::endl;
      JABBuilder jb(targetPtcl,targetPsi,ptclPool);
      success=jb.put(cur);
    }
  return success;
}

bool JastrowBuilder::add_eeI (xmlNodePtr cur)
{
#if OHMMS_DIM ==3
  ReportEngine PRE(ClassName,"add_eeI(xmlNodePtr)");
  PtclPoolType::iterator pit(ptclPool.find(sourceOpt));
  if(pit == ptclPool.end())
  {
    app_error() << "     JastrowBuilder::add_eeI requires a source attribute. "
                << sourceOpt << " is invalid " << std::endl;
    APP_ABORT("  JastrowBuilder::add_eeI");
    return false;
  }
  ParticleSet& sourcePtcl= *((*pit).second);
  eeI_JastrowBuilder jb(targetPtcl, targetPsi, sourcePtcl);
  return jb.put (cur);
#else
  APP_ABORT("  eeI is not valid for OHMMS_DIM != 3 ");
  return true;
#endif
}


bool JastrowBuilder::addTwoBody(xmlNodePtr cur)
{
  ReportEngine PRE(ClassName,"addTwoBody(xmlNodePtr)");
  bool success=false;
  bool useSpline = (targetPtcl.Lattice.BoxBConds[0] && transformOpt == "yes");
  bool ignoreSpin = (spinOpt == "no");
  //convert to lowercase
  tolower(funcOpt);
  if(funcOpt == "pade")
  {
    if (targetPtcl.Lattice.SuperCellEnum != SUPERCELL_OPEN)
    {
      PRE.warning("Pade Jastrow is requested for a periodic system. Please choose other functors.");
      return false;
    }
    PadeJastrowBuilder pbuilder(targetPtcl,targetPsi,ptclPool);
    return pbuilder.put(cur);
  }
  else if((funcOpt == "yukawa") || (funcOpt == "rpa"))
  {
    if(targetPtcl.Lattice.SuperCellEnum == SUPERCELL_OPEN)
    {
      PRE.warning("RPA is requested for an open system. Please choose other functors.");
      return false;
    }
    RPAJastrow* rpajastrow = new RPAJastrow(targetPtcl,targetPsi.is_manager());
    rpajastrow->put(cur);
    targetPsi.addOrbital(rpajastrow,nameOpt);
    return true;
  }
  else if (funcOpt == "bspline" )
  {
    BsplineJastrowBuilder bbuilder(targetPtcl,targetPsi);
    return bbuilder.put(cur);
  }
  else
  {
    app_error() << "Unknown two body function: " << funcOpt << ".\n";
  }
  
  return success;
}



}
