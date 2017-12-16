//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCWaveFunctions/Jastrow/PadeBuilder.h"
#include "QMCWaveFunctions/ComboOrbital.h"
#include "QMCWaveFunctions/Jastrow/PadeConstraints.h"

namespace qmcplusplus
{

PadeBuilder::PadeBuilder(ParticleSet& p, TrialWaveFunction& psi,
                         PtclPoolType& psets):
  OrbitalBuilderBase(p,psi), IgnoreSpin(true),
  ptclPool(psets), sourcePtcl(0)
{
}


bool PadeBuilder::put(xmlNodePtr cur)
{
  const xmlChar* spin=xmlGetProp(cur,(const xmlChar*)"spin");
  if(spin != NULL)
  {
    std::string a((const char*)spin);
    if(a == "yes")
      IgnoreSpin=false;
  }
  std::string functionOpt("pade");
  const xmlChar *ftype = xmlGetProp(cur, (const xmlChar *)"function");
  if(ftype != NULL)
    functionOpt = (const char*) ftype;
  const xmlChar* s = xmlGetProp(cur, (const xmlChar *)"source");
  if(s != NULL)
  {
    std::map<std::string,ParticleSet*>::iterator pa_it(ptclPool.find((const char*)s));
    if(pa_it != ptclPool.end())
    {
      sourcePtcl=(*pa_it).second;
    }
  }
  bool success=false;
  OrbitalConstraintsBase* control=0;
  if(functionOpt == "pade")
  {
    app_log() << "  Pade Jastrow Functions = " << functionOpt << std::endl;
    control = new PadeConstraints(IgnoreSpin);
  }
  else
    if(functionOpt == "scaledpade")
    {
      app_log() << "  Scaled Pade Jastrow Functions = " << functionOpt << std::endl;
      control = new ScaledPadeConstraints(IgnoreSpin);
    }
  if(control==0)
    return false;
  if(!control->put(cur))
  {
    delete control;
    return false;
  }
  ComboOrbital* jcombo=new ComboOrbital(control);
  OrbitalBase* j2=control->createTwoBody(targetPtcl);
  jcombo->Psi.push_back(j2);
  if(sourcePtcl)
    // add one-body term using Zeff and e-e B
  {
    OrbitalBase* j1=control->createOneBody(targetPtcl,*sourcePtcl);
    if(j1)
      jcombo->Psi.push_back(j1);
  }
  targetPsi.addOrbital(jcombo);
  return success;
}
}
