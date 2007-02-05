//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowBuilder.h"
#include "QMCWaveFunctions/ComboOrbital.h"
#include "QMCWaveFunctions/Jastrow/PadeConstraints.h"
#include "QMCWaveFunctions/Jastrow/RPAConstraints.h"
#include "QMCWaveFunctions/Jastrow/JAABuilder.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

  TwoBodyJastrowBuilder::TwoBodyJastrowBuilder(ParticleSet& p, TrialWaveFunction& psi,
      PtclPoolType& psets):
    OrbitalBuilderBase(p,psi), ptclPool(psets)
    { }


  bool TwoBodyJastrowBuilder::put(xmlNodePtr cur) {

    myNode=cur; 

    string functionOpt("pade");
    string transformOpt("no");
    string sourceOpt(targetPtcl.getName());
    string spinOpt("yes");
    OhmmsAttributeSet oAttrib;
    oAttrib.add(functionOpt,"function");
    oAttrib.add(transformOpt,"transform");
    oAttrib.add(sourceOpt,"source");
    oAttrib.add(spinOpt,"spin");
    oAttrib.put(cur);

    bool IgnoreSpin = (spinOpt == "no");


    bool success=false;
    OrbitalConstraintsBase* control=0;

    //@todo automatically set it to yes with PBC
    bool useSpline= (transformOpt == "yes");

    app_log() << "  TwoBodyJastrowBuilder for " << functionOpt << endl;

    if(functionOpt == "pade") 
    {
      app_log() << "    Using analytic Pade Jastrow Functor " <<endl;
      control = new PadeConstraints(targetPtcl,targetPsi,IgnoreSpin);
    } 
    else if(functionOpt == "scaledpade") 
    {
      app_log() << "    Using analytic Scaled Pade Jastrow Functor " <<endl;
      control = new ScaledPadeConstraints(targetPtcl,targetPsi,IgnoreSpin);
    } 
    else if(functionOpt == "rpa") 
    {
      if(useSpline) 
        control = new RPAPBCConstraints(targetPtcl,targetPsi,IgnoreSpin);
      else 
        control = new RPAConstraints(targetPtcl,targetPsi,IgnoreSpin);
    } 
    else //known analytic function
    {
      OrbitalBuilderBase* jbuilder=0;
      jbuilder = new JAABuilder(targetPtcl,targetPsi);
      Children.push_back(jbuilder);
      return jbuilder->put(cur);
    }

    success=control->put(cur);
    if(!success) {
      delete control;
      return false;
    }

    ComboOrbital* jcombo=new ComboOrbital(control);
    control->addTwoBodyPart(jcombo);

    if(sourceOpt != targetPtcl.getName())
    {
      app_log() << "    Adding one-body Jastrow function dependent upon two-body " << functionOpt << endl;
      map<string,ParticleSet*>::iterator pa_it(ptclPool.find(sourceOpt));
      if(pa_it == ptclPool.end()) 
      {
        return false;
      }
      ParticleSet* sourcePtcl= sourcePtcl=(*pa_it).second;
      OrbitalBase* j1=control->createOneBody(*sourcePtcl);
      if(j1) jcombo->Psi.push_back(j1);
    }

    control->addOptimizables(targetPsi.VarList);
    targetPsi.addOrbital(jcombo);
    return success;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
