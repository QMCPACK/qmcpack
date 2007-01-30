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
#include "QMCWaveFunctions/TwoBodyJastrowBuilder.h"
#include "QMCWaveFunctions/ComboOrbital.h"
#include "QMCWaveFunctions/PadeConstraints.h"
#include "QMCWaveFunctions/RPAConstraints.h"
#include "QMCWaveFunctions/WMConstraints.h"
#include "QMCWaveFunctions/JAABuilder.h"
#include "QMCWaveFunctions/NJAABuilder.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

  TwoBodyJastrowBuilder::TwoBodyJastrowBuilder(ParticleSet& p, TrialWaveFunction& psi,
      PtclPoolType& psets):
    OrbitalBuilderBase(p,psi), IgnoreSpin(true),
    ptclPool(psets), sourcePtcl(0){
    }


  bool TwoBodyJastrowBuilder::put(xmlNodePtr cur) {

    string functionOpt("pade");
    string transformOpt("no");
    string sourceOpt("9NONE");
    string spinOpt("yes");
    OhmmsAttributeSet oAttrib;
    oAttrib.add(functionOpt,"function");
    oAttrib.add(transformOpt,"transform");
    oAttrib.add(sourceOpt,"source");
    oAttrib.add(spinOpt,"spin");
    oAttrib.put(cur);

    IgnoreSpin = (spinOpt == "no");

    if(sourceOpt[0] != '9') {
      map<string,ParticleSet*>::iterator pa_it(ptclPool.find(sourceOpt));
      if(pa_it != ptclPool.end()) {
        sourcePtcl=(*pa_it).second;
      }
    }

    bool success=false;
    OrbitalConstraintsBase* control=0;

    //@todo automatically set it to yes with PBC
    bool useSpline= (transformOpt == "yes");

    app_log() << "  TwoBodyJastrowBuilder for " << functionOpt << endl;

    if(functionOpt == "pade") {
      //control = new PadeConstraints(IgnoreSpin);
      //Transform is ignored. Cutoff function is not too good
      if(useSpline) {
        control = new PadeOnGridConstraints(IgnoreSpin);
      } else {
        control = new PadeConstraints(IgnoreSpin);
      }
    } else if(functionOpt == "scaledpade") {
      control = new ScaledPadeConstraints(IgnoreSpin);
    } else if(functionOpt == "rpa") {
      if(useSpline) {
        control = new RPAPBCConstraints(IgnoreSpin);
      } else {
        control = new RPAConstraints(IgnoreSpin);
      }
    } else if(functionOpt == "WM") {
      control = new WMConstraints(IgnoreSpin);
    }

     if(control==0) { //try generic JAABuilder and NJAABuilder
      OrbitalBuilderBase* jbuilder=0;
      if(useSpline) {
        jbuilder = new NJAABuilder(targetPtcl,targetPsi);
      } else {
        jbuilder = new JAABuilder(targetPtcl,targetPsi);
      }
      return jbuilder->put(cur);
    }

    success=control->put(cur);
    if(!control->put(cur)) {
      delete control;
      return false;
    }

    ComboOrbital* jcombo=new ComboOrbital(control);

    control->addTwoBodyPart(targetPtcl, jcombo);

    if(sourcePtcl) { // add one-body term using Zeff and e-e B
      OrbitalBase* j1=control->createOneBody(targetPtcl,*sourcePtcl);
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
