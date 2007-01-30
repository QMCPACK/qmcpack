//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#include "QMCWaveFunctions/Jastrow/PadeBuilder.h"
#include "QMCWaveFunctions/ComboOrbital.h"
#include "QMCWaveFunctions/Jastrow/PadeConstraints.h"

namespace qmcplusplus {

  PadeBuilder::PadeBuilder(ParticleSet& p, TrialWaveFunction& psi,
      PtclPoolType& psets):
    OrbitalBuilderBase(p,psi), IgnoreSpin(true),
    ptclPool(psets), sourcePtcl(0){
    }


  bool PadeBuilder::put(xmlNodePtr cur) {

    const xmlChar* spin=xmlGetProp(cur,(const xmlChar*)"spin");
    if(spin != NULL) {
      string a((const char*)spin);
      if(a == "yes") IgnoreSpin=false;
    }

    string functionOpt("pade");
    const xmlChar *ftype = xmlGetProp(cur, (const xmlChar *)"function");
    if(ftype != NULL) functionOpt = (const char*) ftype;

    const xmlChar* s = xmlGetProp(cur, (const xmlChar *)"source");
    if(s != NULL) {
      map<string,ParticleSet*>::iterator pa_it(ptclPool.find((const char*)s));
      if(pa_it != ptclPool.end()) {
        sourcePtcl=(*pa_it).second;
      }
    }

    bool success=false;
    OrbitalConstraintsBase* control=0;
    if(functionOpt == "pade") {
      app_log() << "  Pade Jastrow Functions = " << functionOpt << endl;
      control = new PadeConstraints(IgnoreSpin);
    } else if(functionOpt == "scaledpade") {
      app_log() << "  Scaled Pade Jastrow Functions = " << functionOpt << endl;
      control = new ScaledPadeConstraints(IgnoreSpin);
    }

    if(control==0) return false;

    if(!control->put(cur)) {
      delete control;
      return false;
    }

    ComboOrbital* jcombo=new ComboOrbital(control);
    OrbitalBase* j2=control->createTwoBody(targetPtcl);
    jcombo->Psi.push_back(j2);

    if(sourcePtcl) { // add one-body term using Zeff and e-e B
      OrbitalBase* j1=control->createOneBody(targetPtcl,*sourcePtcl);
      if(j1) jcombo->Psi.push_back(j1);
    }

    targetPsi.addOrbital(jcombo);
    return success;
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
