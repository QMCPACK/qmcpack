//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCWaveFunctions/Jastrow/JastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/JastrowBasisBuilder.h"
#include "QMCWaveFunctions/Jastrow/AnyConstraints.h"
#include "QMCWaveFunctions/Jastrow/JABBuilder.h"
#include "QMCWaveFunctions/Jastrow/JAAPBCBuilder.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowBuilder.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

  JastrowBuilder::JastrowBuilder(ParticleSet& p, TrialWaveFunction& psi, PtclPoolType& psets):
    OrbitalBuilderBase(p,psi), ptclPool(psets)
  {
    resetOptions();
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

  bool JastrowBuilder::put(xmlNodePtr cur) {

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
      app_warning() << "  WaveFunctionFactory::addJastrowTerm missing type. Ignore " << nameOpt << endl;
      return false;
    }

    if(typeOpt.find("One") < typeOpt.size()) 
      return addOneBody(cur);

    if(typeOpt.find("Two") < typeOpt.size()) 
      return addTwoBody(cur);

    if(typeOpt.find("Three") < typeOpt.size()) 
      return addThreeBody(cur);

    return false;
  }

  bool JastrowBuilder::addOneBody(xmlNodePtr cur) 
  {
    app_log() << "  JastrowBuilder::addOneBody "<< endl;

    if(sourceOpt == targetPtcl.getName()) 
    {
      app_warning() << "  One-Body Jastrow Function needs a source different from " 
        << targetPtcl.getName() << endl;
      return false;
    }

    map<string,ParticleSet*>::iterator pa_it(ptclPool.find(sourceOpt));
    if(pa_it == ptclPool.end()) {
      app_warning() << "  JastrowBuilder::addOneBody failed. " << sourceOpt << " does not exist"
        << endl;
      return false;
    }

    ParticleSet* sourcePtcl= (*pa_it).second;

    if(funcOpt == "any")
    {
      AnyConstraints* control=new AnyConstraints(targetPtcl,targetPsi);
      control->put(cur);
      OrbitalBase* j=control->createOneBody(*sourcePtcl);
      if(j)
      {
        control->addOptimizables(targetPsi.VarList);
        targetPsi.addOrbital(j);
        Children.push_back(control);
        return true;
      }
      else
      {
        delete control;
        return false;
      }
    }
    else
    {
      app_log() << "\n  Using JABBuilder for one-body jatrow with analytic functions" << endl;
      OrbitalBuilderBase* jb = new JABBuilder(targetPtcl,targetPsi,ptclPool);
      Children.push_back(jb);
      return jb->put(cur);
    }
  }

  bool JastrowBuilder::addTwoBody(xmlNodePtr cur) 
  {
    app_log() << "  JastrowBuilder::addTwoBody "<< endl;

    if(funcOpt == "any")
    {
      AnyConstraints* control=new AnyConstraints(targetPtcl,targetPsi);
      control->put(cur);
      OrbitalBase* j=control->createTwoBody();
      if(j)
      {
        control->addOptimizables(targetPsi.VarList);
        targetPsi.addOrbital(j);
        Children.push_back(control);
        return true;
      }
      else
      {
        delete control;
        return false;
      }
    }
    else
    {
      OrbitalBuilderBase* jb = new TwoBodyJastrowBuilder(targetPtcl,targetPsi,ptclPool);
      Children.push_back(jb);
      return jb->put(cur);
    }
  }

  bool JastrowBuilder::addThreeBody(xmlNodePtr cur) 
  {
    app_log() << "  JastrowBuilder::addThreeBody "<< endl;
//    if(jasttype == "Three-Body-Geminal") {
//      app_log() << "\n  creating Three-Body-Germinal Jastrow function " << endl;
//      string source_name("i");
//      const xmlChar* iptr = xmlGetProp(cur, (const xmlChar *)"source");
//      if(iptr != NULL) source_name=(const char*)iptr;
//      PtclPoolType::iterator pit(ptclPool.find(source_name));
//      if(pit != ptclPool.end()) {
//        jbuilder = new ThreeBodyGeminalBuilder(*targetPtcl,*targetPsi,*((*pit).second));
//      }
//    } else if (jasttype == "Three-Body-Pade") {
//      app_log() << "\n  creating Three-Body-Pade Jastrow function " << endl;
//      string source_name("i");
//      const xmlChar* iptr = xmlGetProp(cur, (const xmlChar *)"source");
//      //if(iptr != NULL) source_name=(const char*)iptr;
//      PtclPoolType::iterator pit(ptclPool.find(source_name));
//      if(pit != ptclPool.end()) {
//        jbuilder = new ThreeBodyPadeBuilder(*targetPtcl,*targetPsi,*((*pit).second));
//      }
//    }
    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1669 $   $Date: 2007-01-30 13:21:54 -0600 (Tue, 30 Jan 2007) $
 * $Id: JastrowBuilder.cpp 1669 2007-01-30 19:21:54Z jnkim $ 
 ***************************************************************************/
