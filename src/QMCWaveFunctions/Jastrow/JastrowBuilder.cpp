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
#include "QMCWaveFunctions/Jastrow/PadeConstraints.h"
#include "QMCWaveFunctions/Jastrow/RPAConstraints.h"
#include "QMCWaveFunctions/Jastrow/JAABuilder.h"
#include "QMCWaveFunctions/Jastrow/JABBuilder.h"
#include "QMCWaveFunctions/Jastrow/ThreeBodyGeminal.h"
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
      app_warning() << "  JastrowBuilder::put does not have name "<< endl;
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
      app_warning() << "  One-Body Jastrow Function needs a source different from " << targetPtcl.getName() << endl;
      app_warning() << "  Exit JastrowBuilder::addOneBody." << endl;
      return false;
    }

    map<string,ParticleSet*>::iterator pa_it(ptclPool.find(sourceOpt));
    if(pa_it == ptclPool.end()) {
      app_warning() << "  JastrowBuilder::addOneBody failed. " << sourceOpt << " does not exist" << endl;
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

    bool success=false;
    bool useSpline = (targetPtcl.Lattice.BoxBConds[0] && transformOpt == "yes");
    bool ignoreSpin = (spinOpt == "no");

    OrbitalConstraintsBase* control=0;
    if(funcOpt == "any")
    {
      app_log() << "    Using generic Jastrow Function Builder " <<endl;
      control=new AnyConstraints(targetPtcl,targetPsi);
    }
    else if(funcOpt == "pade") 
    {
      app_log() << "    Using analytic Pade Jastrow Functor " <<endl;
      control = new PadeConstraints(targetPtcl,targetPsi,ignoreSpin);
    } 
    else if(funcOpt == "rpa") 
    {
      if(useSpline) 
        control = new RPAPBCConstraints(targetPtcl,targetPsi,ignoreSpin);
      else 
        control = new RPAConstraints(targetPtcl,targetPsi,ignoreSpin);
    } 
    else if(funcOpt == "scaledpade") 
    {
      app_log() << "    Using analytic Scaled Pade Jastrow Functor " <<endl;
      control = new ScaledPadeConstraints(targetPtcl,targetPsi,ignoreSpin);
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
      app_error() << "   Failed to buld " << nameOpt << "  Jastrow function " << endl;  
      delete control;
      return false;
    }

    OrbitalBase* j2=control->createTwoBody();
    if(j2== 0) {
      app_error() << "     JastrowBuilder::addTwoBody failed to create a two body Jastrow" << endl;
      delete control;
    }

    enum {MULTIPLE=0, LONGRANGE, ONEBODY, TWOBODY, THREEBODY, FOURBODY};

    if(control->JComponent[MULTIPLE])
    {
      //create a combo orbital
      ComboOrbital* jcombo=new ComboOrbital(control);
      jcombo->Psi.push_back(j2);
      if(control->JComponent[ONEBODY] && sourceOpt != targetPtcl.getName())
      {
        app_log() << "    Adding one-body Jastrow function dependent upon two-body " << funcOpt << endl;
        map<string,ParticleSet*>::iterator pa_it(ptclPool.find(sourceOpt));
        if(pa_it != ptclPool.end()) 
        {
          OrbitalBase* j1=control->createOneBody(*((*pa_it).second));
          if(j1) jcombo->Psi.push_back(j1);
        }
      }
      //if(control->JComponent[LONGRANGE])
      //{
      //  app_log() << "    Adding long-range component of " << funcOpt << "  Jastrow function "<< endl;
      //  control->addExtra2ComboOrbital(jcombo);
      //}
      targetPsi.addOrbital(jcombo);
    }
    else
    {
      targetPsi.addOrbital(j2);
    }
    control->addOptimizables(targetPsi.VarList);
    Children.push_back(control);
    return success;
  }

  bool JastrowBuilder::addThreeBody(xmlNodePtr cur) 
  {
    if(sourceOpt == targetPtcl.getName()) 
    {
      app_warning() << "  Three-Body Jastrow Function needs a source different from " << targetPtcl.getName() << endl;
      app_warning() << "  Exit JastrowBuilder::addOneBody." << endl;
      return false;
    }

    PtclPoolType::iterator pit(ptclPool.find(sourceOpt));
    if(pit == ptclPool.end()) {
      app_error() << "     JastrowBuilder::addThreeBody requires a center. " << sourceOpt << " is invalid " << endl;
      return false;
    }

    app_log() << "  JastrowBuilder::addThreeBody source="<< sourceOpt <<  endl;
    xmlNodePtr basisPtr=NULL;
    xmlNodePtr coeffPtr=NULL;
    cur = cur->xmlChildrenNode;
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == basisset_tag) {
        basisPtr=cur;
        //call the BasisSet builder
        //basisSet = gtoBuilder->addBasisSet(cur);
      } else if(cname == "coefficient" || cname == "coefficients") {
        coeffPtr=cur;
      }
      cur=cur->next;
    }

    if(basisPtr == NULL)
    {
      app_error() << "     JastrowBuilder::addThreeBody exit. Missing <basisset/>"<< endl;
      return false;
    }

    ParticleSet* sourcePtcl=(*pit).second;
    BasisSetBuilder* basisBuilder = 
      new JastrowBasisBuilder(targetPtcl,*sourcePtcl,funcOpt,transformOpt == "yes");
    basisBuilder->put(basisPtr);

    if(typeOpt.find("Geminal") < typeOpt.size())
    {
      app_log() << "\n  creating Three-Body-Germinal Jastrow function " << endl;
      ThreeBodyGeminal* J3 = new ThreeBodyGeminal(*sourcePtcl, targetPtcl);
      J3->setBasisSet(basisBuilder->myBasisSet);
      J3->put(coeffPtr,targetPsi.VarList);
      targetPsi.addOrbital(J3);
    }

    //if(jbuilder)
    //{
    //  jbuilder->put(cur);
    //  Children.push_back(jbuilder);
    //  return true;
    //}
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
