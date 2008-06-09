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
#include "QMCWaveFunctions/Jastrow/BsplineJastrowBuilder.h"
#if OHMMS_DIM ==3
#include "QMCWaveFunctions/Jastrow/ThreeBodyGeminal.h"
#include "QMCWaveFunctions/Jastrow/ThreeBodyBlockSparse.h"
#include "QMCWaveFunctions/Jastrow/kSpaceJastrowBuilder.h"
#endif
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

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

    if(typeOpt.find("kSpace") < typeOpt.size())
      return addkSpace(cur);

    return false;
  }

  bool JastrowBuilder::addkSpace(xmlNodePtr cur)
  {
    app_log() << "  JastrowBuilder::addkSpace(xmlNodePtr)" << endl;
    map<string,ParticleSet*>::iterator pa_it(ptclPool.find(sourceOpt));
    if(pa_it == ptclPool.end()) {
      app_warning() << "  JastrowBuilder::addkSpace failed. " 
		    << sourceOpt << " does not exist" << endl;
      return false;
    }
    ParticleSet* sourcePtcl= (*pa_it).second;
    app_log() << "\n  Using kSpaceJastrowBuilder for reciprocal-space Jastrows" << endl;
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

    map<string,ParticleSet*>::iterator pa_it(ptclPool.find(sourceOpt));
    if(pa_it == ptclPool.end()) 
    {
      PRE.warning("JastrowBuilder::addOneBody failed. "+sourceOpt+" does not exist.");
      return false;
    }

    bool success=false;
    ParticleSet* sourcePtcl= (*pa_it).second;
    if(funcOpt == "any")// || funcOpt == "poly")
    {
      OrbitalConstraintsBase* control=0;
      control = new AnyConstraints(targetPtcl,targetPsi);
      control->setReportLevel(ReportLevel);
      //if(funcOpt == "any")
      //  control = new AnyConstraints(targetPtcl,targetPsi);
      //else
      //  control = new PolyConstraints(targetPtcl,targetPsi,true);
      control->put(cur);
      OrbitalBase* j=control->createOneBody(*sourcePtcl);
      if(j)
      {
        control->addOptimizables(targetPsi.VarList);
        targetPsi.addOrbital(j,"J1_WM");
        Children.push_back(control);
        success=true;
      }
      else
      {
        delete control;
      }
    }
    else if (funcOpt == "Bspline" ) 
    {
      app_log() << "\n  Using BsplineBuilder for one-body jastrow with B-spline functions" << endl;
      OrbitalBuilderBase* jb = new BsplineJastrowBuilder (targetPtcl, targetPsi, *sourcePtcl);
      jb->setReportLevel(ReportLevel);
      Children.push_back(jb);
      success= jb->put(cur);
    }
    else
    {
      app_log() << "\n  Using JABBuilder for one-body jastrow with analytic functions" << endl;
      OrbitalBuilderBase* jb = new JABBuilder(targetPtcl,targetPsi,ptclPool);
      jb->setReportLevel(ReportLevel);
      Children.push_back(jb);
      success=jb->put(cur);
    }

    return success;
  }

  bool JastrowBuilder::addTwoBody(xmlNodePtr cur) 
  {
    ReportEngine PRE(ClassName,"addTwoBody(xmlNodePtr)");

    bool success=false;
    bool useSpline = (targetPtcl.Lattice.BoxBConds[0] && transformOpt == "yes");
    bool ignoreSpin = (spinOpt == "no");

    OrbitalConstraintsBase* control=0;
    if(funcOpt == "any")
    {
      control=new AnyConstraints(targetPtcl,targetPsi);
    }
    else if(funcOpt == "pade") 
    {
      if (targetPtcl.Lattice.SuperCellEnum != SUPERCELL_OPEN) {
        PRE.warning("Pade Jastrow is requested for a periodic system. Please choose other functors.");
	return false;
      }	
      control = new PadeConstraints(targetPtcl,targetPsi,ignoreSpin);
    } 
    else if((funcOpt == "Yukawa") || (funcOpt == "RPA")) 
    {
      if(targetPtcl.Lattice.SuperCellEnum == SUPERCELL_OPEN)
      {
        PRE.warning("RPA is requested for an open system. Please choose other functors.");
        return false;
      }
      else 
        control = new RPAPBCConstraints(targetPtcl,targetPsi,ignoreSpin);
    } 
    //else if(funcOpt ==  "poly")
    //{
    //  app_log() << "    Using analytic Polynomial expansion Jastrow Functor " <<endl;
    //  control = new PolyConstraints(targetPtcl,targetPsi,ignoreSpin);
    //}
    else if(funcOpt == "scaledpade") 
    {
      control = new ScaledPadeConstraints(targetPtcl,targetPsi,ignoreSpin);
    } 
    else if (funcOpt == "Bspline" ) 
    {
      OrbitalBuilderBase* sBuilder = new BsplineJastrowBuilder (targetPtcl, targetPsi);
      Children.push_back(sBuilder);
      success=sBuilder->put(cur);

      return success;
    }
    else //try analytic functions
    {
      OrbitalBuilderBase* jb=new JAABuilder(targetPtcl,targetPsi);
      jb->setReportLevel(ReportLevel);
      Children.push_back(jb);
      success=jb->put(cur);
      return success;
    }

    control->setReportLevel(ReportLevel);
    success=control->put(cur);
    if(!success) {
      PRE.error("Failed to build "+nameOpt+" Jastro function");
      delete control;
      return false;
    }

    OrbitalBase* j2=control->createTwoBody();
    if(j2== 0) {
      PRE.error("No two body Jastrow is added.");
      delete control;
      return false;
    }

    enum {MULTIPLE=0, LONGRANGE, ONEBODY, TWOBODY, THREEBODY, FOURBODY};

    if(control->JComponent[MULTIPLE])
    {
      //create a combo orbital
      ComboOrbital* jcombo=new ComboOrbital(control);
      jcombo->Psi.push_back(j2);
      if(control->JComponent[ONEBODY] && sourceOpt != targetPtcl.getName())
      {
        app_log() << ("Adding one-body Jastrow function dependent upon two-body " + funcOpt + " Jastrow function\n");
        map<string,ParticleSet*>::iterator pa_it(ptclPool.find(sourceOpt));
        if(pa_it != ptclPool.end()) 
        {
          OrbitalBase* j1=control->createOneBody(*((*pa_it).second));
          if(j1) jcombo->Psi.push_back(j1);
        }
      }
      if(control->JComponent[LONGRANGE])
      {
        app_log() << ("Adding  long-range component of " + funcOpt + " Jastrow function\n");
        control->addExtra2ComboOrbital(jcombo);
      }
      targetPsi.addOrbital(jcombo,"J2_LR");
    }
    else
    {
      string j2name="J2_"+funcOpt;
      targetPsi.addOrbital(j2,j2name);
    }
    control->addOptimizables(targetPsi.VarList);
    Children.push_back(control);

    return success;
  }

  bool JastrowBuilder::addThreeBody(xmlNodePtr cur) 
  {
#if OHMMS_DIM==3
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
    string diagOnly("no");
    while(cur != NULL) {
      string cname((const char*)(cur->name));
      if(cname == basisset_tag) 
      {
        basisPtr=cur;
      } 
      else if(cname.find("coeff")<cname.size())
      {
        coeffPtr=cur;
        OhmmsAttributeSet oAttrib;
        oAttrib.add(diagOnly,"diagonal");
        oAttrib.put(cur);
      }
      cur=cur->next;
    }

    if(basisPtr == NULL)
    {
      app_error() << "     JastrowBuilder::addThreeBody exit. Missing <basisset/>"<< endl;
      return false;
    }

    ParticleSet* sourcePtcl=(*pit).second;
    JastrowBasisBuilder* basisBuilder = 
      new JastrowBasisBuilder(targetPtcl,*sourcePtcl,funcOpt,transformOpt == "yes");
    basisBuilder->put(basisPtr);

    if(diagOnly == "yes")
    {
      app_log() << "\n  creating Three-Body Jastrow function using only diagnoal blocks." << endl;
      ThreeBodyBlockSparse* J3 = new ThreeBodyBlockSparse(*sourcePtcl, targetPtcl);
      J3->setBasisSet(basisBuilder->myBasisSet);
      J3->put(coeffPtr,targetPsi.VarList);
      J3->setBlocks(basisBuilder->SizeOfBasisPerCenter);
      targetPsi.addOrbital(J3,"J3_block");
    } 
    else
    {
      app_log() << "\n  creating Three-Body Jastrow function using a complete Geminal matrix." << endl;
      ThreeBodyGeminal* J3 = new ThreeBodyGeminal(*sourcePtcl, targetPtcl);
      J3->setBasisSet(basisBuilder->myBasisSet);
      J3->put(coeffPtr,targetPsi.VarList);
      targetPsi.addOrbital(J3,"J3_full");
    }
#else
    app_error() << "  Three-body Jastrow function is not supported for DIM != 3." << endl;
//#error "  Three-body Jastrow is disabled for QMC_DIM != 3\n "
#endif

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
