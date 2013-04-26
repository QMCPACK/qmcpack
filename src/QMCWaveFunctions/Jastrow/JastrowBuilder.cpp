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
#include "QMCWaveFunctions/Jastrow/PadeJastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/RPAJastrow.h"
#include "QMCWaveFunctions/Jastrow/BsplineJastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/kSpaceJastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/singleRPAJastrowBuilder.h"
#if OHMMS_DIM ==3
#include "QMCWaveFunctions/Jastrow/eeI_JastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/ThreeBodyGeminal.h"
#include "QMCWaveFunctions/Jastrow/ThreeBodyBlockSparse.h"
#endif
#include "QMCWaveFunctions/Jastrow/JABBuilder.h"
#include "QMCWaveFunctions/Jastrow/JAABuilder.h"
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
    app_warning() << "  JastrowBuilder::put does not have name "<< endl;
    return false;
  }
  if(typeOpt.find("One") < typeOpt.size())
    return addOneBody(cur);
  if(typeOpt.find("Two") < typeOpt.size())
    return addTwoBody(cur);
  if(typeOpt.find("Three") < typeOpt.size())
    return addThreeBody(cur);
  if(typeOpt.find("eeI") < typeOpt.size())
    return add_eeI(cur);
  if(typeOpt.find("kSpace") < typeOpt.size())
    return addkSpace(cur);
  return false;
}

bool JastrowBuilder::addkSpace(xmlNodePtr cur)
{
  app_log() << "  JastrowBuilder::addkSpace(xmlNodePtr)" << endl;
  map<string,ParticleSet*>::iterator pa_it(ptclPool.find(sourceOpt));
  if(pa_it == ptclPool.end())
  {
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
  //use lowercase, to be handled by parser later
  tolower(funcOpt);
  if (funcOpt == "bspline" )
  {
    app_log() << "\n  Using BsplineBuilder for one-body jastrow with B-spline functions" << endl;
    BsplineJastrowBuilder jb(targetPtcl,targetPsi,*sourcePtcl);
    success=jb.put(cur);
  }
  else
    if (funcOpt == "rpa" )
    {
#if OHMMS_DIM ==3
      app_log() << "\n  Using RPA for one-body jastrow" << endl;
      singleRPAJastrowBuilder jb(targetPtcl, targetPsi, *sourcePtcl);
      success= jb.put(cur);
#else
      APP_ABORT("RPA for one-body jastrow is only available for 3D");
#endif
    }
    else
    {
      app_log() << "\n  Using JABBuilder for one-body jastrow with analytic functions" << endl;
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
                << sourceOpt << " is invalid " << endl;
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
//     OrbitalConstraintsBase* control=0;
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
  else
    if((funcOpt == "yukawa") || (funcOpt == "rpa"))
    {
      if(targetPtcl.Lattice.SuperCellEnum == SUPERCELL_OPEN)
      {
        PRE.warning("RPA is requested for an open system. Please choose other functors.");
        return false;
      }
//         control = new RPAPBCConstraints(targetPtcl,targetPsi,ignoreSpin);
      RPAJastrow* rpajastrow = new RPAJastrow(targetPtcl,targetPsi.is_manager());
      rpajastrow->put(cur);
      targetPsi.addOrbital(rpajastrow,nameOpt);
      return true;
    }
    else
      if (funcOpt == "bspline" )
      {
        BsplineJastrowBuilder bbuilder(targetPtcl,targetPsi);
        return bbuilder.put(cur);
      }
#if QMC_BUILD_LEVEL>2
      else //try other special analytic functions
      {
        app_log() << "\n  Using JAABuilder for two-body jastrow with analytic functions" << endl;
        JAABuilder jb(targetPtcl,targetPsi);
        return jb.put(cur);
      }
#endif
//     control->setReportLevel(ReportLevel);
//     success=control->put(cur);
//     if(!success) {
//       PRE.error("Failed to build "+nameOpt+" Jastrow function");
//       delete control;
//       return false;
//     }
//
//     OrbitalBase* j2=control->createTwoBody();
//     if(j2== 0)
//     {
//       APP_ABORT("JastrowBuild::addTwoBody");
//       delete control;
//       return false;
//     }
//
//     enum {MULTIPLE=0, LONGRANGE, ONEBODY, TWOBODY, THREEBODY, FOURBODY};
//
//     if(control->JComponent[MULTIPLE])
//     {
//       //create a combo orbital
//       ComboOrbital* jcombo=new ComboOrbital(control);
//       jcombo->Psi.push_back(j2);
//       if(control->JComponent[ONEBODY] && sourceOpt != targetPtcl.getName())
//       {
//         app_log() << ("Adding one-body Jastrow function dependent upon two-body " + funcOpt + " Jastrow function\n");
//         map<string,ParticleSet*>::iterator pa_it(ptclPool.find(sourceOpt));
//         if(pa_it != ptclPool.end())
//         {
//           OrbitalBase* j1=control->createOneBody(*((*pa_it).second));
//           if(j1) jcombo->Psi.push_back(j1);
//         }
//       }
//       if(control->JComponent[LONGRANGE])
//       {
//         app_log() << ("Adding  long-range component of " + funcOpt + " Jastrow function\n");
//         control->addExtra2ComboOrbital(jcombo);
//       }
//       targetPsi.addOrbital(jcombo,"J2_LR");
//     }
//     else
//     {
//       string j2name="J2_"+funcOpt;
//       targetPsi.addOrbital(j2,j2name);
//     }
//     Children.push_back(control);
  return success;
}

bool JastrowBuilder::addThreeBody(xmlNodePtr cur)
{
#if OHMMS_DIM==3
  if(sourceOpt == targetPtcl.getName())
  {
    app_warning() << "  Three-Body Jastrow Function needs a source different from " << targetPtcl.getName() << endl;
    APP_ABORT("  JastrowBuilder::addThreeBody");
    return false;
  }
  PtclPoolType::iterator pit(ptclPool.find(sourceOpt));
  if(pit == ptclPool.end())
  {
    app_error() << "     JastrowBuilder::addThreeBody requires a center. " << sourceOpt << " is invalid " << endl;
    APP_ABORT("  JastrowBuilder::addThreeBody");
    return false;
  }
  app_log() << "  JastrowBuilder::addThreeBody source="<< sourceOpt <<  endl;
  xmlNodePtr basisPtr=NULL;
  xmlNodePtr coeffPtr=NULL;
  string diagOnly("no");//default:use diagonal blocks only
  string sameblocks("yes");
  OhmmsAttributeSet tAttrib;
  tAttrib.add(diagOnly,"diagonal");
  tAttrib.add(sameblocks,"sameBlocksForGroup");
  tAttrib.put(cur);
  cur = cur->xmlChildrenNode;
  while(cur != NULL)
  {
    string cname((const char*)(cur->name));
    if(cname == basisset_tag)
    {
      basisPtr=cur;
    }
    else
      if(cname.find("oeff")<cname.size())
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
    J3->put(coeffPtr);
    J3->SameBlocksForGroup=(sameblocks == "yes");
    J3->setBlocks(basisBuilder->SizeOfBasisPerCenter);
    targetPsi.addOrbital(J3,"J3_block");
  }
  else
  {
    app_log() << "\n  creating Three-Body Jastrow function using a complete Geminal matrix." << endl;
    ThreeBodyGeminal* J3 = new ThreeBodyGeminal(*sourcePtcl, targetPtcl);
    J3->setBasisSet(basisBuilder->myBasisSet);
    J3->put(coeffPtr);
    targetPsi.addOrbital(J3,"J3_full");
  }
  delete basisBuilder;
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
