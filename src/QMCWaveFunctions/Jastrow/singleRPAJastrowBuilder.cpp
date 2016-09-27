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
    
    

#include "QMCWaveFunctions/Jastrow/singleRPAJastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/LRBreakupUtilities.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "LongRange/LRRPAHandlerTemp.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"


namespace qmcplusplus
{

bool singleRPAJastrowBuilder::put(xmlNodePtr cur)
{
  return put( cur, 1);
}

bool singleRPAJastrowBuilder::put(xmlNodePtr cur, int addOrbital)
{
  MyName="Jep";
  std::string rpafunc="RPA";
  OhmmsAttributeSet a;
  a.add(MyName,"name");
  a.add(rpafunc,"function");
  a.put(cur);
  ParameterSet params;
  RealType Rs(-1.0);
  RealType Kc(-1.0);
  params.add(Rs,"rs","double");
  params.add(Kc,"kc","double");
  params.put(cur);
  if(Rs<0)
  {
    Rs=tlen;
  }
  if(Kc<0)
  {
    Kc = 1e-6 ;
  };
  if (rpafunc=="RPA")
  {
    myHandler= new LRRPAHandlerTemp<EPRPABreakup<RealType>,LPQHIBasis>(targetPtcl,Kc);
    app_log()<<"  using e-p RPA"<< std::endl;
  }
  else
    if (rpafunc=="dRPA")
    {
      myHandler= new LRRPAHandlerTemp<derivEPRPABreakup<RealType>,LPQHIBasis>(targetPtcl,Kc);
      app_log()<<"  using e-p derivRPA"<< std::endl;
    }
  myHandler->Breakup(targetPtcl,Rs);
//     app_log() << "  Maximum K shell " << myHandler->MaxKshell << std::endl;
//     app_log() << "  Number of k vectors " << myHandler->Fk.size() << std::endl;
  //Add short range part
  Rcut = myHandler->get_rc()-0.1;
  GridType* myGrid = new GridType;
  int npts=static_cast<int>(Rcut/0.01)+1;
  myGrid->set(0,Rcut,npts);
  //create the numerical functor
  nfunc = new FuncType;
  SRA = new ShortRangePartAdapter<RealType>(myHandler);
  SRA->setRmax(Rcut);
  nfunc->initialize(SRA, myGrid);
  J1s = new JneType (*sourcePtcl,targetPtcl);
  for(int ig=0; ig<ng; ig++)
  {
    J1s->addFunc(ig,nfunc);
  }
  app_log()<<" Only Short range part of E-I RPA is implemented"<< std::endl;
  if (addOrbital)
    targetPsi.addOrbital(J1s,MyName);
  return true;
}

OrbitalBase* singleRPAJastrowBuilder::getOrbital()
{
//     RealType Rs=tlen;
//     RealType Kc = 1e-5 ;
//     myHandler= new LRRPAHandlerTemp<derivEPRPABreakup<RealType>,LPQHIBasis>(targetPtcl,Kc);
//     myHandler->Breakup(targetPtcl,Rs);
//
//     //Add short range part
//     Rcut = myHandler->get_rc()-0.1;
//     GridType* myGrid = new GridType;
//     int npts=static_cast<int>(Rcut/0.01)+1;
//     myGrid->set(0,Rcut,npts);
//
//       //create the numerical functor
//     nfunc = new FuncType;
//     SRA = new ShortRangePartAdapter<RealType>(myHandler);
//     SRA->setRmax(Rcut);
//     nfunc->initialize(SRA, myGrid);
//
//     J1s = new JneType (*sourcePtcl,targetPtcl);
//
//     for(int ig=0; ig<ng; ig++) {
//       J1s->addFunc(ig,nfunc);
//     }
  if (J1s==0)
    app_log()<<"  ERROR!! Must singleRPAJastrowBuilder::put() first!"<< std::endl;
  return J1s;
}
}
