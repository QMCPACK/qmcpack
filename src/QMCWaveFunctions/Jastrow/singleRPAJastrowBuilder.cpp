//////////////////////////////////////////////////////////////////
// (c) Copyright 2008- 
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: 
//   Tel:    
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-

#include "QMCWaveFunctions/Jastrow/singleRPAJastrowBuilder.h"
#include "QMCWaveFunctions/Jastrow/LRBreakupUtilities.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "LongRange/LRRPAHandlerTemp.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"


namespace qmcplusplus {

  bool singleRPAJastrowBuilder::put(xmlNodePtr cur)
  {
    MyName="Jep";
    string rpafunc="RPA";
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
        
    if(Rs<0) {
        Rs=tlen;
    }

    if(Kc<0){ 
      Kc = 1e-4 ;
    };


    myHandler= new LRRPAHandlerTemp<EPRPABreakup<RealType>,LPQHIBasis>(targetPtcl,Kc);
    myHandler->Breakup(targetPtcl,Rs);
    
//     app_log() << "  Maximum K shell " << myHandler->MaxKshell << endl;
//     app_log() << "  Number of k vectors " << myHandler->Fk.size() << endl;
    
    
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

    JneType* J1s = new JneType (*sourcePtcl,targetPtcl);

    for(int ig=0; ig<ng; ig++) {
      J1s->addFunc(ig,nfunc);
    }

    app_log()<<" Only Short range part of E-I RPA is implemented"<<endl;
    targetPsi.addOrbital(J1s,MyName);
    return true;
  }
}
