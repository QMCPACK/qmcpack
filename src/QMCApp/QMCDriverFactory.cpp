//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
/**@file QMCDriverFactory.cpp
 * @brief Implments QMCMain operators.
 */
#include "QMCApp/QMCDriverFactory.h"
#include "QMCApp/ParticleSetPool.h"
#include "QMCApp/WaveFunctionPool.h"
#include "QMCApp/HamiltonianPool.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/ConservedEnergy.h"
#include "QMCDrivers/DummyQMC.h"
#include "QMCDrivers/VMC.h"
#include "QMCDrivers/VMCParticleByParticle.h"
#include "QMCDrivers/DMC/DMCFactory.h"
#include "QMCDrivers/QMCOptimize.h"
//#if !defined(QMCPLUSPLUS_RELEASE)
#include "QMCDrivers/VMCMultiple.h"
#include "QMCDrivers/VMCPbyPMultiple.h"
#include "QMCDrivers/ReptationMC.h"
#include "QMCDrivers/RQMCMultiple.h"
#include "QMCDrivers/RQMCMultiple.h"
#include "QMCDrivers/VMCMultipleWarp.h"
#include "QMCDrivers/VMCPbyPMultipleWarp.h"
#include "QMCDrivers/RQMCMultiWarp.h"
//#endif
#include "QMCDrivers/WaveFunctionTester.h"
#include "Utilities/OhmmsInfo.h"
#include <queue>
using namespace std;
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus {

  QMCDriverFactory::QMCDriverFactory(): qmcSystem(0), qmcDriver(0) {
    //create ParticleSetPool
    ptclPool = new ParticleSetPool;

    //create WaveFunctionPool
    psiPool = new WaveFunctionPool;
    psiPool->setParticleSetPool(ptclPool);

    //create HamiltonianPool
    hamPool = new HamiltonianPool;
    hamPool->setParticleSetPool(ptclPool);
    hamPool->setWaveFunctionPool(psiPool);
  }

  QMCDriverFactory::~QMCDriverFactory(){
    delete hamPool;
    delete psiPool;
    delete ptclPool;
  }

  bool
  QMCDriverFactory::setQMCDriver(int curSeries, xmlNodePtr cur) {

    string curName((const char*)cur->name);
    string update_mode("pbyp");
    string qmc_mode("invalid");
    string multi_tag("no");
    string warp_tag("no");
    string append_tag("no");

    OhmmsAttributeSet aAttrib;
    aAttrib.add(qmc_mode,"method");
    aAttrib.add(update_mode,"move");
    aAttrib.add(multi_tag,"multiple");
    aAttrib.add(warp_tag,"warp");
    aAttrib.add(append_tag,"append");
    aAttrib.put(cur);

    bool append_run =(append_tag == "yes");

    bitset<3>  WhatToDo;
    WhatToDo[SPACEWARP_MODE]= (warp_tag == "yes");
    WhatToDo[MULTIPLE_MODE]= (multi_tag == "yes");
    WhatToDo[UPDATE_MODE]= (update_mode == "pbyp");

    QMCRunType newRunType = DUMMY_RUN;
    if(curName == "qmc") { //generic qmc node
      //overwrite the bits
      if(qmc_mode == "vmc") {
        newRunType=VMC_RUN;
        WhatToDo[UPDATE_MODE]=0;
      } else if (qmc_mode == "vmc-ptcl") {
        newRunType=VMC_RUN;
        WhatToDo[UPDATE_MODE]=1;
      } else if (qmc_mode == "vmc-multi") {
        newRunType=VMC_RUN;
        WhatToDo[MULTIPLE_MODE]=1;
        WhatToDo[UPDATE_MODE]=0;
      } else if(qmc_mode == "vmc-warp") {
        newRunType=VMC_RUN;
        WhatToDo[SPACEWARP_MODE]=1;
        WhatToDo[MULTIPLE_MODE]=1;
        WhatToDo[UPDATE_MODE]=0;
      } else if(qmc_mode == "vmc-ptcl-multi") {
        newRunType=VMC_RUN;
        WhatToDo[SPACEWARP_MODE]=1;
        WhatToDo[MULTIPLE_MODE]=1;
        WhatToDo[UPDATE_MODE]=1;
      } else if(qmc_mode == "dmc") {
        newRunType=DMC_RUN;
        WhatToDo[UPDATE_MODE]=0;
      } else if(qmc_mode == "dmc-ptcl") {
        newRunType=DMC_RUN;
        WhatToDo[UPDATE_MODE]=1;
      } else if(qmc_mode == "rmc") {
        newRunType=RMC_RUN;
      } else if(qmc_mode == "rmc-multiple") {
        newRunType=RMC_RUN;
      } else if(qmc_mode == "optimize") {
        newRunType=OPTIMIZE_RUN;
      }
    } else if(curName == "vmc") {
      newRunType=VMC_RUN;
    } else if(curName == "dmc") {
      newRunType=DMC_RUN;
    } else if(curName == "rmc") {
      newRunType=RMC_RUN;
    } else if(curName == "optimize") {
      newRunType=OPTIMIZE_RUN;
    } 

    unsigned long newQmcMode=WhatToDo.to_ulong();
   
    //initialize to 0
    QMCDriver::BranchEngineType* branchEngine=0;

    if(qmcDriver) {
      if(newRunType != curRunType || newQmcMode != curQmcMode) {
        //copy the pointer of the BranchEngine 
        branchEngine=qmcDriver->getBranchEngine();
        //remove the qmcDriver
        delete qmcDriver;
        //set to 0 so that a new driver is created
        qmcDriver = 0;
        //if the current qmc method is different from the previous one, append_run is set to false
        append_run = false;
      } else { 
        app_log() << "  Reusing " << qmcDriver->getEngineName() << endl;
      }
    }

    if(curSeries == 0) append_run = false;

    //carryon
    if(qmcDriver) return append_run;

    //need to create one
    curRunType = newRunType;
    curQmcMode = newQmcMode;
    curQmcModeBits = WhatToDo;
    createQMCDriver(cur);

    //branchEngine has to be transferred to a new QMCDriver
    if(branchEngine) qmcDriver->setBranchEngine(branchEngine);

    return append_run;
  }

  void QMCDriverFactory::createQMCDriver(xmlNodePtr cur) {

    ///////////////////////////////////////////////
    // get primaryPsi and primaryH
    ///////////////////////////////////////////////
    TrialWaveFunction* primaryPsi= 0;
    QMCHamiltonian* primaryH=0;
    queue<TrialWaveFunction*> targetPsi;//FIFO 
    queue<QMCHamiltonian*> targetH;//FIFO

    xmlNodePtr tcur=cur->children;
    while(tcur != NULL) {
      if(xmlStrEqual(tcur->name,(const xmlChar*)"qmcsystem")) {
        const xmlChar* t= xmlGetProp(tcur,(const xmlChar*)"wavefunction");
        targetPsi.push(psiPool->getWaveFunction((const char*)t));
        t= xmlGetProp(tcur,(const xmlChar*)"hamiltonian");
        targetH.push(hamPool->getHamiltonian((const char*)t));
      }
      tcur=tcur->next;
    }

    //mark the first targetPsi and targetH as the primaries
    if(targetH.empty()) {
      primaryPsi=psiPool->getPrimary();
      primaryH=hamPool->getPrimary();
    } else { 
      primaryPsi=targetPsi.front(); targetPsi.pop();
      primaryH=targetH.front();targetH.pop();
    }

    //set primaryH->Primary
    primaryH->setPrimary(true);

    //flux is evaluated only with single-configuration VMC
    if(curRunType == VMC_RUN && !curQmcModeBits[MULTIPLE_MODE]) {
      QMCHamiltonianBase* flux=primaryH->getHamiltonian("Flux");
      if(flux == 0) primaryH->addOperator(new ConservedEnergy,"Flux");
    } else {
      primaryH->remove("Flux");
    }

    //(SPACEWARP_MODE,MULTIPE_MODE,UPDATE_MODE)
    if(curRunType == VMC_RUN) {
      if(curQmcMode == 0) {//(0,0,0)
        qmcDriver = new VMC(*qmcSystem,*primaryPsi,*primaryH);
      } else if(curQmcMode == 1) {//(0,0,1)
        qmcDriver = new VMCParticleByParticle(*qmcSystem,*primaryPsi,*primaryH);
      } else if(curQmcMode == 2) {//(0,1,0)
        qmcDriver = new VMCMultiple(*qmcSystem,*primaryPsi,*primaryH);
      } else if(curQmcMode == 3) {//(0,1,1)
        qmcDriver = new VMCPbyPMultiple(*qmcSystem,*primaryPsi,*primaryH);
      } else if(curQmcMode == 6) {//(1,1,0)
        qmcDriver = new VMCMultipleWarp(*qmcSystem,*primaryPsi,*primaryH, *ptclPool);
      } else if(curQmcMode == 7) {//(1,1,1)
        qmcDriver = new VMCPbyPMultipleWarp(*qmcSystem,*primaryPsi,*primaryH, *ptclPool);
      }
    } else if(curRunType == DMC_RUN) {
      DMCFactory fac(curQmcModeBits[UPDATE_MODE],cur);
      qmcDriver = fac.create(*qmcSystem,*primaryPsi,*primaryH,*hamPool);
    } else if(curRunType == RMC_RUN) {
      if(curQmcModeBits[SPACEWARP_MODE]) {
        qmcDriver = new RQMCMultiWarp(*qmcSystem,*primaryPsi,*primaryH, *ptclPool);
      } else {
        qmcDriver = new RQMCMultiple(*qmcSystem,*primaryPsi,*primaryH);
      }
    } else if(curRunType == OPTIMIZE_RUN) {
      QMCOptimize *opt = new QMCOptimize(*qmcSystem,*primaryPsi,*primaryH);
      //opt->addConfiguration(PrevConfigFile);
      opt->setWaveFunctionNode(psiPool->getWaveFunctionNode("null"));
      qmcDriver=opt;
    } else {
      WARNMSG("Testing wavefunctions. Creating WaveFunctionTester for testing")
      qmcDriver = new WaveFunctionTester(*qmcSystem,*primaryPsi,*primaryH);
      //} else {
      //  WARNMSG("Cannot termine what type of qmc to run. Creating DummyQMC for testing")
      //  qmcDriver = new DummyQMC(*qmcSystem,*primaryPsi,*primaryH);
    }

    if(curQmcModeBits[MULTIPLE_MODE]) {
      while(targetH.size()) {
        qmcDriver->add_H_and_Psi(targetH.front(),targetPsi.front());
        targetH.pop();
        targetPsi.pop(); 
      }
    }
  }

}
/***********************************************************************
 * $RCSfilMCMain.cpp,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
