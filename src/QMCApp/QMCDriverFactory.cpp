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
#if !defined(QMCPLUSPLUS_RELEASE)
#include "QMCDrivers/VMCMultiple.h"
#include "QMCDrivers/VMCPbyPMultiple.h"
#include "QMCDrivers/ReptationMC.h"
#include "QMCDrivers/RQMCMultiple.h"
#include "QMCDrivers/RQMCMultiple.h"
#include "QMCDrivers/VMCMultipleWarp.h"
#include "QMCDrivers/RQMCMultiWarp.h"
#endif
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
  QMCDriverFactory::createQMCDriver(int curSeries, xmlNodePtr cur) {

    string what("invalid");
    string append_tag("no");
    OhmmsAttributeSet aAttrib;
    aAttrib.add(what,"method");
    aAttrib.add(append_tag,"append");
    aAttrib.put(cur);

    bool append_run = (append_tag == "yes");

    if(qmcDriver) {
      if(what != curMethod) {
        delete qmcDriver;
        qmcDriver = 0;
        //if the current qmc method is different from the previous one, append_run is set to false
        append_run = false;
      }
    }

    if(curSeries == 0) append_run = false;

    if(qmcDriver == 0) {
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

      ///////////////////////////////////////////////
      if (what == "vmc"){
        primaryH->addOperator(new ConservedEnergy,"Flux");
        qmcDriver = new VMC(*qmcSystem,*primaryPsi,*primaryH);
        curRunType = VMC_RUN;
      } else if(what == "vmc-ptcl"){
        primaryH->addOperator(new ConservedEnergy,"Flux");
        qmcDriver = new VMCParticleByParticle(*qmcSystem,*primaryPsi,*primaryH);
        curRunType = VMC_RUN;
      } else if(what == "dmc" || what == "dmc-ptcl") {
        DMCFactory fac(what,cur);
        qmcDriver = fac.create(*qmcSystem,*primaryPsi,*primaryH,*hamPool);
        curRunType = DMC_RUN;
      } else if(what == "optimize"){
        primaryH->remove("Flux");
        QMCOptimize *opt = new QMCOptimize(*qmcSystem,*primaryPsi,*primaryH);
        //opt->addConfiguration(PrevConfigFile);
        opt->setWaveFunctionNode(psiPool->getWaveFunctionNode("null"));
        qmcDriver=opt;
        curRunType = OPTIMIZE_RUN;
#if !defined(QMCPLUSPLUS_RELEASE)
      } else if(what == "vmc-multi") {
        qmcDriver = new VMCMultiple(*qmcSystem,*primaryPsi,*primaryH);
        while(targetH.size()) {
          qmcDriver->add_H_and_Psi(targetH.front(),targetPsi.front());
          targetH.pop();
          targetPsi.pop(); 
        }
        curRunType = VMC_RUN;
      } else if(what == "vmc-warp") {
        qmcDriver = new VMCMultipleWarp(*qmcSystem,*primaryPsi,*primaryH,
            *ptclPool);
        while(targetH.size()) {
          qmcDriver->add_H_and_Psi(targetH.front(),targetPsi.front());
          targetH.pop();
          targetPsi.pop(); 
        }
        curRunType = VMC_RUN;
      } else if(what == "vmc-ptcl-multi") {
        qmcDriver = new VMCPbyPMultiple(*qmcSystem,*primaryPsi,*primaryH);
        while(targetH.size()) {
          qmcDriver->add_H_and_Psi(targetH.front(),targetPsi.front());
          targetH.pop();
          targetPsi.pop(); 
        }
        curRunType = VMC_RUN;
      } else if(what == "rmc" || what == "rmc-multi") {
        qmcDriver = new RQMCMultiple(*qmcSystem,*primaryPsi,*primaryH);
        while(targetH.size()) {
          qmcDriver->add_H_and_Psi(targetH.front(),targetPsi.front());
          targetH.pop();
          targetPsi.pop(); 
        }
        curRunType = RMC_RUN;
      } else if(what == "rmc-multi-warp") {
        qmcDriver = new RQMCMultiWarp(*qmcSystem,*primaryPsi,*primaryH, *ptclPool);
        while(targetH.size()) {
          qmcDriver->add_H_and_Psi(targetH.front(),targetPsi.front());
          targetH.pop();
          targetPsi.pop(); 
        }
        curRunType = RMC_RUN;
#endif
      } else if(what == "test") {
        qmcDriver = new WaveFunctionTester(*qmcSystem,*primaryPsi,*primaryH);
        WARNMSG("Testing wavefunctions. Creating WaveFunctionTester for testing")
        curRunType = DUMMY_RUN;
      } else {
        qmcDriver = new DummyQMC(*qmcSystem,*primaryPsi,*primaryH);
        WARNMSG("Cannot termine what type of qmc to run. Creating DummyQMC for testing")
        curRunType = DUMMY_RUN;
      }
    }

    curMethod = what;
    return append_run;
  }
}
/***********************************************************************
 * $RCSfilMCMain.cpp,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
