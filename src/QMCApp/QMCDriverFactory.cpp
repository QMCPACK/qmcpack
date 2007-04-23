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
#include "QMCApp/WaveFunctionPool.h"
#include "QMCApp/HamiltonianPool.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCHamiltonians/ConservedEnergy.h"
#include "QMCDrivers/VMC/VMCFactory.h"
#include "QMCDrivers/DMC/DMCFactory.h"
#include "QMCDrivers/QMCOptimize.h"
#include "QMCDrivers/RQMCMultiple.h"
#if !defined(QMC_COMPLEX)
#include "QMCDrivers/RQMCMultiWarp.h"
#endif
#include "QMCDrivers/WaveFunctionTester.h"
#include "Utilities/OhmmsInfo.h"
#include <queue>
using namespace std;
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"

namespace qmcplusplus {

  ///initialize the static data member
  //ParticleSetPool* QMCDriverFactory::ptclPool = new ParticleSetPool;

  QMCDriverFactory::QMCDriverFactory(): qmcComm(0), qmcSystem(0), qmcDriver(0) 
  {
    ////create ParticleSetPool
    ptclPool = new ParticleSetPool;

    //create WaveFunctionPool
    psiPool = new WaveFunctionPool;
    psiPool->setParticleSetPool(ptclPool);

    //create HamiltonianPool
    hamPool = new HamiltonianPool;
    hamPool->setParticleSetPool(ptclPool);
    hamPool->setWaveFunctionPool(psiPool);
  }

  QMCDriverFactory::~QMCDriverFactory()
  {
    delete hamPool;
    delete psiPool;
    delete ptclPool;
    
    if(qmcComm) delete qmcComm;
  }

  void QMCDriverFactory::putCommunicator(xmlNodePtr cur)
  {
    //this should be done only once
    if(qmcComm) return;
    ParameterSet params;
    int nparts=1;
    params.add(nparts,"groups","int");
    params.add(nparts,"twistAngles","int");
    params.put(cur);
    if(nparts>1) 
    {
      app_log() << "  Communiator groups = " << nparts << endl;
      qmcComm = new Communicate(*OHMMS::Controller,nparts);
    }
  }

  bool QMCDriverFactory::setQMCDriver(int curSeries, xmlNodePtr cur) 
  {

    string curName((const char*)cur->name);
    string update_mode("walker");
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
    if(curName != "qmc") qmc_mode=curName;
    int nchars=qmc_mode.size();
    if(qmc_mode.find("opt") < nchars)
    {
      newRunType=OPTIMIZE_RUN;
    }
    else
    {
      if(qmc_mode.find("vmc")<nchars)
      {
        newRunType=VMC_RUN;
      }
      else if(qmc_mode.find("dmc")<nchars)
      {
        newRunType=DMC_RUN;
      }
      else if(qmc_mode.find("rmc")<nchars)
      {
        newRunType=RMC_RUN;
      }
      if(qmc_mode.find("ptcl")<nchars) WhatToDo[UPDATE_MODE]=1;
      if(qmc_mode.find("mul")<nchars) WhatToDo[MULTIPLE_MODE]=1;
      if(qmc_mode.find("warp")<nchars) WhatToDo[SPACEWARP_MODE]=1;
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

    //carryon with the existing qmcDriver
    if(qmcDriver) return append_run;

    //need to create a qmcDriver
    curRunType = newRunType;
    curQmcMode = newQmcMode;
    curQmcModeBits = WhatToDo;
    createQMCDriver(cur);

    if(qmcComm)
      qmcDriver->setCommunicator(qmcComm);
    else
      qmcDriver->setCommunicator(OHMMS::Controller);

    //branchEngine has to be transferred to a new QMCDriver
    if(branchEngine) qmcDriver->setBranchEngine(branchEngine);

    return append_run;
  }

  void QMCDriverFactory::createQMCDriver(xmlNodePtr cur) 
  {
    ///////////////////////////////////////////////
    // get primaryPsi and primaryH
    ///////////////////////////////////////////////
    TrialWaveFunction* primaryPsi= 0;
    QMCHamiltonian* primaryH=0;
    queue<TrialWaveFunction*> targetPsi;//FIFO 
    queue<QMCHamiltonian*> targetH;//FIFO

    xmlNodePtr tcur=cur->children;
    while(tcur != NULL) 
    {
      if(xmlStrEqual(tcur->name,(const xmlChar*)"qmcsystem")) 
      {
        const xmlChar* t= xmlGetProp(tcur,(const xmlChar*)"wavefunction");
        if(t != NULL) 
        {
          targetPsi.push(psiPool->getWaveFunction((const char*)t));
        } 
        else 
        {
          app_warning() << " qmcsystem does not have wavefunction. Assign 0" << endl;
          targetPsi.push(0);
        }
        t= xmlGetProp(tcur,(const xmlChar*)"hamiltonian");
        if(t != NULL) 
        {
          targetH.push(hamPool->getHamiltonian((const char*)t));
        } 
        else 
        {
          app_warning() << " qmcsystem does not have hamiltonian. Assign 0" << endl;
          targetH.push(0);
        }
      }
      tcur=tcur->next;
    }

    //mark the first targetPsi and targetH as the primaries
    if(targetH.empty()) 
    {
      primaryPsi=psiPool->getPrimary();
      primaryH=hamPool->getPrimary();
    } 
    else 
    { 
      primaryPsi=targetPsi.front(); targetPsi.pop();
      primaryH=targetH.front();targetH.pop();
    }

    //set primaryH->Primary
    primaryH->setPrimary(true);

    //flux is evaluated only with single-configuration VMC
    if(curRunType == VMC_RUN && !curQmcModeBits[MULTIPLE_MODE]) 
    {
      QMCHamiltonianBase* flux=primaryH->getHamiltonian("Flux");
      if(flux == 0) primaryH->addOperator(new ConservedEnergy,"Flux");
    } 
    else 
    {
      primaryH->remove("Flux");
    }

    //(SPACEWARP_MODE,MULTIPE_MODE,UPDATE_MODE)
    if(curRunType == VMC_RUN) 
    {
      VMCFactory fac(curQmcModeBits[UPDATE_MODE],cur);
      qmcDriver = fac.create(*qmcSystem,*primaryPsi,*primaryH,*ptclPool,*hamPool);
    } 
    else if(curRunType == DMC_RUN) 
    {
      DMCFactory fac(curQmcModeBits[UPDATE_MODE],cur);
      qmcDriver = fac.create(*qmcSystem,*primaryPsi,*primaryH,*hamPool);
    } 
    else if(curRunType == RMC_RUN) 
    {
#if defined(QMC_COMPLEX)
      qmcDriver = new RQMCMultiple(*qmcSystem,*primaryPsi,*primaryH);
#else
      if(curQmcModeBits[SPACEWARP_MODE]) 
        qmcDriver = new RQMCMultiWarp(*qmcSystem,*primaryPsi,*primaryH, *ptclPool);
      else 
        qmcDriver = new RQMCMultiple(*qmcSystem,*primaryPsi,*primaryH);
#endif
    } 
    else if(curRunType == OPTIMIZE_RUN)
    {
      QMCOptimize *opt = new QMCOptimize(*qmcSystem,*primaryPsi,*primaryH);
      opt->setWaveFunctionNode(psiPool->getWaveFunctionNode("null"));
      qmcDriver=opt;
    } 
    else 
    {
      WARNMSG("Testing wavefunctions. Creating WaveFunctionTester for testing")
      qmcDriver = new WaveFunctionTester(*qmcSystem,*primaryPsi,*primaryH);
    }

    if(curQmcModeBits[MULTIPLE_MODE]) 
    {
      while(targetH.size()) 
      {
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
