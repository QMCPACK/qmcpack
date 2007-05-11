//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim and Simone Chiesa
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCDrivers/CorrelatedSampling/CSVMC.h"
#include "QMCDrivers/CorrelatedSampling/CSVMCUpdateAll.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/HDFWalkerIO.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "Message/Communicate.h"
#include "Estimators/CSEnergyEstimator.h"
#include "QMCDrivers/DriftOperators.h"

namespace qmcplusplus { 

  /// Constructor.
  CSVMC::CSVMC(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
    QMCDriver(w,psi,h), multiEstimator(0), Mover(0) { 
    RootName = "csvmc";
    QMCType ="CSVMC";
    equilBlocks=-1;
    m_param.add(equilBlocks,"equilBlocks","int");
    cout << "EquilBlocks " << equilBlocks << endl;

    QMCDriverMode.set(QMC_MULTIPLE,1);
    //Add the primary h and psi, extra H and Psi pairs will be added by QMCMain
    add_H_and_Psi(&h,&psi);
  }

  /** allocate internal data here before run() is called
   * @author SIMONE
   *
   * See QMCDriver::process
   */
  bool CSVMC::put(xmlNodePtr q){

    int nPsi=H1.size();

    for(int ipsi=0; ipsi<nPsi; ipsi++) 
      H1[ipsi]->add2WalkerProperty(W);

    Estimators = branchEngine->getEstimatorManager();
    if(Estimators == 0) 
    {
      Estimators = new EstimatorManager(H);
      multiEstimator = new CSEnergyEstimator(H,nPsi);

      Estimators->add(multiEstimator,Estimators->MainEstimatorName);
      branchEngine->setEstimatorManager(Estimators);
    }

    app_log() << "Number of H and Psi " << nPsi << endl;

    H1[0]->setPrimary(true);
    for(int ipsi=1; ipsi<nPsi; ipsi++) H1[ipsi]->setPrimary(false);
    return true;
  }
  
  /** Run the CSVMC algorithm.
   *
   * Similar to VMC::run 
   */
  bool CSVMC::run() {

    if(Mover==0)
    {
      Mover=new CSVMCUpdateAll(W,Psi,H,Random);
      Mover->Psi1=Psi1;
      Mover->H1=H1;
      Mover->multiEstimator=multiEstimator;
    }

    Mover->resetRun(branchEngine,Estimators);
    Mover->initCSWalkers(W.begin(),W.end(), equilBlocks>0);
    Mover->startRun(nBlocks,true);

    IndexType block = 0;
    int nPsi=Psi1.size();
    do {
      IndexType step = 0;
      nAccept = 0; nReject=0;

      Mover->startBlock(nSteps);
      do {
        Mover->advanceWalkers(W.begin(), W.end(),false);
        step++;CurrentStep++;
        Estimators->accumulate(W);
      } while(step<nSteps);

      multiEstimator->evaluateDiff();

      //Modify Norm. 
      if(block < equilBlocks) Mover->updateNorms();

      Mover->stopBlock();
      ++block;

      //record the current configuration
      recordBlock(block);

    } while(block<nBlocks);

    Mover->stopRun();

    //finalize a qmc section
    return finalize(block);
  }
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1593 $   $Date: 2007-01-04 17:23:27 -0600 (Thu, 04 Jan 2007) $
 * $Id: CSVMC.cpp 1593 2007-01-04 23:23:27Z jnkim $ 
 ***************************************************************************/
