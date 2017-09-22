//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/CorrelatedSampling/CSVMC.h"
#include "QMCDrivers/CorrelatedSampling/CSVMCUpdateAll.h"
#include "QMCDrivers/CorrelatedSampling/CSVMCUpdatePbyP.h"
#include "Estimators/CSEnergyEstimator.h"
#include "QMCDrivers/DriftOperators.h"
#if !defined(REMOVE_TRACEMANAGER)
#include "Estimators/TraceManager.h"
#else
typedef int TraceManager;
#endif

namespace qmcplusplus
{

/// Constructor.
CSVMC::CSVMC(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,WaveFunctionPool& ppool):
  QMCDriver(w,psi,h,ppool), multiEstimator(0), Mover(0)
{
  RootName = "csvmc";
  QMCType ="CSVMC";
  equilBlocks=-1;
  m_param.add(equilBlocks,"equilBlocks","int");
  QMCDriverMode.set(QMC_MULTIPLE,1);
  //Add the primary h and psi, extra H and Psi pairs will be added by QMCMain
  add_H_and_Psi(&h,&psi);
}

/** allocate internal data here before run() is called
 * @author SIMONE
 *
 * See QMCDriver::process
 */
bool CSVMC::put(xmlNodePtr q)
{
  int nPsi=H1.size();
  //for(int ipsi=0; ipsi<nPsi; ipsi++)
  //  H1[ipsi]->add2WalkerProperty(W);
  Estimators = branchEngine->getEstimatorManager();
  if(Estimators == 0)
  {
    Estimators = new EstimatorManagerBase(myComm);
    multiEstimator = new CSEnergyEstimator(H,nPsi);
    Estimators->add(multiEstimator,Estimators->MainEstimatorName);
    branchEngine->setEstimatorManager(Estimators);
  }
  H1[0]->setPrimary(true);
  for(int ipsi=1; ipsi<nPsi; ipsi++)
    H1[ipsi]->setPrimary(false);
  return true;
}

/** Run the CSVMC algorithm.
 *
 * Similar to VMC::run
 */
bool CSVMC::run()
{
  if(Mover==0)
  {
    if(QMCDriverMode[QMC_UPDATE_MODE])
    {
      app_log() << "  Using particle-by-particle update " << std::endl;
      Mover=new CSVMCUpdatePbyP(W,Psi,H,Random);
    }
    else
    {
      app_log() << "  Using walker-by-walker update " << std::endl;
      Mover=new CSVMCUpdateAll(W,Psi,H,Random);
    }
    Mover->put(qmcNode);
    Mover->Psi1=Psi1;
    Mover->H1=H1;
    Mover->multiEstimator=multiEstimator;
    Mover->resetRun(branchEngine,Estimators);
  }
  if(QMCDriverMode[QMC_UPDATE_MODE])
    Mover->initCSWalkersForPbyP(W.begin(),W.end(),equilBlocks>0);
  else
    Mover->initCSWalkers(W.begin(),W.end(), equilBlocks>0);
  Mover->startRun(nBlocks,true);
  IndexType block = 0;
  int nPsi=Psi1.size();
  do
  {
    IndexType step = 0;
    nAccept = 0;
    nReject=0;
    Mover->startBlock(nSteps);
    do
    {
      Mover->advanceWalkers(W.begin(), W.end(),false);
      step++;
      CurrentStep++;
      Estimators->accumulate(W);
    }
    while(step<nSteps);
    multiEstimator->evaluateDiff();
    //Modify Norm.
    if(block < equilBlocks)
      Mover->updateNorms();
    Mover->stopBlock();
    ++block;
    //record the current configuration
    recordBlock(block);
    //if(QMCDriverMode[QMC_UPDATE_MODE] && CurrentStep%100 == 0)
    //  Mover->updateCSWalkers(W.begin(),W.end());
  }
  while(block<nBlocks);
  Mover->stopRun();
  //finalize a qmc section
  return finalize(block);
}
}

