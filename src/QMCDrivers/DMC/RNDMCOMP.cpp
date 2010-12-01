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
#include "QMCDrivers/DMC/RNDMCOMP.h"
#include "QMCDrivers/DMC/DMCUpdatePbyP.h" 
#include "QMCDrivers/DMC/DMCUpdateAll.h"
#include "QMCApp/HamiltonianPool.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"
#include "Utilities/Timer.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus { 

  /// Constructor.
  RNDMCOMP::RNDMCOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, HamiltonianPool& hpool,WaveFunctionPool& ppool)
    : QMCDriver(w,psi,h,ppool), CloneManager(hpool), useAlternate("no")
    , BranchInterval(-1),mover_MaxAge(-1)
    {
      RootName = "dmc";
      QMCType ="RNDMCOMP";

      QMCDriverMode.set(QMC_UPDATE_MODE,1);

      m_param.add(KillWalker,"killnode","string");
      m_param.add(BranchInterval,"branchInterval","string");
      m_param.add(mover_MaxAge,"MaxAge","double");
      m_param.add(useAlternate,"alternate", "string");
    }
    
  void RNDMCOMP::resetComponents(xmlNodePtr cur)
  {
    bool klw=(KillWalker=="yes");
    
    m_param.put(cur);
    put(cur);
    branchEngine->setRN(true);
    branchEngine->resetRun(cur);
    Estimators->reset();
    
    
#pragma omp parallel for
      for(int ip=0; ip<NumThreads; ++ip)
      {
        delete Movers[ip]; delete estimatorClones[ip]; delete branchClones[ip];
        estimatorClones[ip]= new EstimatorManager(*Estimators);
        estimatorClones[ip]->setCollectionMode(false);

        branchClones[ip] = new BranchEngineType(*branchEngine);
        if (useAlternate=="yes")
          Movers[ip] = new RNDMCUpdatePbyPAlternate(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]); 
        else
          Movers[ip] = new RNDMCUpdatePbyPCeperley(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]); 
        Movers[ip]->put(qmcNode);
        Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);        
      }
    
  }


  void RNDMCOMP::resetUpdateEngines() 
  {

    ReportEngine PRE("RNDMCOMP","resetUpdateEngines");
 
    Timer init_timer;

    makeClones(W,Psi,H);

    if(Movers.empty()) 
    {
      //load walkers
      W.loadEnsemble();
      for(int ip=1;ip<NumThreads;++ip) wClones[ip]->loadEnsemble(W);

      if (useAlternate=="yes")
        branchEngine->initWalkerController(W,false,false);
        //branchEngine->initWalkerController(W,Tau,false,false);
      else
        branchEngine->initWalkerController(W,false,true);
        //branchEngine->initWalkerController(W,Tau,false,true);
      branchEngine->setRN(true);

      //if(QMCDriverMode[QMC_UPDATE_MODE]) W.clearAuxDataSet();
      Movers.resize(NumThreads,0);
      branchClones.resize(NumThreads,0);
      Rng.resize(NumThreads,0);
      estimatorClones.resize(NumThreads,0);
      FairDivideLow(W.getActiveWalkers(),NumThreads,wPerNode);
      bool klw=(KillWalker=="yes");

      {//log file
        ostringstream o;
        o << "  Initial partition of walkers on a node: ";
        std::copy(wPerNode.begin(),wPerNode.end(),ostream_iterator<int>(o," "));
        o << "\n";
        o << "Killing Walkers at nodes " << (useAlternate!="yes") <<endl;

        o << "Running the released node driver."<<endl;

        app_log() << o.str();
      }

#pragma omp parallel for
      for(int ip=0; ip<NumThreads; ++ip)
      {
        estimatorClones[ip]= new EstimatorManager(*Estimators);
        estimatorClones[ip]->setCollectionMode(false);

        Rng[ip]=new RandomGenerator_t(*RandomNumberControl::Children[ip]);
        hClones[ip]->setRandomGenerator(Rng[ip]);

        branchClones[ip] = new BranchEngineType(*branchEngine);
//         Movers[ip] = new RNDMCUpdatePbyPWithRejectionFast(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]); 
        
        if (useAlternate=="yes")
          Movers[ip] = new RNDMCUpdatePbyPAlternate(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);  
        else
          Movers[ip] = new RNDMCUpdatePbyPCeperley(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]); 
        Movers[ip]->put(qmcNode);
        Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);   
        MCWalkerConfiguration::iterator wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
          Movers[ip]->initWalkersForPbyP(wit, wit_end);
        }
      } 
//     branchEngine->checkParameters(W);
    if(BranchInterval<0) BranchInterval=1;
    {
      ostringstream o;
      if (useAlternate=="yes") o << "  Using Alternate Mover"<<endl;
      else o << "  Using Ceperley Mover"<<endl;
      o << "  BranchInterval = " << BranchInterval << "\n";
      o << "  Steps per block = " << nSteps << "\n";
      o << "  Number of blocks = " << nBlocks << "\n";
      app_log() << endl << o.str() << endl;
    }
        
    app_log() << " RNDMC Engine Initialization = " << init_timer.elapsed() << " secs " << endl;
  }

  bool RNDMCOMP::run() {

    resetUpdateEngines();
    //estimator does not need to collect data
    Estimators->setCollectionMode(true);
    Estimators->start(nBlocks);
    for(int ip=0; ip<NumThreads; ip++) Movers[ip]->startRun(nBlocks,false);

    Timer myclock;
    IndexType block = 0;
    IndexType updatePeriod=(QMCDriverMode[QMC_UPDATE_MODE])?Period4CheckProperties:(nBlocks+1)*nSteps;

    do // block
    {
      Estimators->startBlock(nSteps);
      for(int ip=0; ip<NumThreads; ip++) Movers[ip]->startBlock(nSteps);
      IndexType step = 0;

      for(IndexType step=0; step< nSteps; ++step, CurrentStep+=BranchInterval)
      {
#pragma omp parallel for
        for(int ip=0; ip<NumThreads; ++ip)
        {
          int now=CurrentStep;
          MCWalkerConfiguration::iterator 
            wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);

          for(int interval = 0;interval<BranchInterval-1; ++interval,++now)
            Movers[ip]->advanceWalkers(wit,wit_end,false);

          wClones[ip]->resetCollectables();
          Movers[ip]->advanceWalkers(wit,wit_end,false);

          Movers[ip]->setReleasedNodeMultiplicity(wit,wit_end);

          if(QMCDriverMode[QMC_UPDATE_MODE] && now%updatePeriod == 0) Movers[ip]->updateWalkers(wit, wit_end);
        }//#pragma omp parallel
        
        for(int ip=1;ip<NumThreads; ++ip)
        {
          for(int j=0; j<W.Collectables.size(); ++j)
            W.Collectables[j]+=wClones[ip]->Collectables[j];
        } 
        
        branchEngine->branch(CurrentStep,W, branchClones);
        FairDivideLow(W.getActiveWalkers(),NumThreads,wPerNode);
      }
      
//       #pragma omp parallel for
//       for(int ip=0; ip<NumThreads; ++ip)
//       {
//         MCWalkerConfiguration::iterator wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
//         Movers[ip]->resetbuffers(wit, wit_end);
//       }
      
      Estimators->stopBlock(acceptRatio());
      block++;
      recordBlock(block);

    } while(block<nBlocks && myclock.elapsed()<MaxCPUSecs);

    //for(int ip=0; ip<NumThreads; ip++) Movers[ip]->stopRun();
    for(int ip=0; ip<NumThreads; ip++) 
      *(RandomNumberControl::Children[ip])=*(Rng[ip]);

    Estimators->stop();
    return finalize(block);
  }
  bool 
  RNDMCOMP::put(xmlNodePtr q)
  { 
    //nothing to do
    return true;
  }
}

/***************************************************************************
 * $RCSfile: RNDMCOMP.cpp,v $   $Author: jnkim $
 * $Revision: 1620 $   $Date: 2007-01-14 18:12:23 -0600 (Sun, 14 Jan 2007) $
 * $Id: RNDMCOMP.cpp 1620 2007-01-15 00:12:23Z jnkim $ 
 ***************************************************************************/
