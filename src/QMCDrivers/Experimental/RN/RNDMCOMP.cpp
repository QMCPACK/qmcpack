//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/DMC/RNDMCOMP.h"
#include "QMCDrivers/DMC/DMCUpdatePbyP.h"
#include "QMCDrivers/DMC/DMCUpdateAll.h"
#include "QMCApp/HamiltonianPool.h"
#include "Message/Communicate.h"
#include "Message/OpenMP.h"
#include "Utilities/Timer.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Utilities/ProgressReportEngine.h"

namespace qmcplusplus
{

/// Constructor.
RNDMCOMP::RNDMCOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, TrialWaveFunction& guide, QMCHamiltonian& h, HamiltonianPool& hpool,WaveFunctionPool& ppool)
  : QMCDriver(w,psi,h,ppool), CloneManager(hpool), Guide(guide)
  , BranchInterval(-1),mover_MaxAge(-1), myRNWarmupSteps(0)
{
  RootName = "dmc";
  QMCType ="RNDMCOMP";
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  m_param.add(KillWalker,"killnode","string");
  m_param.add(BranchInterval,"branchInterval","string");
  m_param.add(mover_MaxAge,"MaxAge","double");
  m_param.add(myRNWarmupSteps,"rnwarmupsteps", "int");
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
    delete Movers[ip];
    delete estimatorClones[ip];
    delete branchClones[ip];
    estimatorClones[ip]= new EstimatorManager(*Estimators);
    estimatorClones[ip]->setCollectionMode(false);
    branchClones[ip] = new BranchEngineType(*branchEngine);
    Movers[ip] = new RNDMCUpdatePbyPFast(*wClones[ip],*wgClones[ip],*psiClones[ip],*guideClones[ip],*hClones[ip],*Rng[ip]);
    Movers[ip]->put(qmcNode);
    Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
  }
}


void RNDMCOMP::resetUpdateEngines()
{
  ReportEngine PRE("RNDMCOMP","resetUpdateEngines");
  Timer init_timer;
  makeClones(W,Psi,H);
  makeClones(W,Guide);
  if(Movers.empty())
  {
    W.loadEnsemble(wClones);
    branchEngine->initWalkerController(W,false,true);
    branchEngine->setRN(true);
    //if(QMCDriverMode[QMC_UPDATE_MODE]) W.clearAuxDataSet();
    Movers.resize(NumThreads,0);
    branchClones.resize(NumThreads,0);
    Rng.resize(NumThreads,0);
    estimatorClones.resize(NumThreads,0);
    FairDivideLow(W.getActiveWalkers(),NumThreads,wPerNode);
    bool klw=(KillWalker=="yes");
    {
      //log file
      std::ostringstream o;
      o << "  Initial partition of walkers on a node: ";
      copy(wPerNode.begin(),wPerNode.end(),std::ostream_iterator<int>(o," "));
      o << "\n";
      o << "Killing Walkers at nodes " << (useAlternate!="yes") << std::endl;
      o << "Running the released node driver."<< std::endl;
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
      Movers[ip] = new RNDMCUpdatePbyPFast(*wClones[ip],*wgClones[ip],*psiClones[ip],*guideClones[ip],*hClones[ip],*Rng[ip]);
      Movers[ip]->put(qmcNode);
      Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
      MCWalkerConfiguration::iterator wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
      Movers[ip]->initWalkersForPbyP(wit, wit_end);
    }
  }
//     branchEngine->checkParameters(W);
  if(BranchInterval<0)
    BranchInterval=1;
  {
    std::ostringstream o;
//       if (useAlternate=="yes") o << "  Using Alternate Mover"<< std::endl;
//       else o << "  Using Ceperley Mover"<< std::endl;
    o << "  BranchInterval = " << BranchInterval << "\n";
    o << "  Steps per block = " << nSteps << "\n";
    o << "  Number of blocks = " << nBlocks << "\n";
    app_log() << std::endl << o.str() << std::endl;
  }
//     RealType avg_w(0);
//     RealType n_w(0);
//     RealType logepsilon(0);
// #pragma omp parallel
//       {
//         int ip=omp_get_thread_num();
//
//         for (int step=0; step<nWarmupSteps; ++step)
//         {
//           avg_w=0;
//           n_w=0;
//           for (int prestep=0; prestep<myRNWarmupSteps; ++prestep)
//           {
//             Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
//             #pragma omp single
//             {
//               MCWalkerConfiguration::iterator wit(W.begin()), wit_end(W.end());
//               while (wit!=wit_end)
//               {
//                 avg_w += (*wit)->Weight;
//                 n_w +=1;
//                 wit++;
//               }
//             }
//             #pragma omp barrier
//            }
//            #pragma omp single
//            {
//              avg_w *= 1.0/n_w;
//              RealType w_m = avg_w/(1.0-avg_w);
//              w_m = std::log(0.5+0.5*w_m);
//              if (std::abs(w_m)>0.01)
//                logepsilon += w_m;
//            }
//            #pragma omp barrier
//            Movers[ip]->setLogEpsilon(logepsilon);
//         }
//
//       }
  app_log() << " RNDMC Engine Initialization = " << init_timer.elapsed() << " secs " << std::endl;
}

bool RNDMCOMP::run()
{
  resetUpdateEngines();
  //estimator does not need to collect data
  Estimators->setCollectionMode(true);
  Estimators->start(nBlocks);
  for(int ip=0; ip<NumThreads; ip++)
    Movers[ip]->startRun(nBlocks,false);
  Timer myclock;
  IndexType block = 0;
  IndexType updatePeriod=(QMCDriverMode[QMC_UPDATE_MODE])?Period4CheckProperties:(nBlocks+1)*nSteps;
  do // block
  {
    Estimators->startBlock(nSteps);
    for(int ip=0; ip<NumThreads; ip++)
      Movers[ip]->startBlock(nSteps);
    IndexType step = 0;
    for(IndexType step=0; step< nSteps; ++step, CurrentStep+=BranchInterval)
    {
      #pragma omp parallel for
      for(int ip=0; ip<NumThreads; ++ip)
      {
        int now=CurrentStep;
        MCWalkerConfiguration::iterator
        wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
        for(int interval = 0; interval<BranchInterval-1; ++interval,++now)
          Movers[ip]->advanceWalkers(wit,wit_end,false);
        wClones[ip]->resetCollectables();
        Movers[ip]->advanceWalkers(wit,wit_end,false);
        Movers[ip]->setReleasedNodeMultiplicity(wit,wit_end);
        if(QMCDriverMode[QMC_UPDATE_MODE] && now%updatePeriod == 0)
          Movers[ip]->updateWalkers(wit, wit_end);
      }//#pragma omp parallel
      for(int ip=1; ip<NumThreads; ++ip)
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
  }
  while(block<nBlocks && myclock.elapsed()<MaxCPUSecs);
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

