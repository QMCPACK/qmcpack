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
#include "QMCDrivers/DMC/DMCOMPOPT.h"
#include "QMCDrivers/DMC/DMCUpdatePbyP.h"
//#include "QMCDrivers/DMC/DMCNonLocalUpdate.h"
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
DMCOMPOPT::DMCOMPOPT(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, HamiltonianPool& hpool,WaveFunctionPool& ppool)
  : QMCDriver(w,psi,h,ppool), CloneManager(hpool)
  , KillNodeCrossing(0) ,Reconfiguration("no"), BenchMarkRun("no"), UseFastGrad("yes")
  , BranchInterval(-1),mover_MaxAge(-1), wlen(10), firsttime(true), printderivs("no")
{
  RootName = "opt";
  QMCType ="DMCOPT";
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  m_param.add(KillWalker,"killnode","string");
  m_param.add(BenchMarkRun,"benchmark","string");
  m_param.add(Reconfiguration,"reconfiguration","string");
  m_param.add(BranchInterval,"branchInterval","string");
  m_param.add(NonLocalMove,"nonlocalmove","string");
  m_param.add(NonLocalMove,"nonlocalmoves","string");
  m_param.add(mover_MaxAge,"MaxAge","double");
  m_param.add(UseFastGrad,"fastgrad", "string");
  m_param.add(nTargetSamples,"samples","int");
  m_param.add(printderivs,"printderivs","string");
  m_param.add(wlen,"wlen","int");
}

//   void DMCOMPOPT::resetComponents(xmlNodePtr cur)
//   {
//     m_param.put(cur);
//     put(cur);
// //     branchEngine->resetRun(cur);
// //     branchEngine->checkParameters(W);
//
//
//     //delete Movers[0];
//     for(int ip=0; ip<NumThreads; ++ip)
//     {
//       delete Movers[ip]; delete estimatorClones[ip]; delete branchClones[ip];
//       estimatorClones[ip]= new EstimatorManager(*Estimators);
//       estimatorClones[ip]->setCollectionMode(false);
//       branchClones[ip] = new BranchEngineType(*branchEngine);
//     }
//
// #pragma omp parallel for
//       for(int ip=0; ip<NumThreads; ++ip)
//       {
//         if(QMCDriverMode[QMC_UPDATE_MODE])
//         {
//           if(UseFastGrad == "yes")
//             Movers[ip] = new DMCUpdatePbyPWithRejectionFast(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
//           else
//             Movers[ip] = new DMCUpdatePbyPWithRejection(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
//           Movers[ip]->put(cur);
//           Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
//           Movers[ip]->initWalkersForPbyP(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
//         }
//         else
//         {
//           if(KillNodeCrossing)
//             Movers[ip] = new DMCUpdateAllWithKill(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
//           else
//             Movers[ip] = new DMCUpdateAllWithRejection(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
//           Movers[ip]->put(cur);
//           Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
//           Movers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
//         }
//       }
//       MCWalkerConfiguration::iterator wit(W.begin()), wit_end(W.end());
//       while(wit!=wit_end)
//       {
//         (**wit).resetPropertyHistory();
//         wit++;
//       }
//
//   }


void DMCOMPOPT::resetUpdateEngines()
{
  ReportEngine PRE("DMCOMPOPT","resetUpdateEngines");
  resetComponents(qmcNode);
  bool fixW = (Reconfiguration == "yes");
  Timer init_timer;
//     HACK HACK HACK
// This is so ugly it's probably a crime. It must be fixed.
  if(firsttime)
  {
//       for(int ip=1; ip<NumThreads; ++ip)
//       {
//         delete wClones[ip];
//         delete psiClones[ip];
//         delete hClones[ip];
//       }
//       wClones.resize(0);
    wClones.clear();
    psiClones.clear();
    hClones.clear();
    firsttime=false;
  }
  for (int ip=0; ip<NumThreads; ++ip)
  {
    opt_variables_type dummy;
    psiClones[ip]->checkInVariables(dummy);
    dummy.resetIndex();
    psiClones[ip]->checkOutVariables(dummy);
    dummyOptVars.push_back(dummy);
  }
  NumOptimizables=dummyOptVars[0].size();
  resizeForOpt(NumOptimizables);
  makeClones(W,Psi,H);
  if(Movers.empty())
  {
//       //load walkers
    Eindx = W.addPropertyHistory(wlen);
    W.loadEnsemble(wClones);
    for(int ip=1; ip<NumThreads; ++ip)
      wClones[ip]->addPropertyHistory(wlen);
//       m_param.put(qmcNode);
//       put(qmcNode);
//       //app_log()<<"DMCOMPOPT::resetComponents"<<endl;
//       Estimators->reset();
//       branchEngine->resetRun(qmcNode);
//       branchEngine->checkParameters(W);
    branchEngine->initWalkerController(W,fixW,false);
    //if(QMCDriverMode[QMC_UPDATE_MODE]) W.clearAuxDataSet();
    Movers.resize(NumThreads,0);
    branchClones.resize(NumThreads,0);
    Rng.resize(NumThreads,0);
    estimatorClones.resize(NumThreads,0);
    FairDivideLow(W.getActiveWalkers(),NumThreads,wPerNode);
    {
      //log file
      ostringstream o;
      o << "  Initial partition of walkers on a node: ";
      std::copy(wPerNode.begin(),wPerNode.end(),ostream_iterator<int>(o," "));
      o << "\n";
      if(QMCDriverMode[QMC_UPDATE_MODE])
        o << "  Updates by particle-by-particle moves";
      else
        o << "  Updates by walker moves";
      if(UseFastGrad == "yes")
        o << " using fast gradient version ";
      else
        o << " using full-ratio version ";
      if(KillNodeCrossing)
        o << "\n  Walkers are killed when a node crossing is detected";
      else
        o << "\n  DMC moves are rejected when a node crossing is detected";
      app_log() << o.str() << endl;
    }
    #pragma omp parallel for
    for(int ip=0; ip<NumThreads; ++ip)
    {
      estimatorClones[ip]= new EstimatorManager(*Estimators);
      estimatorClones[ip]->setCollectionMode(false);
      Rng[ip]=new RandomGenerator_t(*RandomNumberControl::Children[ip]);
      hClones[ip]->setRandomGenerator(Rng[ip]);
      branchClones[ip] = new BranchEngineType(*branchEngine);
      if(QMCDriverMode[QMC_UPDATE_MODE])
      {
        if(UseFastGrad == "yes")
          Movers[ip] = new DMCUpdatePbyPWithRejectionFast(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
        else
          Movers[ip] = new DMCUpdatePbyPWithRejection(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
//           Movers[ip]->put(qmcNode);
//           Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
//           Movers[ip]->initWalkersForPbyP(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
      }
      else
      {
        if(KillNodeCrossing)
          Movers[ip] = new DMCUpdateAllWithKill(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
        else
          Movers[ip] = new DMCUpdateAllWithRejection(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
//           Movers[ip]->put(qmcNode);
//           Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
//           Movers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
      }
    }
//       #pragma omp parallel for
//       for (int ip=0; ip<NumThreads; ++ip)
//       {
//         MCWalkerConfiguration::iterator wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
//         while (wit!=wit_end)
//         {
//           Walker_t& thisWalker(**wit);
//           Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
//           wClones[ip]->loadWalker(thisWalker,true);
//           psiClones[ip]->copyFromBuffer(*wClones[ip],w_buffer);
//           psiClones[ip]->updateBuffer(*wClones[ip],w_buffer,true);
//           wit++;
//         }
//       }
  }
  #pragma omp parallel for
  for(int ip=0; ip<NumThreads; ++ip)
  {
    if(QMCDriverMode[QMC_UPDATE_MODE])
    {
      Movers[ip]->put(qmcNode);
      Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
      Movers[ip]->initWalkersForPbyP(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
    }
    else
    {
      Movers[ip]->put(qmcNode);
      Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
      Movers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
    }
  }
  //       adding weight index
  MCWalkerConfiguration::iterator wit(W.begin()), wit_end(W.end());
  (**wit).addPropertyHistory(wlen);
  wit++;
  while(wit!=wit_end)
  {
    (**wit).addPropertyHistory(wlen);
    wit++;
  }
  std::vector<IndexType> samples_th(omp_get_max_threads(),0);
  myPeriod4WalkerDump=(Period4WalkerDump>0)?Period4WalkerDump:(nBlocks+1)*nSteps;
  samples_this_node = nTargetSamples/myComm->size();
  if (nTargetSamples%myComm->size() > myComm->rank())
    samples_this_node+=1;
  int samples_each_thread = samples_this_node/omp_get_max_threads();
  for (int ip=0; ip<omp_get_max_threads(); ++ip)
    samples_th[ip]=samples_each_thread;
  if(samples_this_node%omp_get_max_threads())
    for (int ip=0; ip < samples_this_node%omp_get_max_threads(); ++ip)
      samples_th[ip] +=1;
  int ndumps = std::max(samples_this_node/W.getActiveWalkers() + 1,2);
  myPeriod4WalkerDump = nBlocks*nSteps/ndumps;
  app_log() << "  Samples are dumped every " << myPeriod4WalkerDump << " steps " << endl;
  app_log() << "  Total Sample Size =" << nTargetSamples << endl;
  app_log() << "  Nodes Sample Size =" << samples_this_node << endl;
  for (int ip=0; ip<NumThreads; ++ip)
    app_log()  << "    Sample size for thread " <<ip<<" = " << samples_th[ip] << endl;
  #pragma omp critical
  for(int ip=0; ip<NumThreads; ++ip)
  {
    wClones[ip]->clearEnsemble();
    wClones[ip]->setNumSamples(samples_th[ip]);
  }
  t= Movers[0]->getTau();
  clearComponentMatrices();
  branchEngine->checkParameters(W);
  int mxage=mover_MaxAge;
  if(fixW)
  {
    if(BranchInterval<0)
      BranchInterval=nSteps;
    mxage=(mover_MaxAge<0)?0:mover_MaxAge;
    for(int ip=0; ip<Movers.size(); ++ip)
      Movers[ip]->MaxAge=mxage;
  }
  else
  {
    if(BranchInterval<0)
      BranchInterval=1;
    int miage=(QMCDriverMode[QMC_UPDATE_MODE])?1:5;
    mxage=(mover_MaxAge<0)?miage:mover_MaxAge;
    for(int ip=0; ip<Movers.size(); ++ip)
      Movers[ip]->MaxAge=mxage;
  }
  {
    ostringstream o;
    if(fixW)
      o << "  Fixed population using reconfiguration method\n";
    else
      o << "  Fluctuating population\n";
    o << "  Persisent walkers are killed after " << mxage << " MC sweeps\n";
    o << "  BranchInterval = " << BranchInterval << "\n";
    o << "  Steps per block = " << nSteps << "\n";
    o << "  Number of blocks = " << nBlocks << "\n";
    app_log() << o.str() << endl;
  }
  app_log() << "  DMC Engine Initialization = " << init_timer.elapsed() << " secs " << endl;
}

bool DMCOMPOPT::run()
{
  bool variablePop = (Reconfiguration == "no");
  resetUpdateEngines();
  //estimator does not need to collect data
  Estimators->setCollectionMode(true);
  Estimators->start(nBlocks);
  for(int ip=0; ip<NumThreads; ip++)
    Movers[ip]->startRun(nBlocks,false);
  Timer myclock;
  IndexType block = 0;
  IndexType updatePeriod=(QMCDriverMode[QMC_UPDATE_MODE])?Period4CheckProperties:(nBlocks+1)*nSteps;
  int nsampls(0);
  int g_nsampls(0);
  do // block
  {
    Estimators->startBlock(nSteps);
    for(int ip=0; ip<NumThreads; ip++)
      Movers[ip]->startBlock(nSteps);
    IndexType step = 0;
    for(IndexType step=0; step< nSteps; ++step, CurrentStep+=BranchInterval)
    {
      #pragma omp parallel
      {
        int ip=omp_get_thread_num();
        int now=CurrentStep;
        MCWalkerConfiguration::iterator wit(W.begin()+wPerNode[ip]), wit_first(W.begin()+wPerNode[ip]), wit2(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);

        for(int interval = 0; interval<BranchInterval-1; ++interval,++now)
        {
          Movers[ip]->advanceWalkers(wit,wit_end,false);
          while(wit2!=wit_end)
          {
            (**wit2).addPropertyHistoryPoint(Eindx,(**wit2).getPropertyBase()[LOCALENERGY]);
            wit2++;
          }
          wit2=wit_first;
          if (Period4WalkerDump&&((now+1)%myPeriod4WalkerDump==0))
          {
            wClones[ip]->saveEnsemble(wit2,wit_end);
            #pragma omp master
            nsampls+=W.getActiveWalkers();
          }
        }

        wClones[ip]->resetCollectables();

        Movers[ip]->advanceWalkers(wit,wit_end,false);
        wit2=wit_first;
        while(wit2!=wit_end)
        {
          (**wit2).addPropertyHistoryPoint(Eindx,(**wit2).getPropertyBase()[LOCALENERGY]);
          wit2++;
        }
        wit2=wit_first;
        if (myPeriod4WalkerDump&&((now+1)%myPeriod4WalkerDump==0))
        {
          wClones[ip]->saveEnsemble(wit2,wit_end);
          #pragma omp master
          nsampls+=W.getActiveWalkers();
        }

        Movers[ip]->setMultiplicity(wit,wit_end);

        if(QMCDriverMode[QMC_UPDATE_MODE] && now%updatePeriod == 0)
          Movers[ip]->updateWalkers(wit, wit_end);

      }//#pragma omp parallel
    }
//       branchEngine->debugFWconfig();
    #pragma omp parallel for
    for (int ip=0; ip<NumThreads; ++ip)
    {
      MCWalkerConfiguration::iterator wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
      while (wit!=wit_end)
      {
        Walker_t& thisWalker(**wit);
        Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
        wClones[ip]->loadWalker(thisWalker,true);
        psiClones[ip]->copyFromBuffer(*wClones[ip],w_buffer);
        vector<RealType> Dsaved(NumOptimizables,0);
        vector<RealType> HDsaved(NumOptimizables,0);
        psiClones[ip]->evaluateDerivatives(*wClones[ip],dummyOptVars[ip],Dsaved,HDsaved,true);//SH like deriv style
        #pragma omp critical
        fillComponentMatrices(Dsaved,HDsaved,thisWalker);
        wit++;
      }
    }
    branchEngine->branch(CurrentStep,W, branchClones);
    if(variablePop)
      FairDivideLow(W.getActiveWalkers(),NumThreads,wPerNode);
    Estimators->stopBlock(acceptRatio());
    block++;
    recordBlock(block);
    g_nsampls=nsampls;
    myComm->allreduce(g_nsampls);
  }
  while(((block<nBlocks) || (g_nsampls<nTargetSamples)) && myclock.elapsed()<MaxCPUSecs);
  //for(int ip=0; ip<NumThreads; ip++) Movers[ip]->stopRun();
  for(int ip=0; ip<NumThreads; ip++)
    *(RandomNumberControl::Children[ip])=*(Rng[ip]);
  //       adding weight index
  MCWalkerConfiguration::iterator wit(W.begin()), wit_end(W.end());
  while(wit!=wit_end)
  {
    (**wit).deletePropertyHistory();
    wit++;
  }
  Estimators->stop();
  return finalize(block);
}
}

/***************************************************************************
 * $RCSfile: DMCOMPOPT.cpp,v $   $Author: jnkim $
 * $Revision: 1620 $   $Date: 2007-01-14 18:12:23 -0600 (Sun, 14 Jan 2007) $
 * $Id: DMCOMPOPT.cpp 1620 2007-01-15 00:12:23Z jnkim $
 ***************************************************************************/
