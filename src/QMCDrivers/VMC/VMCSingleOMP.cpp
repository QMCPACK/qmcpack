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
#include "QMCDrivers/VMC/VMCSingleOMP.h"
#include "QMCDrivers/VMC/VMCUpdatePbyP.h"
#include "QMCDrivers/VMC/VMCUpdateAll.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Message/OpenMP.h"
#include "Message/CommOperators.h"
#include "tau/profiler.h"
#include <qmc_common.h>
//#define ENABLE_VMC_OMP_MASTER
#include "ADIOS/ADIOS_profile.h"

namespace qmcplusplus
{

/// Constructor.
VMCSingleOMP::VMCSingleOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
                           HamiltonianPool& hpool, WaveFunctionPool& ppool):
  QMCDriver(w,psi,h,ppool),  CloneManager(hpool),
  UseDrift("yes") //, logoffset(2.0), logepsilon(0)
{
  RootName = "vmc";
  QMCType ="VMCSingleOMP";
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  QMCDriverMode.set(QMC_WARMUP,0);
  m_param.add(UseDrift,"useDrift","string");
  m_param.add(UseDrift,"usedrift","string");
  m_param.add(UseDrift,"use_drift","string");

  prevSteps=nSteps;
  prevStepsBetweenSamples=nStepsBetweenSamples;
}

bool VMCSingleOMP::run()
{
  resetRun();
  //start the main estimator
  Estimators->start(nBlocks);
  for (int ip=0; ip<NumThreads; ++ip)
    Movers[ip]->startRun(nBlocks,false);
  Traces->startRun(nBlocks,traceClones);
  const bool has_collectables=W.Collectables.size();
  for (int block=0; block<nBlocks; ++block)
  {
    double comp_start = MPI_Wtime(); 
    #pragma omp parallel
    {
      int ip=omp_get_thread_num();
      //IndexType updatePeriod=(QMCDriverMode[QMC_UPDATE_MODE])?Period4CheckProperties:(nBlocks+1)*nSteps;
      IndexType updatePeriod=(QMCDriverMode[QMC_UPDATE_MODE])?Period4CheckProperties:0;
      //assign the iterators and resuse them
      MCWalkerConfiguration::iterator wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
      Movers[ip]->startBlock(nSteps);
      int now_loc=CurrentStep;
      RealType cnorm=1.0/static_cast<RealType>(wPerNode[ip+1]-wPerNode[ip]);
      for (int step=0; step<nSteps; ++step)
      {
        Movers[ip]->set_step(now_loc);
        //collectables are reset, it is accumulated while advancing walkers
        wClones[ip]->resetCollectables();
        Movers[ip]->advanceWalkers(wit,wit_end,false);
        if(has_collectables)
          wClones[ip]->Collectables *= cnorm;
        Movers[ip]->accumulate(wit,wit_end);
        ++now_loc;
        //if (updatePeriod&& now_loc%updatePeriod==0) Movers[ip]->updateWalkers(wit,wit_end);
        if (Period4WalkerDump&& now_loc%Period4WalkerDump==0)
          wClones[ip]->saveEnsemble(wit,wit_end);
//           if(storeConfigs && (now_loc%storeConfigs == 0))
//             ForwardWalkingHistory.storeConfigsForForwardWalking(*wClones[ip]);
      }
      Movers[ip]->stopBlock(false);
    }//end-of-parallel for
    //Estimators->accumulateCollectables(wClones,nSteps);
    CurrentStep+=nSteps;
    Estimators->stopBlock(estimatorClones);
    app_log()<<block<<" comp time "<<MPI_Wtime() - comp_start<<endl; 
    Traces->write_buffers(traceClones, block);
    if(storeConfigs)
      recordBlock(block);
  }//block
  Estimators->stop(estimatorClones);
  for (int ip=0; ip<NumThreads; ++ip)
    Movers[ip]->stopRun2();
  Traces->stopRun();
  //copy back the random states
  for (int ip=0; ip<NumThreads; ++ip)
    *(RandomNumberControl::Children[ip])=*(Rng[ip]);
  ///write samples to a file
  bool wrotesamples=DumpConfig;
  if(DumpConfig)
  {
    wrotesamples=W.dumpEnsemble(wClones,wOut,myComm->size(),nBlocks);
    if(wrotesamples)
      app_log() << "  samples are written to the config.h5" << endl;
  }
  //finalize a qmc section
  return finalize(nBlocks,!wrotesamples);
}

void VMCSingleOMP::resetRun()
{
  ////only VMC can overwrite this
  if(nTargetPopulation>0)
    branchEngine->iParam[SimpleFixedNodeBranch::B_TARGETWALKERS]=static_cast<int>(std::ceil(nTargetPopulation));
  makeClones(W,Psi,H);
  FairDivideLow(W.getActiveWalkers(),NumThreads,wPerNode);
  app_log() << "  Initial partition of walkers ";
  std::copy(wPerNode.begin(),wPerNode.end(),ostream_iterator<int>(app_log()," "));
  app_log() << endl;

  bool movers_created=false;
  if (Movers.empty())
  {
    movers_created=true;
    Movers.resize(NumThreads,0);
    branchClones.resize(NumThreads,0);
    estimatorClones.resize(NumThreads,0);
    traceClones.resize(NumThreads,0);
    Rng.resize(NumThreads,0);
#if !defined(BGP_BUG)
    #pragma omp parallel for
#endif
    for(int ip=0; ip<NumThreads; ++ip)
    {
      ostringstream os;
      estimatorClones[ip]= new EstimatorManager(*Estimators);//,*hClones[ip]);
      estimatorClones[ip]->resetTargetParticleSet(*wClones[ip]);
      estimatorClones[ip]->setCollectionMode(false);
      traceClones[ip] = Traces->makeClone();
      Rng[ip]=new RandomGenerator_t(*(RandomNumberControl::Children[ip]));
      hClones[ip]->setRandomGenerator(Rng[ip]);
      branchClones[ip] = new BranchEngineType(*branchEngine);
      //         if(reweight=="yes")
      //         {
      //           if (ip== 0) app_log() << "  WFMCUpdateAllWithReweight"<<endl;
      //           Movers[ip]=new WFMCUpdateAllWithReweight(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip],weightLength,Eindex);
      //         }
      //         else
      //           if (reweight=="psi")
      //           {
      //             os << "  Sampling Psi to increase number of walkers near nodes"<<endl;
      //             if (QMCDriverMode[QMC_UPDATE_MODE]) Movers[ip]=new VMCUpdatePbyPSamplePsi(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
      //             else Movers[ip]=new VMCUpdateAllSamplePsi(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
      //           }
      //           else
      if (QMCDriverMode[QMC_UPDATE_MODE])
      {
        //             if (UseDrift == "rn")
        //             {
        //               os <<"  PbyP moves with RN, using VMCUpdatePbyPSampleRN"<<endl;
        //               Movers[ip]=new VMCUpdatePbyPSampleRN(*wClones[ip],*psiClones[ip],*guideClones[ip],*hClones[ip],*Rng[ip]);
        //               Movers[ip]->setLogEpsilon(logepsilon);
        //               // Movers[ip]=new VMCUpdatePbyPWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
        //             }
        //             else
        if (UseDrift == "yes")
        {
          os <<"  PbyP moves with drift, using VMCUpdatePbyPWithDriftFast"<<endl;
          Movers[ip]=new VMCUpdatePbyPWithDriftFast(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
          // Movers[ip]=new VMCUpdatePbyPWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
        }
        else
        {
          os <<"  PbyP moves with |psi^2|, using VMCUpdatePbyP"<<endl;
          Movers[ip]=new VMCUpdatePbyP(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
        }
        //Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
      }
      else
      {
        //             if (UseDrift == "rn")
        //             {
        //               os <<"  walker moves with RN, using VMCUpdateAllSampleRN"<<endl;
        //               Movers[ip]=new VMCUpdateAllSampleRN(*wClones[ip],*psiClones[ip],*guideClones[ip],*hClones[ip],*Rng[ip]);
        //               Movers[ip]->setLogEpsilon(logepsilon);
        //             }
        //             else
        if (UseDrift == "yes")
        {
          os <<"  walker moves with drift, using VMCUpdateAllWithDriftFast"<<endl;
          Movers[ip]=new VMCUpdateAllWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
        }
        else
        {
          os <<"  walker moves with |psi|^2, using VMCUpdateAll"<<endl;
          Movers[ip]=new VMCUpdateAll(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
        }
        //Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
      }
      Movers[ip]->nSubSteps=nSubSteps;
      if(ip==0)
        app_log() << os.str() << endl;
    }
  }
  else
  {
#if !defined(BGP_BUG)
    #pragma omp parallel for
#endif
    for(int ip=0; ip<NumThreads; ++ip)
    {
      traceClones[ip]->transfer_state_from(*Traces);
    }
  }
  app_log() << "  Total Sample Size   =" << nTargetSamples << endl;
  app_log() << "  Walker distribution on root = ";
  std::copy(wPerNode.begin(),wPerNode.end(),ostream_iterator<int>(app_log()," "));
  app_log() << endl;
  //app_log() << "  Sample Size per node=" << samples_this_node << endl;
  //for (int ip=0; ip<NumThreads; ++ip)
  //  app_log()  << "    Sample size for thread " <<ip<<" = " << samples_th[ip] << endl;
  app_log().flush();
#if !defined(BGP_BUG)
  #pragma omp parallel for
#endif
  for(int ip=0; ip<NumThreads; ++ip)
  {
    //int ip=omp_get_thread_num();
    Movers[ip]->put(qmcNode);
    Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip],traceClones[ip]);
    if (QMCDriverMode[QMC_UPDATE_MODE])
      Movers[ip]->initWalkersForPbyP(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
    else
      Movers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
//       if (UseDrift != "rn")
//       {
    for (int prestep=0; prestep<nWarmupSteps; ++prestep)
      Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
    if (nWarmupSteps && QMCDriverMode[QMC_UPDATE_MODE])
      Movers[ip]->updateWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
//       }
  }

  if(movers_created)
  {
    size_t before=qmc_common.memory_allocated;
    app_log() << "  Anonymous Buffer size per walker " << W[0]->DataSet.size() << endl;
    qmc_common.memory_allocated+=W.getActiveWalkers()*W[0]->DataSet.size()*sizeof(OHMMS_PRECISION);
    qmc_common.print_memory_change("VMCSingleOMP::resetRun",before);
  }
  //     //JNKIM: THIS IS BAD AND WRONG
//     if (UseDrift == "rn")
//     {
//       RealType avg_w(0);
//       RealType n_w(0);
// #pragma omp parallel
//       {
//         int ip=omp_get_thread_num();
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
//           }
//
//         for (int prestep=0; prestep<nWarmupSteps; ++prestep)
//           Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
//
//         if (nWarmupSteps && QMCDriverMode[QMC_UPDATE_MODE])
//           Movers[ip]->updateWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
//       }
//     }
  for(int ip=0; ip<NumThreads; ++ip)
    wClones[ip]->clearEnsemble();
  if(nSamplesPerThread)
    for(int ip=0; ip<NumThreads; ++ip)
      wClones[ip]->setNumSamples(nSamplesPerThread);
  nWarmupSteps=0;
  //Used to debug and benchmark opnemp
  //#pragma omp parallel for
  //    for(int ip=0; ip<NumThreads; ip++)
  //    {
  //      Movers[ip]->benchMark(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],ip);
  //    }
}

bool
VMCSingleOMP::put(xmlNodePtr q)
{
  //grep minimumTargetWalker
  int target_min=-1;
  ParameterSet p;
  p.add(target_min,"minimumtargetwalkers","int"); //p.add(target_min,"minimumTargetWalkers","int"); 
  p.add(target_min,"minimumsamples","int"); //p.add(target_min,"minimumSamples","int");
  p.put(q);

  app_log() << "\n<vmc function=\"put\">"
    << "\n  qmc_counter=" << qmc_common.qmc_counter << "  my_counter=" << MyCounter<< endl;
  if(qmc_common.qmc_counter && MyCounter)
  {
    nSteps=prevSteps;
    nStepsBetweenSamples=prevStepsBetweenSamples;
  }
  else
  {
    int nw=W.getActiveWalkers();
    //compute samples and overwrite steps for the given samples
    int Nthreads = omp_get_max_threads();
    int Nprocs=myComm->size();
    //target samples set by samples or samplesperthread/dmcwalkersperthread
    nTargetPopulation=std::max(nTargetPopulation,nSamplesPerThread*Nprocs*Nthreads);
    nTargetSamples=static_cast<int>(std::ceil(nTargetPopulation));

    if(nTargetSamples)
    {
      int nwtot=nw*Nprocs;  //total number of walkers used by this qmcsection
      nTargetSamples=std::max(nwtot,nTargetSamples);
      if(target_min>0) 
      { 
        nTargetSamples=std::max(nTargetSamples,target_min);
        nTargetPopulation=std::max(nTargetPopulation,static_cast<double>(target_min));
      }
      nTargetSamples=((nTargetSamples+nwtot-1)/nwtot)*nwtot; // nTargetSamples are always multiples of total number of walkers
      nSamplesPerThread=nTargetSamples/Nprocs/Nthreads;
      int ns_target=nTargetSamples*nStepsBetweenSamples; //total samples to generate
      int ns_per_step=Nprocs*nw;  //total samples per step
      nSteps=std::max(nSteps,(ns_target/ns_per_step+nBlocks-1)/nBlocks);
      Period4WalkerDump=nStepsBetweenSamples=(ns_per_step*nSteps*nBlocks)/nTargetSamples;
    }
    else
    {
      Period4WalkerDump = nStepsBetweenSamples=(nBlocks+1)*nSteps; //some positive number, not used
      nSamplesPerThread=0;
    }
  }
  prevSteps=nSteps;
  prevStepsBetweenSamples=nStepsBetweenSamples;

  app_log() << "  time step      = " << Tau << endl;
  app_log() << "  blocks         = " << nBlocks << endl;
  app_log() << "  steps          = " << nSteps << endl;
  app_log() << "  substeps       = " << nSubSteps << endl;
  app_log() << "  current        = " << CurrentStep << endl;
  app_log() << "  target samples = " << nTargetPopulation << endl;
  app_log() << "  walkers/mpi    = " << W.getActiveWalkers() << endl << endl;
  app_log() << "  stepsbetweensamples = " << nStepsBetweenSamples << endl;

  m_param.get(app_log());

  if(DumpConfig)
  {
    app_log() << "  DumpConfig==true Configurations are dumped to config.h5 with a period of " << Period4CheckPoint << " blocks" << endl;
  }
  else
  {
    app_log() << "  DumpConfig==false Nothing (configurations, state) will be saved." << endl;
  }

  if (Period4WalkerDump>0)
    app_log() << "  Walker Samples are dumped every " << Period4WalkerDump << " steps." << endl;

  app_log() << "</vmc>" << endl;
  app_log().flush();

  return true;
}
}

/***************************************************************************
 * $RCSfile: VMCSingleOMP.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: VMCSingleOMP.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $
 ***************************************************************************/
