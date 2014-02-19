#include "QMCDrivers/RMC/RMCSingleOMP.h"
#include "QMCDrivers/RMC/RMCUpdatePbyP.h"
#include "QMCDrivers/RMC/RMCUpdateAll.h"
#include "QMCDrivers/DriftOperators.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Message/OpenMP.h"
#include "Message/CommOperators.h"
#include "tau/profiler.h"
#include "Particle/Reptile.h"


namespace qmcplusplus
{

/// Constructor.
RMCSingleOMP::RMCSingleOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
                           HamiltonianPool& hpool, WaveFunctionPool& ppool):
  QMCDriver(w,psi,h,ppool),  CloneManager(hpool), rescaleDrift("no"), beta(-1), beads(-1),resizeReptile(0)
{
  RootName = "rmc";
  QMCType ="RMCSingleOMP";
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  QMCDriverMode.set(QMC_WARMUP,0);
  m_param.add(rescaleDrift,"drift","string");
  m_param.add(beta,"beta","double");
  m_param.add(beads,"beads","int");
  m_param.add(resizeReptile,"resize","int");
  
 // if (w.reptile_new==0){w.reptile_new = new Reptile_new(&w);};
  
 // w.reptile_new->loadFromWalkerList(w.WalkerList);
//  app_log()<<"  reptile_new.size()="<<w.reptile_new->size()<<endl;
//  app_log()<<"  WalkerList.size()="<<w.WalkerList.size()<<endl;
  
  Action.resize(3);
  Action[0]=w.addProperty("ActionBackward");
  Action[1]=w.addProperty("ActionForward");
  Action[2]=w.addProperty("ActionLocal");
//     app_log()<<Action[0]<<" "<<Action[1]<<" "<<Action[2]<<endl;
  TransProb.resize(2);
  TransProb[0]=w.addProperty("TransProbBackward");
  TransProb[1]=w.addProperty("TransProbForward");
//     app_log()<<TransProb[0]<<" "<<TransProb[1]<<endl;
}

bool RMCSingleOMP::run()
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
      #pragma omp parallel
      {
        int ip=omp_get_thread_num();
        //IndexType updatePeriod=(QMCDriverMode[QMC_UPDATE_MODE])?Period4CheckProperties:(nBlocks+1)*nSteps;
        IndexType updatePeriod=(QMCDriverMode[QMC_UPDATE_MODE])?Period4CheckProperties:0;
        //assign the iterators and resuse them
        MCWalkerConfiguration::iterator wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
//         MCWalkerConfiguration::iterator wit(first);
        Movers[ip]->startBlock(nSteps);
        int now_loc=CurrentStep;

        RealType cnorm=1.0/static_cast<RealType>(wPerNode[ip+1]-wPerNode[ip]);



        for (int step=0; step<nSteps; ++step)
          {
            //collectables are reset, it is accumulated while advancing walkers
            wClones[ip]->resetCollectables();
//           wit=first;
            Movers[ip]->advanceWalkers(wit,wit_end,false);
            if(has_collectables)
              wClones[ip]->Collectables *= cnorm;
//           wit=first;
            Movers[ip]->accumulate(wit,wit_end);
            
            ++now_loc;
            //if (updatePeriod&& now_loc%updatePeriod==0) Movers[ip]->updateWalkers(wit,wit_end);
            if (Period4WalkerDump&& now_loc%myPeriod4WalkerDump==0)
              wClones[ip]->saveEnsemble(wit,wit_end);
              
             if (ip==0) branchEngine->collect(CurrentStep, W, branchClones);  //Ray Clay:  For now, collects and syncs based on first reptile.  Need a better way to do this.
          }
       // wClones[ip]->reptile->calcTauScaling();
       // wClones[ip]->reptile->calcERun();
		//branchEngine->collect(CurrentStep, W, branchClones);
       // wClones[ip]->reptile->resetR2Avg();
        Movers[ip]->stopBlock(false);
        //  if (block > 2*wClones[ip]->reptile->nbeads){
        //	 wClones[ip]->reptile->printState();
        //	 ((RMCUpdateAllWithDrift*)Movers[ip])->checkReptile(wit,wit_end);

        //	if (block>1) branchEngine->updateReptileStats(*wClones[ip]);

        //}
      }//end-of-parallel for
      //Estimators->accumulateCollectables(wClones,nSteps);
      CurrentStep+=nSteps;
     // branchEngine->collect(CurrentStep, W, branchClones);
      Estimators->stopBlock(estimatorClones);
      //why was this commented out? Are checkpoints stored some other way?
      if(storeConfigs)
        recordBlock(block);
    }//block
  Estimators->stop(estimatorClones);
  //copy back the random states
  for (int ip=0; ip<NumThreads; ++ip)
    *(RandomNumberControl::Children[ip])=*(Rng[ip]);
  //finalize a qmc section
  return finalize(nBlocks);
}

void RMCSingleOMP::resetRun()
{
  if(resizeReptile<2)
    {
//     compute number of beads you need. either by looking at the beads parameter, or by looking at beta and tau.
      if(beads<1)
        beads=beta/Tau;
      else
        beta=beads*Tau;
      //here we resize the number of walkers if it is wrong.
      //new walkers go on the end of W
      int nwtot= beads*NumThreads;
      FairDivideLow(nwtot,NumThreads,wPerNode);
      if(W.getActiveWalkers()-nwtot !=0)
        addWalkers(nwtot-W.getActiveWalkers());
        
       //   W.reptile_new->loadFromWalkerList(W.WalkerList);
 // app_log()<<"  reptile_new.size()="<<W.reptile_new->size()<<endl;
  //app_log()<<"  WalkerList.size()="<<W.WalkerList.size()<<endl;
    }
  else
    {
      ///Norm
      app_log()<<"resizing the reptile not yet implemented."<<endl;
    }
  makeClones(W,Psi,H);
  myPeriod4WalkerDump=(Period4WalkerDump>0)?Period4WalkerDump:(nBlocks+1)*nSteps;
  if (Movers.empty())
    {
      Movers.resize(NumThreads,0);
      branchClones.resize(NumThreads,0);
      estimatorClones.resize(NumThreads,0);
      traceClones.resize(NumThreads,0);
      Rng.resize(NumThreads,0);
      //   W.loadEnsemble(wClones);
      branchEngine->initReptile(W);
#if !defined(BGP_BUG)
      #pragma omp parallel for
#endif
      for(int ip=0; ip<NumThreads; ++ip)
        {
          ostringstream os;
          estimatorClones[ip]= new EstimatorManager(*Estimators);//,*hClones[ip]);
          estimatorClones[ip]->resetTargetParticleSet(*wClones[ip]);
          estimatorClones[ip]->setCollectionMode(false);
          Rng[ip]=new RandomGenerator_t(*(RandomNumberControl::Children[ip]));
          traceClones[ip] = Traces->makeClone();
          hClones[ip]->setRandomGenerator(Rng[ip]);
          branchClones[ip] = new BranchEngineType(*branchEngine);
            // branchClones[ip]->initReptile(W);
          if (QMCDriverMode[QMC_UPDATE_MODE])
            {
              os <<"  PbyP moves with drift, using RMCUpdatePbyPWithDriftFast"<<endl;
              Movers[ip]=new RMCUpdatePbyPWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip],Action,TransProb);

            }
          else
            {
              os <<"  walker moves with drift, using RMCUpdateAllWithDriftFast"<<endl;
              Movers[ip]=new RMCUpdateAllWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip],Action,TransProb);
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
  app_log().flush();
#if !defined(BGP_BUG)
  #pragma omp parallel for
#endif
  for(int ip=0; ip<NumThreads; ++ip)
    { 
      Movers[ip]->put(qmcNode);
      Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip], traceClones[ip]);
      wClones[ip]->reptile = new Reptile(*wClones[ip], W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
      //app_log()<<"Thread # "<<ip<<endl;
     // printf(" Thread# %d  WalkerList.size()=%d \n",ip,wClones[ip]->WalkerList.size());
	  
      // wClones[ip]->reptile->printState();
      wClones[ip]->activeBead= 0;
      wClones[ip]->direction = +1;

      if (QMCDriverMode[QMC_UPDATE_MODE])
        {
          Movers[ip]->initWalkersForPbyP(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
          Movers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
        }
      else
        {
          Movers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
          //  wClones[ip]->reptile = new Reptile(*wClones[ip], W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
        }

      //set up initial action and transprob.
      MCWalkerConfiguration::iterator wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);

      //IndexType initsteps = wClones[ip]->reptile->nbeads * 2;


/// thermalization Norm
//         for (int prestep=0; prestep<nWarmupSteps; ++prestep)
//           Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
//
//         if (nWarmupSteps && QMCDriverMode[QMC_UPDATE_MODE])
//           Movers[ip]->updateWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
    }
  branchEngine->checkParameters(W);
}

bool
RMCSingleOMP::put(xmlNodePtr q)
{
  //nothing to add
  return true;
}
}
