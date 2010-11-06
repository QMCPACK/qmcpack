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
#include "QMCDrivers/WFMC/WFMCUpdateAll.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Message/OpenMP.h"
#include "Message/CommOperators.h"
//#define ENABLE_VMC_OMP_MASTER

namespace qmcplusplus
  {

  /// Constructor.
  VMCSingleOMP::VMCSingleOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
                             HamiltonianPool& hpool):
      QMCDriver(w,psi,h),  CloneManager(hpool),
      myWarmupSteps(0),UseDrift("yes"), logoffset(2.0), logepsilon(0)
  {
    RootName = "vmc";
    QMCType ="VMCSingleOMP";
    QMCDriverMode.set(QMC_UPDATE_MODE,1);
    QMCDriverMode.set(QMC_WARMUP,0);
    m_param.add(UseDrift,"useDrift","string"); m_param.add(UseDrift,"usedrift","string"); m_param.add(UseDrift,"use_drift","string");
    m_param.add(logepsilon,"logepsilon","double");
    m_param.add(logoffset,"logoffset","double");
    m_param.add(myRNWarmupSteps,"rnwarmupsteps","int");
    m_param.add(myWarmupSteps,"warmupSteps","int"); m_param.add(myWarmupSteps,"warmupsteps","int"); m_param.add(myWarmupSteps,"warmup_steps","int");
    m_param.add(nTargetSamples,"targetWalkers","int"); m_param.add(nTargetSamples,"targetwalkers","int"); m_param.add(nTargetSamples,"target_walkers","int");
  }

  bool VMCSingleOMP::run()
  {
    resetRun();

    //start the main estimator
    Estimators->start(nBlocks);

    for (int ip=0; ip<NumThreads; ++ip) Movers[ip]->startRun(nBlocks,false);

    for (int block=0;block<nBlocks; ++block)
      {
#pragma omp parallel for
        for (int ip=0; ip<NumThreads; ++ip)
          {
            //IndexType updatePeriod=(QMCDriverMode[QMC_UPDATE_MODE])?Period4CheckProperties:(nBlocks+1)*nSteps;
            IndexType updatePeriod=(QMCDriverMode[QMC_UPDATE_MODE])?Period4CheckProperties:0;
            //assign the iterators and resuse them
            MCWalkerConfiguration::iterator wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);

            Movers[ip]->startBlock(nSteps);
            int now_loc=CurrentStep;
            //rest the collectables and keep adding
            wClones[ip]->resetCollectables();
            for (int step=0; step<nSteps;++step)
              {
                Movers[ip]->advanceWalkers(wit,wit_end,false);
                Movers[ip]->accumulate(wit,wit_end);
                ++now_loc;
                if (updatePeriod&& now_loc%updatePeriod==0) Movers[ip]->updateWalkers(wit,wit_end);
                if (Period4WalkerDump&& now_loc%myPeriod4WalkerDump==0) wClones[ip]->saveEnsemble(wit,wit_end);                
                if(storeConfigs && (now_loc%storeConfigs == 0)) 
                {
                  ForwardWalkingHistory.storeConfigsForForwardWalking(*wClones[ip]);
                }
              }
            Movers[ip]->stopBlock(false);
          }//end-of-parallel for

        Estimators->accumulateCollectables(wClones,nSteps);

        CurrentStep+=nSteps;
        Estimators->stopBlock(estimatorClones);
        //why was this commented out? Are checkpoints stored some other way?
        if(storeConfigs) recordBlock(block);
      }//block

    Estimators->stop(estimatorClones);

    //copy back the random states
    for (int ip=0; ip<NumThreads; ++ip)
      *(RandomNumberControl::Children[ip])=*(Rng[ip]);

    //finalize a qmc section
    return finalize(nBlocks);
  }

  void VMCSingleOMP::resetRun()
  {
    makeClones(W,Psi,H);

    std::vector<IndexType> samples_th(omp_get_max_threads(),0);
    myPeriod4WalkerDump=(Period4WalkerDump>0)?Period4WalkerDump:(nBlocks+1)*nSteps;
    
    int samples_this_node = nTargetSamples/myComm->size();
    if (nTargetSamples%myComm->size() > myComm->rank()) samples_this_node+=1;
    
    int samples_each_thread = samples_this_node/omp_get_max_threads();
    for (int ip=0; ip<omp_get_max_threads(); ++ip) samples_th[ip]=samples_each_thread; 
    
    if(samples_this_node%omp_get_max_threads())
      for (int ip=0; ip < samples_this_node%omp_get_max_threads(); ++ip) samples_th[ip] +=1;
    
    app_log() << "  Samples are dumped every " << myPeriod4WalkerDump << " steps " << endl;
    app_log() << "  Total Sample Size =" << nTargetSamples << endl;  
    app_log() << "  Nodes Sample Size =" << samples_this_node << endl;  
    for (int ip=0; ip<NumThreads; ++ip)
      app_log()  << "    Sample size for thread " <<ip<<" = " << samples_th[ip] << endl;
    app_log() << "  Warmup Steps " << myWarmupSteps << endl;

    
    if (Movers.empty())
      {
        Movers.resize(NumThreads,0);
        branchClones.resize(NumThreads,0);
        estimatorClones.resize(NumThreads,0);
        Rng.resize(NumThreads,0);
        int nwtot=(W.getActiveWalkers()/NumThreads)*NumThreads;
        FairDivideLow(nwtot,NumThreads,wPerNode);

        app_log() << "  Initial partition of walkers ";
        std::copy(wPerNode.begin(),wPerNode.end(),ostream_iterator<int>(app_log()," "));
        app_log() << endl;

#pragma omp parallel for
        for (int ip=0; ip<NumThreads; ++ip)
        {
          ostringstream os;
          estimatorClones[ip]= new EstimatorManager(*Estimators);//,*hClones[ip]);
          estimatorClones[ip]->resetTargetParticleSet(*wClones[ip]);
          estimatorClones[ip]->setCollectionMode(false);
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
            if (UseDrift == "rn")
            {
              os <<"  PbyP moves with RN, using VMCUpdatePbyPSampleRN"<<endl;
              Movers[ip]=new VMCUpdatePbyPSampleRN(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
              Movers[ip]->setLogEpsilon(logepsilon);
              // Movers[ip]=new VMCUpdatePbyPWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
            }
            else if (UseDrift == "yes")
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
            if (UseDrift == "rn")
            {
              os <<"  walker moves with RN, using VMCUpdateAllSampleRN"<<endl;
              Movers[ip]=new VMCUpdateAllSampleRN(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
              Movers[ip]->setLogEpsilon(logepsilon);
            }
            else if (UseDrift == "yes")
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

          if(ip==0) app_log() << os.str() << endl;
        }
      }
        
#pragma omp parallel
      {
        int ip=omp_get_thread_num();
        Movers[ip]->put(qmcNode);
        Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);

        if (QMCDriverMode[QMC_UPDATE_MODE])
          Movers[ip]->initWalkersForPbyP(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
        else
          Movers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);


        for (int prestep=0; prestep<myWarmupSteps; ++prestep)
          Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);

        if (myWarmupSteps && QMCDriverMode[QMC_UPDATE_MODE])
          Movers[ip]->updateWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
        
#pragma omp critical
        {
          wClones[ip]->clearEnsemble();
          wClones[ip]->setNumSamples(samples_th[ip]);
        }
      }
    
    
    if ((UseDrift == "rn")&&(logepsilon==0.0))
    {
      long double psi2_sum(0.0);
      long double psi4_sum(0.0);
      RealType psi2_0_0(0.0);
      RealType nw_local(0.0);
      
      if (myRNWarmupSteps==0) myRNWarmupSteps=myWarmupSteps;
      #pragma omp parallel
      {
        int ip=omp_get_thread_num();
        Movers[ip]->setLogEpsilon(logepsilon);
        for (int prestep=0; prestep<myRNWarmupSteps; ++prestep)
        {
          Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
          #pragma omp flush
          if((ip==0)&&(prestep>5))
          {
            MCWalkerConfiguration::iterator wit(W.begin()), wit_end(W.end());
            if (prestep==0)
            {
              if (myComm->rank()==0) psi2_0_0=2.0* (*wit)->getPropertyBase()[LOGPSI];
              myComm->bcast(psi2_0_0);
            }
            nw_local += wit_end-wit;
            while (wit!=wit_end)
            {
              long double t(expl(2.0*(*wit)->getPropertyBase()[LOGPSI]-psi2_0_0));
              psi2_sum += t;
              psi4_sum += t*t;
              wit++;
            }
//             app_log()<<logl(psi2_sum/nw_local)<<"  "<<logl(sqrtl((psi4_sum/nw_local - (psi2_sum/nw_local * psi2_sum/nw_local))/(1+prestep)) )<<endl;
          }
          #pragma omp barrier
          
        }
      }

//           we need to set the logepsilon to something reasonable
        myComm->allreduce(psi2_sum);
        myComm->allreduce(psi4_sum);
        myComm->allreduce(nw_local);

        long double p2(psi2_sum*nw_local);
        RealType eps0 = logl(p2) + psi2_0_0 - logoffset;
        app_log()<<" Using logepsilon "<<eps0<<endl;
        app_log()<<" log<Psi^2>= "<<logl(p2)+ psi2_0_0<<endl;
        long double sigp2(psi4_sum*nw_local-p2*p2);
        app_log()<<" logError  = "<< logl(sqrtl(sigp2*nw_local))<<endl;
        for (int ip=0; ip<NumThreads; ++ip) Movers[ip]->setLogEpsilon(eps0);
        #pragma omp parallel
        {
          int ip=omp_get_thread_num();
          for (int prestep=0; prestep<myRNWarmupSteps; ++prestep)
          {
            Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
          }
        }
        
//         bool badEPS(true);
//         RealType stabilizer=0;
//         RealType eps0 = logl(psi2_sum*nw_local) + psi2_0_0;
// 
//         int Accept(0), Reject(0);
//           for (int ip=0; ip<NumThreads; ++ip) Accept+=Movers[ip]->nAccept;
//           for (int ip=0; ip<NumThreads; ++ip) Reject+=Movers[ip]->nReject;
//           myComm->allreduce(Accept);
//           myComm->allreduce(Reject);
//           RealType accpt=(1.0*Accept)/(Accept+Reject);
//         app_log()<<"Pre: "<<accpt<<endl;
//         while (badEPS)
//         {
//           eps0 = logl(psi2_sum*nw_local) + psi2_0_0 + stabilizer;
//             #pragma omp parallel
//             {
//               int ip=omp_get_thread_num();
//               Movers[ip]->setLogEpsilon(eps0);
//               for (int prestep=0; prestep<myWarmupSteps; ++prestep)
//               {
//                 Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
//                 if(ip==0)
//                 {
//                   MCWalkerConfiguration::iterator wit(W.begin()), wit_end(W.end());
//                   nw_local += wit_end-wit;
//                   while (wit!=wit_end)
//                   {
//                     psi2_sum += expl(2.0*(*wit)->getPropertyBase()[LOGPSI]-psi2_0_0);
//                     wit++;
//                   }
//                 }
//                 #pragma omp barrier
//               }
//             }
//           myComm->allreduce(psi2_sum);
//           myComm->allreduce(nw_local);
//         
//           for (int ip=0; ip<NumThreads; ++ip) 
//             {
//               Movers[ip]->nAccept=0;
//               Movers[ip]->nReject=0;
//               Accept=0;
//               Reject=0;
//             }
//             #pragma omp parallel
//             {
//               int ip=omp_get_thread_num();
//               for (int prestep=0; prestep<myRNWarmupSteps; ++prestep)
//                 Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
//             }
//             for (int ip=0; ip<NumThreads; ++ip) Accept+=Movers[ip]->nAccept;
//             for (int ip=0; ip<NumThreads; ++ip) Reject+=Movers[ip]->nReject;
//             myComm->allreduce(Accept);
//             myComm->allreduce(Reject);
//             RealType accptE=(1.0*Accept)/(Accept+Reject);
//             app_log()<<Accept<<"  "<<Reject<<endl;
//             app_log()<<stabilizer<<"  "<<accptE<<endl;
//             if (accptE<rn_accept_target) badEPS=false;
//             else if (accptE>99) stabilizer -= 8.0;
//             else if (accptE>rn_accept_target+0.3) stabilizer -=1.0;
//             else if (accptE>rn_accept_target+0.2) stabilizer -=0.1;
//             else stabilizer -= 0.01;
//         }
//         app_log()<<"  Using LogEpsilon of: "<<eps0<<endl;
    }
    else if (UseDrift == "rn")
    {
      app_log()<<"  Using LogEpsilon of: "<<logepsilon<<endl;
#pragma omp parallel
      {
        int ip=omp_get_thread_num();
        Movers[ip]->setLogEpsilon(logepsilon);
        for (int prestep=0; prestep<myRNWarmupSteps; ++prestep)
        {
          Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
        }
      }
    }
    #pragma omp parallel
    {
      int ip=omp_get_thread_num();
      if (myWarmupSteps && QMCDriverMode[QMC_UPDATE_MODE])
        Movers[ip]->updateWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
    }
            
    myWarmupSteps=0;
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
    //nothing to add
    return true;
  }
}

/***************************************************************************
 * $RCSfile: VMCSingleOMP.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: VMCSingleOMP.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $
 ***************************************************************************/
