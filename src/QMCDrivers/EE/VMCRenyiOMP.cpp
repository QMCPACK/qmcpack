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
#include "QMCDrivers/EE/VMCRenyiOMP.h"
#include "QMCDrivers/EE/VMCMultiRenyiOMP.h"
#include "QMCDrivers/EE/RenyiUpdatePbyP.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Message/OpenMP.h"
#include "Message/CommOperators.h"
#include "tau/profiler.h"
//#define ENABLE_VMC_OMP_MASTER

namespace qmcplusplus
  {

  /// Constructor.
  VMCRenyiOMP::VMCRenyiOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
                             HamiltonianPool& hpool, WaveFunctionPool& ppool):
      QMCDriver(w,psi,h,ppool),  CloneManager(hpool),
      myWarmupSteps(0),UseDrift("yes"), computeEE("sphere"),vsize(1),EEN(2),Nmax(20),Nmin(0)
  {
    RootName = "vmc";
    QMCType ="VMCRenyiOMP";
    QMCDriverMode.set(QMC_UPDATE_MODE,1);
    QMCDriverMode.set(QMC_WARMUP,0);
    m_param.add(UseDrift,"useDrift","string"); m_param.add(UseDrift,"usedrift","string"); m_param.add(UseDrift,"use_drift","string");
    m_param.add(computeEE,"computeEE","string");
    m_param.add(vsize,"EEsize","double");
    m_param.add(EEN,"EEN","int");
    m_param.add(Nmin,"Nmin","int");
    m_param.add(Nmax,"Nmax","int");
    m_param.add(myWarmupSteps,"warmupSteps","int"); m_param.add(myWarmupSteps,"warmupsteps","int"); m_param.add(myWarmupSteps,"warmup_steps","int");
    m_param.add(nTargetSamples,"targetWalkers","int"); m_param.add(nTargetSamples,"targetwalkers","int"); m_param.add(nTargetSamples,"target_walkers","int");
  }

  bool VMCRenyiOMP::run()
  {
    resetRun();
    if(myComm->rank()==0)
      {   
        ee_dat.str("");
        ee_dat<<RootName.c_str()<<".ee.dat";

        file_out.open(ee_dat.str().c_str(),fstream::out | fstream::trunc);
        file_out<<"# indx";
        file_out<<"  EE_S";
        for(int i(Nmin);i<Nmax+1;i++)
          file_out<<"  N"<<i;
        file_out<<endl;
        file_out.close();
      }

    //start the main estimator
    Estimators->start(nBlocks);
    for (int ip=0; ip<NumThreads; ++ip) RenyiMovers[ip]->startRun(nBlocks,false);

    const bool has_collectables=W.Collectables.size();

    hpmStart(QMC_VMC_0_EVENT,"vmc::main");
    for (int block=0;block<nBlocks; ++block)
    {
#pragma omp parallel
      {
        int ip=omp_get_thread_num();
        MCWalkerConfiguration::iterator wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
        RenyiMovers[ip]->startBlock(nSteps);
        for (int step=0; step<nSteps;++step)
          RenyiMovers[ip]->advanceWalkers(wit,wit_end,false);
        wit=W.begin()+wPerNode[ip];
        RenyiMovers[ip]->accumulate(wit,wit_end);
        RenyiMovers[ip]->stopBlock(false);
      }
#pragma omp barrier
      RealType avgsgn(0);
      std::vector<RealType> n_stats(Nmax+1-Nmin,0);
      for (int ip=0; ip<NumThreads; ++ip)
        avgsgn+=RenyiMovers[ip]->get_stats(n_stats);
      myComm->allreduce(avgsgn);
      myComm->allreduce(n_stats);
      RealType nrm=1.0/(NumThreads*myComm->size());
      
      if(myComm->rank()==0)
      {
        file_out.open(ee_dat.str().c_str(),fstream::out | fstream::app);
        file_out<<block<<" ";
        file_out<<avgsgn*nrm<<" ";
        for(int i(0);i<(Nmax+1-Nmin);i++)
          file_out<<n_stats[i]*nrm<<" ";
        file_out<<endl;
        file_out.close();
      }
      
//       for (int ip=0; ip<NumThreads; ++ip)
//         if(RenyiMovers[ip]->regions[W.getTotalNum()+1]==0)
//           RenyiMovers[ip]->print_all();
        
      for (int ip=0; ip<NumThreads; ++ip)
        RenyiMovers[ip]->clear_stats();
      Estimators->stopBlock(estimatorClones);
    }
    hpmStop(QMC_VMC_0_EVENT);

    Estimators->stop(estimatorClones);
    
    

    //copy back the random states
    for (int ip=0; ip<NumThreads; ++ip)
      *(RandomNumberControl::Children[ip])=*(Rng[ip]);

    //finalize a qmc section
    return finalize(nBlocks);
  }

  void VMCRenyiOMP::resetRun()
  {
    makeClones(W,Psi,H);

    std::vector<IndexType> samples_th(omp_get_max_threads(),0);
    app_log() << "  Warmup Steps " << myWarmupSteps << endl;
    
    if (RenyiMovers.empty())
      {
        RenyiMovers.resize(NumThreads,0);
        branchClones.resize(NumThreads,0);
        estimatorClones.resize(NumThreads,0);
        Rng.resize(NumThreads,0);
        
        int nwtot=2*EEN*NumThreads;
        if(nwtot<W.getActiveWalkers())
        {
          app_log() << "  walker number not right. Removing "<< W.getActiveWalkers()-nwtot <<" walkers. ";
          W.destroyWalkers(W.getActiveWalkers()-nwtot);
          vector<int> nwoff(myComm->size()+1,0);
          for(int ip=0; ip<myComm->size(); ip++) nwoff[ip+1]=nwoff[ip]+nwtot;
          W.setGlobalNumWalkers(nwoff[myComm->size()]);
          W.setWalkerOffsets(nwoff);
        }
        else if(nwtot>W.getActiveWalkers())
        {
          app_log() << "  walker number not right. Adding "<< nwtot-W.getActiveWalkers() <<" walkers. ";
          W.createWalkers(nwtot-W.getActiveWalkers());
          vector<int> nwoff(myComm->size()+1,0);
          for(int ip=0; ip<myComm->size(); ip++) nwoff[ip+1]=nwoff[ip]+nwtot;
          W.setGlobalNumWalkers(nwoff[myComm->size()]);
          W.setWalkerOffsets(nwoff);
        }
        FairDivideLow(nwtot,NumThreads,wPerNode);

        app_log() << "  Initial partition of walkers ";
        std::copy(wPerNode.begin(),wPerNode.end(),ostream_iterator<int>(app_log()," "));
        app_log() << endl;

        for(int ip=0; ip<NumThreads; ++ip)
        {
          ostringstream os;
          estimatorClones[ip]= new EstimatorManager(*Estimators);//,*hClones[ip]);
          estimatorClones[ip]->resetTargetParticleSet(*wClones[ip]);
          estimatorClones[ip]->setCollectionMode(false);
          Rng[ip]=new RandomGenerator_t(*(RandomNumberControl::Children[ip]));
          hClones[ip]->setRandomGenerator(Rng[ip]);

          branchClones[ip] = new BranchEngineType(*branchEngine);


          if (QMCDriverMode[QMC_UPDATE_MODE])
          {
            if (UseDrift != "yes")
            {
              os <<"  PbyP moves with drift, using RenyiUpdatePbyP"<<endl;
              RenyiMovers[ip]=new RenyiUpdatePbyP(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip],EEN);
              // Movers[ip]=new VMCUpdatePbyPWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
            }
            else
              APP_ABORT("NOT IMPLEMENTED");
//             else
//             {
//               os <<"  PbyP moves with |psi^2|, using VMCUpdatePbyP"<<endl;
//               Movers[ip]=new VMCUpdatePbyP(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
//             }
            //Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
          }
          else
          {
           APP_ABORT("NOT IMPLEMENTED");
//             if (UseDrift == "yes")
//             {
//               os <<"  walker moves with drift, using VMCUpdateAllWithDrift"<<endl;
//               Movers[ip]=new VMCUpdateAllWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
//             }
//             else
//             {
//               os <<"  walker moves with |psi|^2, using VMCUpdateAll"<<endl;
//               Movers[ip]=new VMCUpdateAll(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
//             }
            //Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
          }

          if(ip==0) app_log() << os.str() << endl;
        }
      }
      
    ParticleSet::ParticlePos_t L(1);
    L.setUnit(PosUnit::LatticeUnit);
    for (int i=0;i<DIM;i++) L[0][i]=1;
    if(W.Lattice.SuperCellEnum!=SUPERCELL_OPEN)
      W.convert2Cart(L);

    ParticleSet::ParticlePos_t C(1);
    C.setUnit(PosUnit::LatticeUnit);
    for (int i=0;i<DIM;i++) C[0][i]=0.5;
    if(W.Lattice.SuperCellEnum!=SUPERCELL_OPEN)
      W.convert2Cart(C);    
    
#if !defined(BGP_BUG)
#pragma omp parallel for
#endif
    for(int ip=0; ip<NumThreads; ++ip)
    {
      //int ip=omp_get_thread_num();
      RenyiMovers[ip]->put(qmcNode);
      RenyiMovers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
    
      MCWalkerConfiguration::iterator  wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
      RenyiMovers[ip]->check_region(wit,wit_end,vsize,computeEE,L,C,Nmax,Nmin);

      if (QMCDriverMode[QMC_UPDATE_MODE])
        RenyiMovers[ip]->initWalkersForPbyP(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
      else
        RenyiMovers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
    }
#pragma omp barrier
//     for(int ip=0; ip<NumThreads; ++ip)
//     {
//       wClones[ip]->clearEnsemble();
//       wClones[ip]->setNumSamples(samples_th[ip]);
//     }

//     myWarmupSteps=0;
    
  }

  bool
  VMCRenyiOMP::put(xmlNodePtr q)
  {
    //nothing to add
    return true;
  }
}

/***************************************************************************
 * $RCSfile: VMCRenyiOMP.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: VMCRenyiOMP.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $
 ***************************************************************************/
