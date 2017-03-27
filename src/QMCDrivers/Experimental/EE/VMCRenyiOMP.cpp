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
  myWarmupSteps(0),UseDrift("yes"), computeEE("sphere"),vsize(1),EEN(2),Nmax(20),Nmin(0),plotSwapAmplitude(false),grid_spacing(0)
{
  RootName = "vmc";
  QMCType ="VMCRenyiOMP";
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  QMCDriverMode.set(QMC_WARMUP,0);
  m_param.add(UseDrift,"useDrift","string");
  m_param.add(UseDrift,"usedrift","string");
  m_param.add(UseDrift,"use_drift","string");
  m_param.add(computeEE,"computeEE","string");
  m_param.add(vsize,"EEsize","double");
  m_param.add(EEN,"EEN","int");
  m_param.add(Nmin,"Nmin","int");
  m_param.add(Nmax,"Nmax","int");
  m_param.add(plotSwapAmplitude,"plot","int");
  m_param.add(grid_spacing,"grid","int");
  m_param.add(myWarmupSteps,"warmupSteps","int");
  m_param.add(myWarmupSteps,"warmupsteps","int");
  m_param.add(myWarmupSteps,"warmup_steps","int");
  m_param.add(nTargetSamples,"targetWalkers","int");
  m_param.add(nTargetSamples,"targetwalkers","int");
  m_param.add(nTargetSamples,"target_walkers","int");
}

bool VMCRenyiOMP::run()
{
  resetRun();
  int total_grid=2*grid_spacing+1;
  int cntr=(total_grid+1)/2;
  int n_grid_pe(grid_spacing*grid_spacing+1);
  if(myComm->rank()==0)
  {
    ee_dat.str("");
    ee_dat<<RootName.c_str()<<".ee.dat";
    file_out.open(ee_dat.str().c_str(),fstream::out | fstream::trunc);
    file_out<<"# indx";
    file_out<<"  EE_S";
    for(int i(Nmin); i<Nmax+1; i++)
      file_out<<"  N"<<i;
    file_out<< std::endl;
    file_out.close();
    if(plotSwapAmplitude)
    {
      pe_dat.str("");
      pe_dat<<RootName.c_str()<<".pe.dat";
      file_out.open(pe_dat.str().c_str(),fstream::out | fstream::trunc);
      file_out<<"# indx";
      for(int i(0); i<n_grid_pe; i++)
        file_out<<"  S_"<<i;
      file_out<< std::endl;
    }
  }
  std::vector<Matrix<RealType> > averageSwaps(NumThreads,Matrix<RealType>(total_grid,total_grid));
  //start the main estimator
  Estimators->start(nBlocks);
  for (int ip=0; ip<NumThreads; ++ip)
    RenyiMovers[ip]->startRun(nBlocks,false);
  const bool has_collectables=W.Collectables.size();
  hpmStart(QMC_VMC_0_EVENT,"vmc::main");
  for (int block=0; block<nBlocks; ++block)
  {
    #pragma omp parallel
    {
      int ip=omp_get_thread_num();
      MCWalkerConfiguration::iterator wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
      RenyiMovers[ip]->startBlock(nSteps);
      for (int step=0; step<nSteps; ++step)
        RenyiMovers[ip]->advanceWalkers(wit,wit_end,false);
      wit=W.begin()+wPerNode[ip];
      RenyiMovers[ip]->accumulate(wit,wit_end);
      if(plotSwapAmplitude)
      {
        wit=W.begin()+wPerNode[ip];
        RenyiMovers[ip]->plotSwapAmplitude(wit,wit_end,averageSwaps[ip]);
      }
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
    if(plotSwapAmplitude)
    {
      for (int ip=1; ip<NumThreads; ++ip)
        averageSwaps[0]+=averageSwaps[ip];
      myComm->allreduce(averageSwaps[0]);
      averageSwaps[0]*=nrm;
      std::vector<RealType> pe_r(n_grid_pe,0);
      std::vector<int> n_pe_r(n_grid_pe,0);
      for (int i(0); i<total_grid; i++)
        for (int j(0); j<total_grid; j++)
        {
          int z=(i-cntr)*(i-cntr)+(j-cntr)*(j-cntr);
          if (z<n_grid_pe)
          {
            pe_r[z]+=averageSwaps[0](i,j);
            n_pe_r[z]+=1;
          }
        }
      for (int i(0); i<n_grid_pe; i++)
        if(n_pe_r[i]>0)
          pe_r[i]*=1.0/n_pe_r[i];
      if(myComm->rank()==0)
      {
        file_out.open(pe_dat.str().c_str(),fstream::out | fstream::app);
        file_out<<block<<" ";
        for(int i(0); i<n_grid_pe; i++)
          file_out<<pe_r[i]<<" ";
        file_out<< std::endl;
        file_out.close();
      }
      for (int ip=0; ip<NumThreads; ++ip)
        averageSwaps[ip]=0;
    }
    if(myComm->rank()==0)
    {
      file_out.open(ee_dat.str().c_str(),fstream::out | fstream::app);
      file_out<<block<<" ";
      file_out<<avgsgn*nrm<<" ";
      for(int i(0); i<(Nmax+1-Nmin); i++)
        file_out<<n_stats[i]*nrm<<" ";
      file_out<< std::endl;
      file_out.close();
    }
//       for (int ip=0; ip<NumThreads; ++ip)
//         RenyiMovers[ip]->print_all();
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
  app_log() << "  Warmup Steps " << myWarmupSteps << std::endl;
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
      std::vector<int> nwoff(myComm->size()+1,0);
      for(int ip=0; ip<myComm->size(); ip++)
        nwoff[ip+1]=nwoff[ip]+nwtot;
      W.setGlobalNumWalkers(nwoff[myComm->size()]);
      W.setWalkerOffsets(nwoff);
    }
    else
      if(nwtot>W.getActiveWalkers())
      {
        app_log() << "  walker number not right. Adding "<< nwtot-W.getActiveWalkers() <<" walkers. ";
        W.createWalkers(nwtot-W.getActiveWalkers());
        std::vector<int> nwoff(myComm->size()+1,0);
        for(int ip=0; ip<myComm->size(); ip++)
          nwoff[ip+1]=nwoff[ip]+nwtot;
        W.setGlobalNumWalkers(nwoff[myComm->size()]);
        W.setWalkerOffsets(nwoff);
      }
    FairDivideLow(nwtot,NumThreads,wPerNode);
    app_log() << "  Initial partition of walkers ";
    copy(wPerNode.begin(),wPerNode.end(),std::ostream_iterator<int>(app_log()," "));
    app_log() << std::endl;
    for(int ip=0; ip<NumThreads; ++ip)
    {
      std::ostringstream os;
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
          os <<"  PbyP moves without drift, using RenyiUpdatePbyP"<< std::endl;
          RenyiMovers[ip]=new RenyiUpdatePbyP(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip],EEN);
        }
        else
          APP_ABORT("NOT IMPLEMENTED");
      }
      else
      {
        APP_ABORT("NOT IMPLEMENTED");
//             if (UseDrift != "yes")
//             {
//               os <<"  Walkers moves, using RenyiUpdatePbyP"<< std::endl;
//               RenyiMovers[ip]=new RenyiUpdateAll(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip],EEN);
//             }
//             else
//             {
//               os <<"  Walkers moves with drift, using RenyiUpdatePbyP"<< std::endl;
//               RenyiMovers[ip]=new RenyiUpdateAllWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip],EEN);
//             }
      }
      if(ip==0)
        app_log() << os.str() << std::endl;
    }
  }
  ParticleSet::ParticlePos_t L(1);
  L.setUnit(PosUnit::LatticeUnit);
  for (int i=0; i<DIM; i++)
    L[0][i]=1;
  if(W.Lattice.SuperCellEnum!=SUPERCELL_OPEN)
    W.convert2Cart(L);
  ParticleSet::ParticlePos_t C(1);
  C.setUnit(PosUnit::LatticeUnit);
  for (int i=0; i<DIM; i++)
    C[0][i]=0.5;
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
    RenyiMovers[ip]->check_region(wit,wit_end,vsize,computeEE,L,C,Nmax,Nmin,QMCDriverMode[QMC_UPDATE_MODE]);
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

