//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/VMC/VMC_CUDA.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Utilities/RandomGenerator.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/CommOperators.h"
#include "QMCDrivers/DriftOperators.h"
#include "type_traits/scalar_traits.h"
#include "Utilities/RunTimeManager.h"
#include "qmc_common.h"

namespace qmcplusplus
{

/// Constructor.
VMCcuda::VMCcuda(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                 QMCHamiltonian& h,WaveFunctionPool& ppool):
  QMCDriver(w,psi,h,ppool),UseDrift("yes"),
  myPeriod4WalkerDump(0), GEVtype("mixed"), w_alpha(0.0), w_beta(0.0), forOpt(false)
{
  RootName = "vmc";
  QMCType ="VMCcuda";
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  QMCDriverMode.set(QMC_WARMUP,0);
  m_param.add(UseDrift,"useDrift","string");
  m_param.add(UseDrift,"usedrift","string");
  m_param.add(nTargetSamples,"targetWalkers","int");
  m_param.add(w_beta,"beta","double");
  m_param.add(w_alpha,"alpha","double");
  m_param.add(GEVtype,"GEVMethod","string");
}

bool VMCcuda::checkBounds (std::vector<PosType> &newpos,
                           std::vector<bool> &valid)
{
  for (int iw=0; iw<newpos.size(); iw++)
  {
    PosType red = W.Lattice.toUnit(newpos[iw]);
    valid[iw] = W.Lattice.isValid(red);
  }
  return true;
}

void VMCcuda::advanceWalkers()
{
  int nat = W.getTotalNum();
  int nw  = W.getActiveWalkers();
  std::vector<PosType>   delpos(nw);
  std::vector<PosType>   newpos(nw);
  std::vector<ValueType> ratios(nw);
#ifdef QMC_COMPLEX
  std::vector<RealType> ratios_real(nw);
#endif
  std::vector<GradType>  newG(nw);
  std::vector<ValueType> newL(nw);
  std::vector<Walker_t*> accepted(nw);

  for (int isub=0; isub<nSubSteps; isub++)
  {
    for(int iat=0; iat<nat; ++iat)
    {
      //create a 3N-Dimensional Gaussian with variance=1
      makeGaussRandomWithEngine(delpos,Random);
      for(int iw=0; iw<nw; ++iw)
      {
        GradType G = W[iw]->G[iat];
        newpos[iw]=W[iw]->R[iat] + m_sqrttau*delpos[iw];
        ratios[iw] = 1.0;
#ifdef QMC_COMPLEX
        ratios_real[iw] = 1.0;
#endif
      }
      W.proposeMove_GPU(newpos, iat);
      Psi.ratio(W,iat,ratios,newG, newL);
#ifdef QMC_COMPLEX
      Psi.convertRatiosFromComplexToReal(ratios, ratios_real);
#endif
      accepted.clear();
      std::vector<bool> acc(nw, true);
      if (W.UseBoundBox)
        checkBounds (newpos, acc);
      for(int iw=0; iw<nw; ++iw)
      {
#ifdef QMC_COMPLEX
        if(acc[iw] && ratios_real[iw]*ratios_real[iw] > Random())
#else
        if(acc[iw] && ratios[iw]*ratios[iw] > Random())
#endif
        {
          accepted.push_back(W[iw]);
          nAccept++;
          W[iw]->R[iat] = newpos[iw];
          acc[iw] = true;
        }
        else
        {
          acc[iw] = false;
          nReject++;
        }
      }
      W.acceptMove_GPU(acc);
      if (accepted.size())
        Psi.update(accepted,iat);
    }
  }
}

bool VMCcuda::run()
{
  if (UseDrift == "yes")
    return runWithDrift();
  resetRun();
  IndexType block = 0;
  IndexType nAcceptTot = 0;
  IndexType nRejectTot = 0;
  IndexType updatePeriod= (QMCDriverMode[QMC_UPDATE_MODE])
                          ? Period4CheckProperties
                          : (nBlocks+1)*nSteps;
  int nat = W.getTotalNum();
  int nw  = W.getActiveWalkers();
  std::vector<RealType>  LocalEnergy(nw);
  Matrix<ValueType> lapl(nw, nat);
  Matrix<GradType>  grad(nw, nat);
  double Esum;

  LoopTimer vmc_loop;
  RunTimeControl runtimeControl(RunTimeManager, MaxCPUSecs);
  bool enough_time_for_next_iteration = true;

  // First do warmup steps
  for (int step=0; step<nWarmupSteps; step++)
    advanceWalkers();

  // Then accumulate statistics
  do
  {
    vmc_loop.start();
    IndexType step = 0;
    nAccept = nReject = 0;
    Esum = 0.0;
    Estimators->startBlock(nSteps);
    do
    {
      ++step;
      ++CurrentStep;
      W.resetCollectables();
      advanceWalkers();
      Psi.gradLapl(W, grad, lapl);
      H.evaluate (W, LocalEnergy);
      if (myPeriod4WalkerDump && (CurrentStep % myPeriod4WalkerDump)==0)
        W.saveEnsemble();
      Estimators->accumulate(W);
    }
    while(step<nSteps);
    if ( nBlocksBetweenRecompute && (1+block)%nBlocksBetweenRecompute == 0 ) Psi.recompute(W);
    if(forOpt)
    {
      d_logpsi_dalpha=0.0;
      d_hpsioverpsi_dalpha=0.0;
      W.copyWalkerGradToGPU();
      Psi.evaluateDerivatives(W, dummy, d_logpsi_dalpha, d_hpsioverpsi_dalpha);
      std::vector<RealType> g_stats(5,0);
      for (int ip=0; ip<nw; ip++)
      {
        RealType E_L = LocalEnergy[ip];
        RealType E_L2= E_L*E_L;
        sE +=E_L;
        sE2+=E_L2;
        sE4+=E_L2*E_L2;
        sW += 1;
        sN+=1;
        for (int i=0; i<numParams; i++)
        {
          RealType di  = d_logpsi_dalpha(ip,i);
          RealType hdi = d_hpsioverpsi_dalpha(ip,i);
          //             vectors
          D_E[i]+= di*E_L;
          HD[i]+=  hdi;
          HD2[i]+= E_L*(hdi+di*E_L);
          D[i]+=   di;
          for (int j=0; j<numParams; j++)
          {
            RealType dj  = d_logpsi_dalpha(ip,j);
            RealType hdj = d_hpsioverpsi_dalpha(ip,j);
            Olp(i,j) += di*dj;
            Ham(i,j) += di*(hdj+dj*E_L);
            Ham2(i,j)+= (hdj+dj*E_L)*(hdi+di*E_L);
          }
        }
      }
      g_stats[0]=sE;
      g_stats[1]=sE2;
      g_stats[2]=sE4;
      g_stats[3]=sW;
      g_stats[4]=sN;
      myComm->allreduce(g_stats);
      RealType nrm = 1.0/g_stats[3];
      E_avg = nrm*g_stats[0];
      V_avg = nrm*g_stats[1]-E_avg*E_avg;
    }
    // std::vector<RealType> logPsi(W.WalkerList.size(), 0.0);
    // Psi.evaluateLog(W, logPsi);
    double accept_ratio = (double)nAccept/(double)(nAccept+nReject);
    Estimators->stopBlock(accept_ratio);
    nAcceptTot += nAccept;
    nRejectTot += nReject;
    ++block;
    recordBlock(block);

    vmc_loop.stop();
    enough_time_for_next_iteration = runtimeControl.enough_time_for_next_iteration(vmc_loop);
    // Rank 0 decides whether the time limit was reached
    myComm->bcast(enough_time_for_next_iteration);
    if (!enough_time_for_next_iteration)
    {
      app_log() << runtimeControl.time_limit_message("VMC", block);
    }

  }
  while(block<nBlocks && enough_time_for_next_iteration);
  //Mover->stopRun();
  //finalize a qmc section
  if (!myComm->rank())
  {
    std::cerr << "At the end of VMC" << std::endl;
    gpu::cuda_memory_manager.report();
  }
  return finalize(block);
}

void VMCcuda::advanceWalkersWithDrift()
{
  int nat = W.getTotalNum();
  int nw  = W.getActiveWalkers();
  std::vector<PosType>   delpos(nw);
  std::vector<PosType>   newpos(nw);
  std::vector<ValueType> ratios(nw);
#ifdef QMC_COMPLEX
  std::vector<RealType>  ratios_real(nw);
#endif
  std::vector<GradType>  oldG(nw), newG(nw);
  std::vector<ValueType> oldL(nw), newL(nw);
  std::vector<Walker_t*> accepted(nw);

  for (int isub=0; isub<nSubSteps; isub++)
  {
    for(int iat=0; iat<nat; iat++)
    {
      Psi.calcGradient (W, iat, oldG);
      //create a 3N-Dimensional Gaussian with variance=1
      makeGaussRandomWithEngine(delpos,Random);
      Psi.addGradient(W, iat, oldG);
      for(int iw=0; iw<nw; iw++)
      {
        PosType dr;
        delpos[iw] *= m_sqrttau;
        getScaledDrift(m_tauovermass,oldG[iw],dr);
        newpos[iw] = W[iw]->R[iat] + delpos[iw] + dr;
        ratios[iw] = 1.0;
#ifdef QMC_COMPLEX
        ratios_real[iw] = 1.0;
#endif
      }
      W.proposeMove_GPU(newpos, iat);
      Psi.calcRatio(W,iat,ratios,newG, newL);
      accepted.clear();
      std::vector<bool> acc(nw, true);
      if (W.UseBoundBox)
        checkBounds (newpos, acc);
      std::vector<RealType> logGf_v(nw);
      std::vector<RealType> rand_v(nw);
      for(int iw=0; iw<nw; ++iw)
      {
        logGf_v[iw] = -m_oneover2tau * dot(delpos[iw], delpos[iw]);
        rand_v[iw] = Random();
      }
      Psi.addRatio(W, iat, ratios, newG, newL);
#ifdef QMC_COMPLEX
      Psi.convertRatiosFromComplexToReal(ratios, ratios_real);
#endif
      for(int iw=0; iw<nw; ++iw)
      {
        PosType drNew;
        getScaledDrift(m_tauovermass,newG[iw],drNew);
        drNew += newpos[iw] - W[iw]->R[iat];
        // if (dot(drNew, drNew) > 25.0)
        //   std::cerr << "Large drift encountered!  Drift = " << drNew << std::endl;
        RealType logGb = -m_oneover2tau * dot(drNew, drNew);
        RealType x = logGb - logGf_v[iw];
#ifdef QMC_COMPLEX
        RealType prob = ratios_real[iw]*ratios_real[iw]*std::exp(x);
#else
        RealType prob = ratios[iw]*ratios[iw]*std::exp(x);
#endif
        if(acc[iw] && rand_v[iw] < prob)
        {
          accepted.push_back(W[iw]);
          nAccept++;
          W[iw]->R[iat] = newpos[iw];
          acc[iw] = true;
        }
        else
        {
          acc[iw] = false;
          nReject++;
        }
      }
      W.acceptMove_GPU(acc);
      if (accepted.size())
        Psi.update(accepted,iat);
    }
  }
}

bool VMCcuda::runWithDrift()
{
  resetRun();
  IndexType block = 0;
  IndexType nAcceptTot = 0;
  IndexType nRejectTot = 0;
  int nat = W.getTotalNum();
  int nw  = W.getActiveWalkers();
  std::vector<RealType>  LocalEnergy(nw);
  Matrix<ValueType> lapl(nw, nat);
  Matrix<GradType>  grad(nw, nat);

  // First, do warmup steps
  for (int step=0; step<nWarmupSteps; step++)
    advanceWalkersWithDrift();

  // Now do data collection steps
  do
  {
    IndexType step = 0;
    nAccept = nReject = 0;
    Estimators->startBlock(nSteps);
    do
    {
      step++;
      CurrentStep++;
      W.resetCollectables();
      advanceWalkersWithDrift();
      Psi.gradLapl(W, grad, lapl);
      H.evaluate (W, LocalEnergy);
      if (myPeriod4WalkerDump && (CurrentStep % myPeriod4WalkerDump)==0)
        W.saveEnsemble();
      Estimators->accumulate(W);
    }
    while(step<nSteps);
    if ( nBlocksBetweenRecompute && (1+block)%nBlocksBetweenRecompute == 0 ) Psi.recompute(W);
    if(forOpt)
    {
      d_logpsi_dalpha=0.0;
      d_hpsioverpsi_dalpha=0.0;
      W.copyWalkerGradToGPU();
      Psi.evaluateDerivatives(W, dummy, d_logpsi_dalpha, d_hpsioverpsi_dalpha);
      std::vector<RealType> g_stats(5,0);
      for (int ip=0; ip<nw; ip++)
      {
        RealType E_L = LocalEnergy[ip];
        RealType E_L2= E_L*E_L;
        sE +=E_L;
        sE2+=E_L2;
        sE4+=E_L2*E_L2;
        sW += 1;
        sN+=1;
        for (int i=0; i<numParams; i++)
        {
          RealType di  = d_logpsi_dalpha(ip,i);
          RealType hdi = d_hpsioverpsi_dalpha(ip,i);
          //             vectors
          D_E[i]+= di*E_L;
          HD[i]+=  hdi;
          HD2[i]+= E_L*(hdi+di*E_L);
          D[i]+=   di;
          for (int j=0; j<numParams; j++)
          {
            RealType dj  = d_logpsi_dalpha(ip,j);
            RealType hdj = d_hpsioverpsi_dalpha(ip,j);
            Olp(i,j) += di*dj;
            Ham(i,j) += di*(hdj+dj*E_L);
            Ham2(i,j)+= (hdj+dj*E_L)*(hdi+di*E_L);
          }
        }
      }
      g_stats[0]=sE;
      g_stats[1]=sE2;
      g_stats[2]=sE4;
      g_stats[3]=sW;
      g_stats[4]=sN;
      myComm->allreduce(g_stats);
      RealType nrm = 1.0/g_stats[3];
      E_avg = nrm*g_stats[0];
      V_avg = nrm*g_stats[1]-E_avg*E_avg;
      //for (int i=0; i<numParams; i++) app_log()<<HD[i]<<" ";
      //  app_log()<< std::endl;
    }
//      Psi.recompute(W);
    double accept_ratio = (double)nAccept/(double)(nAccept+nReject);
    Estimators->stopBlock(accept_ratio);
    nAcceptTot += nAccept;
    nRejectTot += nReject;
    ++block;
    recordBlock(block);
  }
  while(block<nBlocks);
  //finalize a qmc section
  if (!myComm->rank())
  {
    std::cerr << "At the end of VMC with drift" << std::endl;
    gpu::cuda_memory_manager.report();
  }
  return finalize(block);
}





void VMCcuda::resetRun()
{
  SpeciesSet tspecies(W.getSpeciesSet());
  int massind=tspecies.addAttribute("mass");
  RealType mass = tspecies(massind,0);
  RealType oneovermass = 1.0/mass;
  RealType oneoversqrtmass = std::sqrt(oneovermass);
  m_oneover2tau = 0.5*mass/Tau;
  m_sqrttau = std::sqrt(Tau/mass);
  m_tauovermass = Tau/mass;
  if (!myComm->rank())
  {
    std::cerr << "Before allocating GPU buffer" << std::endl;
    gpu::cuda_memory_manager.report();
  }
  // Compute the size of data needed for each walker on the GPU card
  PointerPool<Walker_t::cuda_Buffer_t > pool;
  app_log() << "Starting VMCcuda::resetRun() " << std::endl;
  Psi.reserve (pool);
  app_log() << "Each walker requires " << pool.getTotalSize() * sizeof(CudaValueType)
            << " bytes in GPU memory.\n";
  // Now allocate memory on the GPU card for each walker
  // for (int iw=0; iw<W.WalkerList.size(); iw++) {
  //   Walker_t &walker = *(W.WalkerList[iw]);
  //   walker.resizeCuda(pool.getTotalSize());
  //   // pool.allocate(walker.cuda_DataSet);
  // }
  W.allocateGPU(pool.getTotalSize());
  if (!myComm->rank())
  {
    std::cerr << "After allocating GPU buffer" << std::endl;
    gpu::cuda_memory_manager.report();
  }
  app_log() << "Successfully allocated walkers.\n";
  W.copyWalkersToGPU();
  W.updateLists_GPU();
  std::vector<RealType> logPsi(W.WalkerList.size(), 0.0);
  //Psi.evaluateLog(W, logPsi);
  Psi.recompute(W, true);
  Estimators->start(nBlocks, true);
  // Compute sample dumping frequency
  if (nTargetSamples>(myComm->size()*W.WalkerList.size()))
  {
    int nnodes = myComm->size();
    int nw = W.WalkerList.size();
    int samples_per_node = (nTargetSamples+nnodes-1)/nnodes;
    int dumps_per_node   = (samples_per_node+nw-1) / nw;
    myPeriod4WalkerDump = Period4WalkerDump;
    app_log() << "  Dumping walker ensemble every " << myPeriod4WalkerDump
              << " steps.\n";
  }
  else
  {
    myPeriod4WalkerDump=nBlocks*nSteps;
    nTargetSamples=myComm->size()*W.WalkerList.size();
  }
  W.clearEnsemble();
  int samples_this_node = nTargetSamples/myComm->size();
  if (nTargetSamples%myComm->size() > myComm->rank()) samples_this_node+=1;
  app_log() << "  Node zero will generate " << samples_this_node << " samples.\n";
  W.setNumSamples(samples_this_node);
  if(forOpt)
  {
    dummy.clear();
    Psi.checkInVariables(dummy);
    dummy.resetIndex();
    Psi.checkOutVariables(dummy);
    numParams = dummy.size();
    resizeForOpt(numParams);
    int nw = W.getActiveWalkers();
    d_logpsi_dalpha.resize(nw, numParams);
    d_hpsioverpsi_dalpha.resize(nw, numParams);
    d_logpsi_dalpha=0.0;
    d_hpsioverpsi_dalpha=0.0;
  }
}

bool
VMCcuda::put(xmlNodePtr q)
{
  app_log() << "\n<vmc function=\"put\">"
    << "\n  qmc_counter=" << qmc_common.qmc_counter << "  my_counter=" << MyCounter<< std::endl;
  if(qmc_common.qmc_counter && MyCounter)
  {
    nSteps=prevSteps;
    nStepsBetweenSamples=prevStepsBetweenSamples;
  }
  else
  {
    int nw=W.getActiveWalkers();
    //compute samples and overwrite steps for the given samples
    int Nthreads = 1;
    int Nprocs=myComm->size();
    //target samples set by samples or samplesperthread/dmcwalkersperthread
    nTargetPopulation=std::max(nTargetPopulation,nSamplesPerThread*Nprocs*Nthreads);
    nTargetSamples=static_cast<int>(std::ceil(nTargetPopulation));
    if(nTargetSamples)
    {
      int nwtot=nw*Nprocs;  //total number of walkers used by this qmcsection
      nTargetSamples=std::max(nwtot,nTargetSamples);
      nTargetSamples=((nTargetSamples+nwtot-1)/nwtot)*nwtot; // nTargetSamples are always multiples of total number of walkers
      nSamplesPerThread=nTargetSamples/Nprocs/Nthreads;
      int ns_target=nTargetSamples*nStepsBetweenSamples; //total samples to generate
      int ns_per_step=Nprocs*nw;  //total samples per step
      nSteps=std::max(nSteps,(ns_target/ns_per_step+nBlocks-1)/nBlocks);
      Period4WalkerDump=nStepsBetweenSamples=ns_per_step*nSteps*nBlocks/nTargetSamples;
    }
    else
    {
      Period4WalkerDump = nStepsBetweenSamples=(nBlocks+1)*nSteps; //some positive number, not used
      nSamplesPerThread=0;
    }
  }
  prevSteps=nSteps;
  prevStepsBetweenSamples=nStepsBetweenSamples;

  app_log() << "  time step      = " << Tau << std::endl;
  app_log() << "  blocks         = " << nBlocks << std::endl;
  app_log() << "  steps          = " << nSteps << std::endl;
  app_log() << "  substeps       = " << nSubSteps << std::endl;
  app_log() << "  current        = " << CurrentStep << std::endl;
  app_log() << "  target samples = " << nTargetPopulation << std::endl;
  app_log() << "  walkers/mpi    = " << W.getActiveWalkers() << std::endl << std::endl;
  app_log() << "  stepsbetweensamples = " << nStepsBetweenSamples << std::endl;

  if(DumpConfig)
    app_log() << "  DumpConfig==true Configurations are dumped to config.h5 with a period of " << Period4CheckPoint << " blocks" << std::endl;
  else
    app_log() << "  DumpConfig==false Nothing (configurations, state) will be saved." << std::endl;
  if (Period4WalkerDump>0)
    app_log() << "  Walker Samples are dumped every " << Period4WalkerDump << " steps." << std::endl;
  app_log() << "</vmc>" << std::endl;
  app_log().flush();
  //nothing to add
  return true;
}

VMCcuda::RealType VMCcuda::fillOverlapHamiltonianMatrices(Matrix<RealType>& LeftM, Matrix<RealType>& RightM)
{
  RealType b1,b2;
  if (GEVtype=="H2")
  {
    b1=w_beta;
    b2=w_alpha;
  }
  else
  {
    b2=w_beta;
    b1=0;
  }
  std::vector<RealType> g_stats(5,0);
  g_stats[0]=sE;
  g_stats[1]=sE2;
  g_stats[2]=sE4;
  g_stats[3]=sW;
  g_stats[4]=sN;
  myComm->allreduce(g_stats);
  RealType g_nrm = 1.0/g_stats[3];
  E_avg = g_nrm*g_stats[0];
  RealType E_avg2=E_avg*E_avg;
  RealType E2_avg = g_nrm*g_stats[1];
  V_avg = E2_avg-E_avg2;
  myComm->allreduce(Ham2);
  Ham2*=g_nrm;
  myComm->allreduce(Ham);
  Ham *=g_nrm;
  myComm->allreduce(Olp);
  Olp *=g_nrm;
  myComm->allreduce(D_E);
  myComm->allreduce(D);
  myComm->allreduce(HD);
  myComm->allreduce(HD2);
  for (int i=0; i<numParams; i++)
  {
    D_E[i]*=g_nrm;
    D[i]*=g_nrm;
    HD[i]*=g_nrm;
    HD2[i]*=g_nrm;
  }
  for (int i=0; i<numParams; i++)
    for (int j=0; j<numParams; j++)
      Ham(i,j) += -D[i]*(HD[j]+ D_E[j] - D[j]*E_avg) - D[j]*D_E[i];
  for (int i=0; i<numParams; i++)
    for (int j=0; j<numParams; j++)
      Olp(i,j) -= D[i]*D[j];
  for (int i=0; i<numParams; i++)
    for (int j=0; j<numParams; j++)
      Ham2(i,j) += D[i]*D[j]*E2_avg - D[i]*HD2[j] - D[j]*HD2[i]  - E_avg*(Ham(i,j)+Ham(j,i))+ 2.0*E_avg2*Olp(i,j) ;
  RealType b1_rat = b1/E_avg2;
  for (int i=1; i<numParams+1; i++)
    for (int j=1; j<numParams+1; j++)
    {
      LeftM(i,j) = (1-b2)*Ham(i-1,j-1) + b2*(Ham2(i-1,j-1) - E_avg2*Olp(i-1,j-1)) ;
      RightM(i,j) = (1.0-b1)*Olp(i-1,j-1) + b1_rat*Ham2(i-1,j-1);
    }
  RightM(0,0)= 1.0-b1 + b1_rat*E_avg2;
  LeftM(0,0)=(1-b2)*E_avg+b2*V_avg;
  for (int i=1; i<numParams+1; i++)
  {
    RightM(0,i)= RightM(i,0) = b1_rat*(HD2[i-1] - E_avg*(HD[i-1]+ 2.0*(D_E[i-1]-D[i-1]*E_avg)) - D[i-1]*E2_avg);
//       RightM(0,i)= RightM(i,0) = b1_rat*(HD2[i-1] -  E_avg*HD[i-1] -  D[i-1]*E2_avg  - 2.0*E_avg*(D_E[i-1]-D[i-1]*E_avg ));
    LeftM(i,0) = (1-b2)*(D_E[i-1]-E_avg*D[i-1])         +b2*(HD2[i-1] - E_avg*(HD[i-1]+ 2.0*(D_E[i-1]-D[i-1]*E_avg)) - D[i-1]*E2_avg);
    LeftM(0,i) = (1-b2)*(HD[i-1]+D_E[i-1]-E_avg*D[i-1]) +b2*(HD2[i-1] - E_avg*(HD[i-1]+ 2.0*(D_E[i-1]-D[i-1]*E_avg)) - D[i-1]*E2_avg);
  }
  /*
      for (int i=0; i<numParams; i++) app_log()<<D[i]<<" ";
      app_log()<< std::endl;
      for (int i=0; i<numParams; i++) app_log()<<D_E[i]<<" ";
      app_log()<< std::endl;
      for (int i=0; i<numParams; i++) app_log()<<HD[i]<<" ";
      app_log()<< std::endl;
      for (int i=0; i<numParams; i++) app_log()<<HD2[i]<<" ";
      app_log()<< std::endl;

      for (int i=0; i<numParams+1; i++)
      {
        for (int j=0; j<numParams+1; j++)
          app_log()<<LeftM(i,j)<<" ";
        app_log()<< std::endl;
      }
      app_log()<< std::endl;
      for (int i=0; i<numParams+1; i++)
      {
        for (int j=0; j<numParams+1; j++)
          app_log()<<RightM(i,j)<<" ";
        app_log()<< std::endl;
      }
      app_log()<< std::endl;
  */
  return 1.0;
}

}

