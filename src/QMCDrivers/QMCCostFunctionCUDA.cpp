//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    same
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/QMCCostFunctionCUDA.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Particle/HDFWalkerInputCollect.h"
#include "Message/CommOperators.h"

namespace qmcplusplus
{

QMCCostFunctionCUDA::QMCCostFunctionCUDA
( MCWalkerConfiguration& w, TrialWaveFunction& psi,
  QMCHamiltonian& h, HamiltonianPool& hpool):
  QMCCostFunctionBase(w,psi,h), CloneManager(hpool)
{
}


/** Clean up the vector */
QMCCostFunctionCUDA::~QMCCostFunctionCUDA()
{
  delete_iter(RecordsOnNode.begin(),RecordsOnNode.end());
  delete_iter(DerivRecords.begin(),DerivRecords.end());
  delete_iter(HDerivRecords.begin(),HDerivRecords.end());
}


/**  Perform the correlated sampling algorthim.
 */
QMCCostFunctionCUDA::Return_t QMCCostFunctionCUDA::correlatedSampling(bool needDerivs)
{
  Return_t wgt_tot=0.0;
  int nw = W.getActiveWalkers();
  //#pragma omp parallel reduction(+:wgt_tot)
  Return_t eloc_new;
  //int totalElements=W.getTotalNum()*OHMMS_DIM;
  MCWalkerConfiguration::iterator it(W.begin());
  MCWalkerConfiguration::iterator it_end(W.end());
  int iw=0;
  int numParams = OptVariablesForPsi.size();
  int numPtcl   = W.getTotalNum();
  std::vector<RealType> logpsi_new(nw), logpsi_fixed(nw), KE(nw);
  RealMatrix_t d_logpsi_dalpha, d_hpsioverpsi_dalpha;
  GradMatrix_t  fixedG(nw, numPtcl);
  ValueMatrix_t fixedL(nw, numPtcl);
  GradMatrix_t  newG(nw, numPtcl);
  ValueMatrix_t newL(nw, numPtcl);
  if (needDerivs)
  {
    d_logpsi_dalpha.resize(nw, numParams);
    d_hpsioverpsi_dalpha.resize(nw, numParams);
  }
  Psi.evaluateOptimizableLog(W, logpsi_new, newG, newL);
  //    RealType factor = samplePsi2 ? 2.0 : 1.0;
  RealType factor = 2.0;
  // Add optimizable and non-optimizable gradients and Laplacians
  for (int iw=0; iw<nw; iw++)
  {
    ParticleSet::Walker_t& walker = *(W[iw]);
    for (int iat=0; iat<numPtcl; iat++)
    {
      walker.G[iat] = newG(iw, iat) + dlogPsi_fixed(iw, iat);
      walker.L[iat] = newL(iw, iat) + d2logPsi_fixed(iw, iat);
    }
    // W.R = walker.R;
    // W.update();
    // Return_t* restrict saved= &(Records(iw,0));
    // Psi.evaluateDeltaLog(W, saved[LOGPSI_FIXED], saved[LOGPSI_FREE],
    // 			   *dLogPsi[iw],*d2LogPsi[iw]);
    // Psi.evaluateDeltaLog(W);
    // fprintf (stderr, "iw=%4d CPU deltaLog = %12.6e  dLog = [%12.6e %12.6e %12.6e]  d2Log=%12.6e\n",
    // 	       iw, saved[LOGPSI_FREE], W.G[0][0], W.G[0][1], W.G[0][2],
    // 	       W.L[0]);
    // fprintf (stderr, "iw=%4d GPU deltaLog = %12.6e  dLog = [%12.6e %12.6e %12.6e]  d2Log=%12.6e\n",
    // 	       iw, logpsi_new[iw], newG(iw,0)[0], newG(iw,0)[1], newG(iw,0)[2],
    // 	       newL(iw,0));
    // fprintf (stderr, "iw=%4d CPU deltaLog = %12.6e  dLog = [%12.6e %12.6e %12.6e]  d2Log=%12.6e\n",
    // 	       iw, saved[LOGPSI_FREE], (*dLogPsi[iw])[0][0],(*dLogPsi[iw])[0][1],(*dLogPsi[iw])[0][2],
    // 	       (*d2LogPsi[iw])[0]);
    // fprintf (stderr, "iw=%4d GPU deltaLog = %12.6e  dLog = [%12.6e %12.6e %12.6e]  d2Log=%12.6e\n",
    // 	       iw, logpsi_new[iw],dlogPsi_fixed(iw,0)[0],dlogPsi_fixed(iw,0)[1],dlogPsi_fixed(iw,0)[2],
    // 	       d2logPsi_fixed(iw,0));
  }
  W.copyWalkersToGPU(needDerivs);
  H_KE.evaluate (W, KE);
  if (needDerivs)
    Psi.evaluateDerivatives(W, OptVariablesForPsi,
                            d_logpsi_dalpha, d_hpsioverpsi_dalpha);
  wgt_tot = 0.0;
  for (int iw=0; iw<nw; iw++)
  {
    ParticleSet::Walker_t& walker = *(W[iw]);
    Return_t* restrict saved= &(Records(iw,0));
    RealType weight = factor*(logpsi_new[iw] - saved[LOGPSI_FREE]);
    RealType eloc_new = KE[iw] + saved[ENERGY_FIXED];
    saved[ENERGY_NEW] = eloc_new;
    saved[REWEIGHT]   = weight;
    wgt_tot += weight;
    if (needDerivs)
      for (int ip=0; ip<NumOptimizables; ip++)
      {
        TempDerivRecords[iw][ip]  =      d_logpsi_dalpha(iw,ip);
        TempHDerivRecords[iw][ip] = d_hpsioverpsi_dalpha(iw,ip);
      }
  }
  //this is MPI barrier
  //OHMMS::Controller->barrier();
  //collect the total weight for normalization and apply maximum weight
  myComm->allreduce(wgt_tot);
  
  Return_t wgtnorm = wgt_tot/NumSamples;
  wgt_tot=0.0;
  for (int iw=0; iw<nw; iw++)
  {
    Return_t* restrict saved = &(Records(iw,0));
    saved[REWEIGHT] = std::min(std::exp(saved[REWEIGHT]-wgtnorm), std::numeric_limits<Return_t>::max()*0.1 );
    wgt_tot+= saved[REWEIGHT];
  }
  myComm->allreduce(wgt_tot);

  wgtnorm =(wgt_tot==0)?1:(1.0*NumSamples)/wgt_tot;
  wgt_tot=0.0;
  for (int iw=0; iw<nw; iw++)
  {
    Return_t* restrict saved = Records[iw];
    saved[REWEIGHT] = std::min(saved[REWEIGHT]*wgtnorm,MaxWeight) ;
//app_log()<<saved[REWEIGHT]<< std::endl;
    wgt_tot+= saved[REWEIGHT];
  }
  myComm->allreduce(wgt_tot);
  for (int i=0; i<SumValue.size(); i++)
    SumValue[i]=0.0;
  CSWeight=wgt_tot=1.0/wgt_tot;
  for (int iw=0; iw<nw; iw++)
  {
    const Return_t* restrict saved = &(Records(iw,0));
    Return_t weight=saved[REWEIGHT]*wgt_tot;
    Return_t eloc_new=saved[ENERGY_NEW];
    Return_t delE=std::pow(std::abs(eloc_new-EtargetEff),PowerE);
    SumValue[SUM_E_BARE]    += eloc_new;
    SumValue[SUM_ESQ_BARE]  += eloc_new*eloc_new;
    SumValue[SUM_ABSE_BARE] += delE;
    SumValue[SUM_E_WGT]     += eloc_new*saved[REWEIGHT];
    SumValue[SUM_ESQ_WGT]   += eloc_new*eloc_new*saved[REWEIGHT];
    SumValue[SUM_ABSE_WGT]  += delE*saved[REWEIGHT];
    SumValue[SUM_WGT]       += saved[REWEIGHT];
    SumValue[SUM_WGTSQ]     += saved[REWEIGHT]*saved[REWEIGHT];
  }
  //collect everything
  myComm->allreduce(SumValue);
  RealType effective_walkers = SumValue[SUM_WGT]*SumValue[SUM_WGT]/SumValue[SUM_WGTSQ];
  return SumValue[SUM_WGT]*SumValue[SUM_WGT]/SumValue[SUM_WGTSQ];
}

void
QMCCostFunctionCUDA::getConfigurations(const std::string& aroot)
{
  app_log() << "   Loading configuration from MCWalkerConfiguration::SampleStack " << std::endl;
  app_log() << "    number of walkers before load " << W.getActiveWalkers() << std::endl;
  if (H_KE.size()==0)
  {
    H_KE.addOperator(H.getHamiltonian("Kinetic"),"Kinetic");
    if (includeNonlocalH != "no")
    {
      if (includeNonlocalH=="yes")
        includeNonlocalH="NonLocalECP";
      QMCHamiltonianBase* a=H.getHamiltonian(includeNonlocalH);
      if(a)
      {
        H_KE.addOperator(a,includeNonlocalH);
      }
    }
    H_KE.addObservables(W);
  }

  outputManager.pause();
  W.loadEnsemble();
  outputManager.resume();
  app_log() << "    number of walkers after load: "
            << W.getActiveWalkers() << std::endl;
  if(dLogPsi.size() != W.getActiveWalkers())
  {
    delete_iter(dLogPsi.begin(),dLogPsi.end());
    delete_iter(d2LogPsi.begin(),d2LogPsi.end());
    int nptcl=W.getTotalNum();
    int nwtot=W.getActiveWalkers();
    dLogPsi.resize(nwtot,0);
    d2LogPsi.resize(nwtot,0);
    for(int i=0; i<nwtot; ++i)
      dLogPsi[i] =new ParticleGradient_t (nptcl);
    for(int i=0; i<nwtot; ++i)
      d2LogPsi[i]=new ParticleLaplacian_t(nptcl);
  }
  PointerPool<Walker_t::cuda_Buffer_t > pool;
  // Reserve memory only for optimizable parts of the wavefunction
  if (includeNonlocalH != "no")
    Psi.reserve (pool, false);
  else
    Psi.reserve (pool, true);
  app_log() << "Each walker requires " << pool.getTotalSize() * sizeof(CudaRealType)
            << " bytes in GPU memory.\n";
  // for (int iw=0; iw<W.WalkerList.size(); iw++) {
  //   Walker_t &walker = *(W.WalkerList[iw]);
  //   walker.resizeCuda(pool.getTotalSize());
  //   //pool.allocate(walker.cuda_DataSet);
  // }
  W.allocateGPU(pool.getTotalSize());
  W.copyWalkersToGPU();
  W.updateLists_GPU();
  app_log() << "Successfully allocated walkers.\n";
}



/** evaluate everything before optimization */
void
QMCCostFunctionCUDA::checkConfigurations()
{
  RealType et_tot=0.0;
  RealType e2_tot=0.0;
  int numWalkers=W.getActiveWalkers();
  Records.resize(numWalkers,NUMPROPERTIES);
  int nptcl = W.getTotalNum();
  dlogPsi_opt.resize (numWalkers, nptcl);
  d2logPsi_opt.resize(numWalkers, nptcl);
  dlogPsi_fixed.resize (numWalkers, nptcl);
  d2logPsi_fixed.resize(numWalkers, nptcl);
  //if (needGrads)
  {
    TempHDerivRecords.resize(numWalkers,std::vector<Return_t>(NumOptimizables,0));
    TempDerivRecords.resize(numWalkers,std::vector<Return_t>(NumOptimizables,0));
    LogPsi_Derivs.resize(numWalkers, NumOptimizables);
    LocE_Derivs.resize  (numWalkers, NumOptimizables);
  }
  Return_t e0=0.0;
  Return_t e2=0.0;
  std::vector<RealType> logPsi_free(numWalkers,0), d2logPsi_free(numWalkers,0);
  std::vector<RealType> logPsi_fixed(numWalkers,0);
  std::vector<GradType> dlogPsi_free(numWalkers);
  //Psi.evaluateDeltaLog(W, logPsi_free);
  if (includeNonlocalH != "no")
    Psi.evaluateDeltaLog(W, logPsi_fixed, logPsi_free, dlogPsi_fixed, d2logPsi_fixed);
  else
    Psi.evaluateOptimizableLog (W, logPsi_free, dlogPsi_opt, d2logPsi_opt);
  if (needGrads)
    Psi.evaluateDerivatives (W, OptVariables, LogPsi_Derivs, LocE_Derivs);
  int nat = W.getTotalNum();
  MCWalkerConfiguration::iterator it(W.begin());
  MCWalkerConfiguration::iterator it_end(W.end());
  for(int iw=0; it!=it_end; ++it,++iw)
  {
    Walker_t& w(**it);
    RealType *prop = w.getPropertyBase();
    Return_t* restrict saved= &(Records(iw,0));
    saved[ENERGY_TOT] = saved[ENERGY_NEW] = prop[LOCALENERGY];
    saved[ENERGY_FIXED] = prop[LOCALPOTENTIAL];
    saved[LOGPSI_FIXED] = prop[LOGPSI] - logPsi_free[iw];
    saved[LOGPSI_FREE]  = logPsi_free[iw];
    e0 += prop[LOCALENERGY];
    e2 += prop[LOCALENERGY] * prop[LOCALENERGY];
    if (needGrads)
      for (int ip=0; ip<OptVariables.size(); ip++)
      {
        TempDerivRecords[iw][ip] = LogPsi_Derivs(iw,ip);
        TempHDerivRecords[iw][ip] = LocE_Derivs(iw,ip);
      }
    if (includeNonlocalH != "no")
      for (int iat=0; iat<nat; iat++)
      {
        dlogPsi_opt(iw, iat)  = w.G[iat] - dlogPsi_fixed(iw,iat);
        d2logPsi_opt(iw, iat) = w.L[iat] - d2logPsi_fixed(iw,iat);
      }
    else
      for (int iat=0; iat<nat; iat++)
      {
        dlogPsi_fixed(iw, iat)  = w.G[iat] - dlogPsi_opt(iw,iat);
        d2logPsi_fixed(iw, iat) = w.L[iat] - d2logPsi_opt(iw,iat);
      }
    w.Weight = 1.0;
    // DEBUG
    // W.R = w.R;
    // W.update();
    // Psi.evaluateDeltaLog (W);
    // for (int iat=0; iat<nat; iat++) {
    // 	fprintf (stderr, "CPU: Grad=[%12.6e, %12.6e, %12.6e]  Lap=%12.6e\n",
    // 		 W.G[iat][0],W.G[iat][1],W.G[iat][2],W.L[iat]);
    // 	fprintf (stderr, "GPU: Grad=[%12.6e, %12.6e, %12.6e]  Lap=%12.6e\n",
    // 		 dlogPsi_opt(iw, iat)[0], dlogPsi_opt(iw, iat)[1],  dlogPsi_opt(iw, iat)[2],
    // 		 d2logPsi_opt(iw, iat));
    // 	// dlogPsi_fixed(iw, iat)  = (*dLogPsi[iw])[iat];
    // 	// d2logPsi_fixed(iw, iat) = (*d2LogPsi[iw])[iat];
    // }
  }
  // it=W.begin();
  // for(int iw=0; it!=it_end; ++it,++iw)  {
  //   Walker_t& thisWalker(**it);
  //   W.R=thisWalker.R;
  //   W.update();
  //   Return_t* restrict saved= &(Records(iw,0));
  //   Psi.evaluateDeltaLog(W, saved[LOGPSI_FIXED], saved[LOGPSI_FREE], *dLogPsi[iw],*d2LogPsi[iw]);
  //   // Return_t x= H.evaluate(W);
  //   // e0 += saved[ENERGY_TOT] = x;
  //   // e2 += x*x;
  //   // saved[ENERGY_FIXED] = H.getLocalPotential();
  //   for (int iat=0; iat<nat; iat++) {
  // 	fprintf (stderr, "CPU: Grad=[%12.6e, %12.6e, %12.6e]  Lap=%12.6e\n",
  // 		 (*dLogPsi[iw])[iat][0],  (*dLogPsi[iw])[iat][1], (*dLogPsi[iw])[iat][2],
  // 		 (*d2LogPsi[iw])[iat]);
  // 	fprintf (stderr, "GPU: Grad=[%12.6e, %12.6e, %12.6e]  Lap=%12.6e\n",
  // 		 dlogPsi_fixed(iw, iat)[0], dlogPsi_fixed(iw, iat)[1],  dlogPsi_fixed(iw, iat)[2],
  // 		 d2logPsi_fixed(iw, iat));
  // 	// dlogPsi_fixed(iw, iat)  = (*dLogPsi[iw])[iat];
  // 	// d2logPsi_fixed(iw, iat) = (*d2LogPsi[iw])[iat];
  //   }
  // }
  et_tot+=e0;
  e2_tot+=e2;
  //Need to sum over the processors
  std::vector<Return_t> etemp(3);
  etemp[0]=et_tot;
  etemp[1]=static_cast<Return_t>(W.getActiveWalkers());
  etemp[2]=e2_tot;
  myComm->allreduce(etemp);
  Etarget = static_cast<Return_t>(etemp[0]/etemp[1]);
  NumSamples = static_cast<int>(etemp[1]);
  app_log() << "  VMC Eavg = " << Etarget << std::endl;
  app_log() << "  VMC Evar = " << etemp[2]/etemp[1]-Etarget*Etarget << std::endl;
  app_log() << "  Total weights = " << etemp[1] << std::endl;
  app_log().flush();
  setTargetEnergy(Etarget);
  ReportCounter=0;

  //collect SumValue for computedCost
  NumWalkersEff = etemp[1];
  SumValue[SUM_WGT] = etemp[1];
  SumValue[SUM_WGTSQ] = etemp[1];
  SumValue[SUM_E_WGT] = etemp[0];
  SumValue[SUM_ESQ_WGT] = etemp[2];
  SumValue[SUM_E_BARE] = etemp[0];
  SumValue[SUM_ESQ_BARE] = etemp[2];
  SumValue[SUM_ABSE_BARE] = 0.0;
}

void QMCCostFunctionCUDA::resetPsi(bool final_reset)
{
  if(OptVariables.size() < OptVariablesForPsi.size())
  {
    for(int i=0; i<equalVarMap.size(); ++i)
      OptVariablesForPsi[equalVarMap[i][0]]=OptVariables[equalVarMap[i][1]];
  }
  else
    for(int i=0; i<OptVariables.size(); ++i)
      OptVariablesForPsi[i]=OptVariables[i];
  //cout << "######### QMCCostFunctionCUDA::resetPsi " << std::endl;
  //OptVariablesForPsi.print(std::cout);
  //cout << "-------------------------------------- " << std::endl;
  Psi.resetParameters(OptVariablesForPsi);
  if (final_reset)
  {
    Psi.stopOptimization();
    for (int i=0; i<psiClones.size(); ++i)
      psiClones[i]->stopOptimization();
  }
}


void
QMCCostFunctionCUDA::GradCost(std::vector<Return_t>& PGradient, const std::vector<Return_t>& PM, Return_t FiniteDiff)
{
  if (std::abs(FiniteDiff) > 1.0e-10)
  {
    QMCTraits::RealType dh=1.0/(2.0*FiniteDiff);
    for (int i=0; i<NumOptimizables ; i++)
    {
      for (int j=0; j<NumOptimizables; j++)
        OptVariables[j]=PM[j];
      OptVariables[i] = PM[i]+ FiniteDiff;
      QMCTraits::RealType CostPlus = this->Cost();
      OptVariables[i] = PM[i]- FiniteDiff;
      QMCTraits::RealType CostMinus = this->Cost();
      PGradient[i]= (CostPlus-CostMinus)*dh;
    }
  }
  else
  {
    for (int j=0; j<NumOptimizables; j++)
      OptVariables[j]=PM[j];
    resetPsi();
    //evaluate new local energies and derivatives
    NumWalkersEff=correlatedSampling(true);
    //Estimators::accumulate has been called by correlatedSampling
    curAvg_w = SumValue[SUM_E_WGT]/SumValue[SUM_WGT];
    Return_t curAvg2_w = curAvg_w*curAvg_w;
    curVar_w = SumValue[SUM_ESQ_WGT]/SumValue[SUM_WGT]-curAvg_w*curAvg_w;
    std::vector<Return_t> EDtotals(NumOptimizables,0.0);
    std::vector<Return_t> EDtotals_w(NumOptimizables,0.0);
    std::vector<Return_t> E2Dtotals_w(NumOptimizables,0.0);
    std::vector<Return_t> URV(NumOptimizables,0.0);
    std::vector<Return_t> HD_avg(NumOptimizables,0.0);
    Return_t wgtinv = 1.0/SumValue[SUM_WGT];
    Return_t delE_bar;
    int nw=W.getActiveWalkers();
    for (int iw=0; iw<nw; iw++)
    {
      const Return_t* restrict saved = &(Records(iw,0));
      Return_t weight=saved[REWEIGHT]*wgtinv;
      Return_t eloc_new=saved[ENERGY_NEW];
      delE_bar += weight*std::pow(std::abs(eloc_new-EtargetEff),PowerE);
      std::vector<Return_t> &Dsaved = TempDerivRecords[iw];
      std::vector<Return_t> &HDsaved= TempHDerivRecords[iw];
      for (int pm=0; pm<NumOptimizables; pm++)
        HD_avg[pm]+= HDsaved[pm];
    }
    myComm->allreduce(HD_avg);
    myComm->allreduce(delE_bar);
    for (int pm=0; pm<NumOptimizables; pm++)
      HD_avg[pm] *= 1.0/static_cast<Return_t>(NumSamples);
    for (int iw=0; iw<nw; iw++)
    {
      const Return_t* restrict saved = &(Records(iw,0));
      Return_t weight=saved[REWEIGHT]*wgtinv;
      Return_t eloc_new=saved[ENERGY_NEW];
      Return_t delta_l = (eloc_new-curAvg_w);
      bool ltz(true);
      if (eloc_new-EtargetEff<0)
        ltz=false;
      Return_t delE=std::pow(std::abs(eloc_new-EtargetEff),PowerE);
      Return_t ddelE = PowerE*std::pow(std::abs(eloc_new-EtargetEff),PowerE-1);
      std::vector<Return_t> &Dsaved=  TempDerivRecords[iw];
      std::vector<Return_t> &HDsaved= TempHDerivRecords[iw];
      for (int pm=0; pm<NumOptimizables; pm++)
      {
        EDtotals_w[pm] += weight*(HDsaved[pm] + 2.0*Dsaved[pm]*delta_l );
        URV[pm] += 2.0*(eloc_new*HDsaved[pm] - curAvg*HD_avg[pm]);
        if (ltz)
          EDtotals[pm]+= weight*( 2.0*Dsaved[pm]*(delE-delE_bar) + ddelE*HDsaved[pm]);
        else
          EDtotals[pm] += weight*( 2.0*Dsaved[pm]*(delE-delE_bar) - ddelE*HDsaved[pm]);
      }
    }
    myComm->allreduce(EDtotals);
    myComm->allreduce(EDtotals_w);
    myComm->allreduce(URV);
    Return_t smpinv=1.0/static_cast<Return_t>(NumSamples);
    for (int iw=0; iw<nw; iw++)
    {
      const Return_t* restrict saved = &(Records(iw,0));
      Return_t weight=saved[REWEIGHT]*wgtinv;
      Return_t eloc_new=saved[ENERGY_NEW];
      Return_t delta_l = (eloc_new-curAvg_w);
      Return_t sigma_l = delta_l*delta_l;
      std::vector<Return_t> &Dsaved=  TempDerivRecords[iw];
      std::vector<Return_t> &HDsaved= TempHDerivRecords[iw];
      for (int pm=0; pm<NumOptimizables; pm++)
        E2Dtotals_w[pm] += weight*2.0*( Dsaved[pm]*(sigma_l-curVar_w) + delta_l*(HDsaved[pm]-EDtotals_w[pm]) );
    }
    myComm->allreduce(E2Dtotals_w);
    for (int pm=0; pm<NumOptimizables; pm++)
      URV[pm] *=smpinv;
    for (int j=0; j<NumOptimizables; j++)
    {
      PGradient[j] = 0.0;
      if (std::abs(w_var) > 1.0e-10)
        PGradient[j] += w_var*E2Dtotals_w[j];
      if (std::abs(w_en)  > 1.0e-10)
        PGradient[j] += w_en*EDtotals_w[j];
      if (std::abs(w_w)   > 1.0e-10)
        PGradient[j] += w_w*URV[j];
      if (std::abs(w_abs) > 1.0e-10)
        PGradient[j] += w_abs*EDtotals[j];
    }
    IsValid=true;
    //         if ((CSWeight/wgtinv) < MinNumWalkers)
    if (NumWalkersEff < MinNumWalkers*NumSamples)
    {
      ERRORMSG("CostFunction-> Number of Effective Walkers is too small " << NumWalkersEff<< "Minimum required"<<MinNumWalkers*NumSamples)
      //ERRORMSG("Going to stop now.")
      IsValid=false;
    }
  }
}


QMCCostFunctionCUDA::Return_t QMCCostFunctionCUDA::fillOverlapHamiltonianMatrices(Matrix<Return_t>& Left, Matrix<Return_t>& Right)
{
  RealType b1,b2;
  if (GEVType=="H2")
  {
    b1=w_beta;
    b2=0;
  }
  else
  {
    b2=w_beta;
    b1=0;
  }
  Right=0.0;
  Left=0.0;
  //Is this really needed here?
  //resetPsi();
  //Return_t NWE = NumWalkersEff=correlatedSampling(true);
  curAvg_w = SumValue[SUM_E_WGT]/SumValue[SUM_WGT];
  Return_t curAvg2_w = SumValue[SUM_ESQ_WGT]/SumValue[SUM_WGT];
  RealType H2_avg = 1.0/(curAvg_w*curAvg_w);
  RealType V_avg = curAvg2_w - curAvg_w*curAvg_w;
  std::vector<Return_t> D_avg(NumParams(),0);
  Return_t wgtinv = 1.0/SumValue[SUM_WGT];
  int nw=W.getActiveWalkers();
  for (int iw=0; iw<nw; iw++)
  {
    const Return_t* restrict saved= &(Records(iw,0));
    Return_t weight=saved[REWEIGHT]*wgtinv;
    const std::vector<Return_t> &Dsaved = TempDerivRecords[iw];
    for (int pm=0; pm<NumParams(); pm++)
    {
      D_avg[pm]+= Dsaved[pm]*weight;
    }
  }
  myComm->allreduce(D_avg);
  for (int iw=0; iw<nw; iw++)
  {
    const Return_t* restrict saved = &Records(iw,0);
    Return_t weight=saved[REWEIGHT]*wgtinv;
    Return_t eloc_new=saved[ENERGY_NEW];
    const std::vector<Return_t> &Dsaved = TempDerivRecords[iw];
    const std::vector<Return_t> &HDsaved= TempHDerivRecords[iw];
    /*
    for (int pm=0; pm<NumParams(); pm++)
    {
      Return_t wfe = (HDsaved[pm] + Dsaved[pm]*(eloc_new - curAvg_w) )*weight;
      Return_t wfm = (HDsaved[pm] - 2.0*Dsaved[pm]*(eloc_new - curAvg_w) )*weight;
      Return_t wfd = (Dsaved[pm]-D_avg[pm])*weight;
      //                 H2
      Right(0,pm+1) += b1*H2_avg*wfe*(eloc_new);
      Right(pm+1,0) += b1*H2_avg*wfe*(eloc_new);
      Return_t vterm = HDsaved[pm]*(eloc_new-curAvg_w)+(eloc_new*eloc_new-curAvg2_w)*Dsaved[pm]-2.0*curAvg_w*Dsaved[pm]*(eloc_new - curAvg_w);
      //                 Variance
      Left(0,pm+1) += b2*vterm*weight;
      Left(pm+1,0) += b2*vterm*weight;
      //                 Hamiltonian
      Left(0,pm+1) += (1-b2)*wfe;
      Left(pm+1,0) += (1-b2)*wfd*(eloc_new-curAvg_w);
      for (int pm2=0; pm2<NumParams(); pm2++)
      {
        //                   Hamiltonian
        Left(pm+1,pm2+1) += (1-b2)*wfd*(HDsaved[pm2]+ Dsaved[pm2]*(eloc_new-curAvg_w));
        //                   Overlap
        RealType ovlij=wfd*(Dsaved[pm2]-D_avg[pm2]);
        Right(pm+1,pm2+1) += ovlij;
        //                   Variance
        RealType varij=wfm*(HDsaved[pm2] - 2.0*Dsaved[pm2]*(eloc_new - curAvg_w));
        Left(pm+1,pm2+1) +=  b2*(varij+V_avg*ovlij);
        //                   H2
        Right(pm+1,pm2+1) += b1*H2_avg*varij; //(b1*H2_avg)*wfe*(HDsaved[pm2]+ Dsaved[pm2]*(eloc_new - curAvg_w));
      }
    } */
    #pragma omp parallel for
    for (int pm=0; pm<NumParams(); pm++)
    {
      Return_t wfe = (HDsaved[pm] +(Dsaved[pm]-D_avg[pm])*eloc_new)*weight;
      Return_t wfd = (Dsaved[pm]-D_avg[pm])*weight;
      Return_t vterm =  HDsaved[pm]*(eloc_new-curAvg_w)+(Dsaved[pm]-D_avg[pm])*eloc_new*(eloc_new-2.0*curAvg_w);
      //                Return_t vterm = (HDsaved[pm]+(Dsaved[pm]-D_avg[pm])*eloc_new -curAvg_w)*(eloc_new-curAvg_w);
      //                 H2
      Right(0,pm+1) += b1*H2_avg*vterm*weight;
      Right(pm+1,0) += b1*H2_avg*vterm*weight;
      //                 Variance
      Left(0,pm+1) += b2*vterm*weight;
      Left(pm+1,0) += b2*vterm*weight;
      //                 Hamiltonian
      Left(0,pm+1) += (1-b2)*wfe;
      Left(pm+1,0) += (1-b2)*wfd*eloc_new;
      for (int pm2=0; pm2<NumParams(); pm2++)
      {
        //                Hamiltonian
        Left(pm+1,pm2+1) += (1-b2)*wfd*(HDsaved[pm2] + (Dsaved[pm2]-D_avg[pm2])*eloc_new);
        //                Overlap
        RealType ovlij=wfd*(Dsaved[pm2]-D_avg[pm2]);
        Right(pm+1,pm2+1) += ovlij;
        //                Variance
        RealType varij=weight*(HDsaved[pm] - 2.0*(Dsaved[pm]-D_avg[pm])*eloc_new)*(HDsaved[pm2] - 2.0*(Dsaved[pm2]-D_avg[pm2])*eloc_new);
        //                  RealType varij=weight*(HDsaved[pm] +(Dsaved[pm]-D_avg[pm])*eloc_new-curAvg_w)*
        //                                      (HDsaved[pm2] + (Dsaved[pm2]-D_avg[pm2])*eloc_new-curAvg_w);
        Left(pm+1,pm2+1) +=  b2*(varij+V_avg*ovlij);
//                H2
        Right(pm+1,pm2+1) += b1*H2_avg*varij;
      }
    }
  }
  myComm->allreduce(Right);
  myComm->allreduce(Left);
  Left(0,0) = (1-b2)*curAvg_w + b2*V_avg;
  Right(0,0) = 1.0+b1*H2_avg*V_avg;
  if (GEVType=="H2")
    return H2_avg;
  return 1.0;
}

}
