//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/QMCCostFunctionOMP.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Particle/HDFWalkerInputCollect.h"
#include "Message/CommOperators.h"
//#define QMCCOSTFUNCTION_DEBUG


namespace qmcplusplus
{

QMCCostFunctionOMP::QMCCostFunctionOMP(MCWalkerConfiguration& w,
                                       TrialWaveFunction& psi, QMCHamiltonian& h, HamiltonianPool& hpool):
  QMCCostFunctionBase(w,psi,h), CloneManager(hpool)
{
  CSWeight=1.0;
  app_log()<<" Using QMCCostFunctionOMP::QMCCostFunctionOMP"<< std::endl;
}


/** Clean up the vector */
QMCCostFunctionOMP::~QMCCostFunctionOMP()
{
  delete_iter(H_KE_Node.begin(),H_KE_Node.end());
  delete_iter(RngSaved.begin(),RngSaved.end());
  delete_iter(RecordsOnNode.begin(),RecordsOnNode.end());
  delete_iter(DerivRecords.begin(),DerivRecords.end());
  delete_iter(HDerivRecords.begin(),HDerivRecords.end());
}

void QMCCostFunctionOMP::GradCost(std::vector<Return_t>& PGradient, const std::vector<Return_t>& PM, Return_t FiniteDiff)
{
  if (FiniteDiff > 0)
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
    Return_t delE_bar=0;
    for (int ip=0; ip<NumThreads; ip++)
    {
      int nw=wClones[ip]->numSamples();
      for (int iw=0; iw<nw; iw++)
      {
        const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
        Return_t weight=saved[REWEIGHT]*wgtinv;
        Return_t eloc_new=saved[ENERGY_NEW];
        delE_bar += weight*std::pow(std::abs(eloc_new-EtargetEff),PowerE);
        const Return_t* HDsaved= (*HDerivRecords[ip])[iw];
        for (int pm=0; pm<NumOptimizables; pm++)
          HD_avg[pm]+= HDsaved[pm];
      }
    }
    myComm->allreduce(HD_avg);
    myComm->allreduce(delE_bar);
    for (int pm=0; pm<NumOptimizables; pm++)
      HD_avg[pm] *= 1.0/static_cast<Return_t>(NumSamples);
    for (int ip=0; ip<NumThreads; ip++)
    {
      int nw=wClones[ip]->numSamples();
      for (int iw=0; iw<nw; iw++)
      {
        const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
        Return_t weight=saved[REWEIGHT]*wgtinv;
        Return_t eloc_new=saved[ENERGY_NEW];
        Return_t delta_l = (eloc_new-curAvg_w);
        bool ltz(true);
        if (eloc_new-EtargetEff<0)
          ltz=false;
        Return_t delE=std::pow(std::abs(eloc_new-EtargetEff),PowerE);
        Return_t ddelE = PowerE*std::pow(std::abs(eloc_new-EtargetEff),PowerE-1);
        const Return_t* Dsaved= (*DerivRecords[ip])[iw];
        const Return_t* HDsaved= (*HDerivRecords[ip])[iw];
        for (int pm=0; pm<NumOptimizables; pm++)
        {
          EDtotals_w[pm] += weight*(HDsaved[pm] + 2.0*Dsaved[pm]*delta_l);
          URV[pm] += 2.0*(eloc_new*HDsaved[pm] - curAvg*HD_avg[pm]);
          if (ltz)
            EDtotals[pm]+= weight*(2.0*Dsaved[pm]*(delE-delE_bar) + ddelE*HDsaved[pm]);
          else
            EDtotals[pm] += weight*(2.0*Dsaved[pm]*(delE-delE_bar) - ddelE*HDsaved[pm]);
        }
      }
    }
    myComm->allreduce(EDtotals);
    myComm->allreduce(EDtotals_w);
    myComm->allreduce(URV);
    Return_t smpinv=1.0/static_cast<Return_t>(NumSamples);
    for (int ip=0; ip<NumThreads; ip++)
    {
      int nw=wClones[ip]->numSamples();
      for (int iw=0; iw<nw; iw++)
      {
        const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
        Return_t weight=saved[REWEIGHT]*wgtinv;
        Return_t eloc_new=saved[ENERGY_NEW];
        Return_t delta_l = (eloc_new-curAvg_w);
        Return_t sigma_l = delta_l*delta_l;
        const Return_t* Dsaved= (*DerivRecords[ip])[iw];
        const Return_t* HDsaved= (*HDerivRecords[ip])[iw];
        for (int pm=0; pm<NumOptimizables; pm++)
        {
          E2Dtotals_w[pm] += weight*2.0*(Dsaved[pm]*(sigma_l-curVar_w) + delta_l*(HDsaved[pm]-EDtotals_w[pm]));
        }
      }
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


void QMCCostFunctionOMP::getConfigurations(const std::string& aroot)
{
  //makeClones(W,Psi,H);
  if (H_KE_Node.empty())
  {
    app_log() << "  QMCCostFunctionOMP is created with " << NumThreads << " threads." << std::endl;
    //make H_KE_Node
    H_KE_Node.resize(NumThreads,0);
    RecordsOnNode.resize(NumThreads,0);
    DerivRecords.resize(NumThreads,0);
    HDerivRecords.resize(NumThreads,0);
  }
  
  app_log()<<"  Using Nonlocal PP in Opt: " << includeNonlocalH << std::endl;
  outputManager.pause();
  //#pragma omp parallel for
  for (int ip=0; ip<NumThreads; ++ip)
  {
    if (H_KE_Node[ip]==0)
    {
      H_KE_Node[ip]= new QMCHamiltonian;
      H_KE_Node[ip]->addOperator(hClones[ip]->getHamiltonian("Kinetic"),"Kinetic");
      if (includeNonlocalH!="no")
      {
        QMCHamiltonianBase* a=hClones[ip]->getHamiltonian(includeNonlocalH);
        if(a)
        {
          app_log()<<" Found non-local Hamiltonian element named "<<includeNonlocalH<< std::endl;
          H_KE_Node[ip]->addOperator(a,includeNonlocalH);
        }
        else
          app_log()<<" Did not find non-local Hamiltonian element named "<<includeNonlocalH<< std::endl;
      }
    }
  }

  //load samples from SampleStack
  outputManager.resume();
  app_log() << "   Number of samples loaded to each thread : ";
  wPerNode[0]=0;
  for (int ip=0; ip<NumThreads; ++ip)
  {
    wPerNode[ip+1] = wPerNode[ip] + wClones[ip]->numSamples();
    app_log() << wClones[ip]->numSamples() <<  " " ;
  }
  app_log() << std::endl;
  app_log().flush();

  if (dLogPsi.size() != wPerNode[NumThreads])
  {
    delete_iter(dLogPsi.begin(),dLogPsi.end());
    delete_iter(d2LogPsi.begin(),d2LogPsi.end());
    int nptcl=W.getTotalNum();
    int nwtot=wPerNode[NumThreads];
    dLogPsi.resize(nwtot);
    d2LogPsi.resize(nwtot);
    for (int i=0; i<nwtot; ++i)
      dLogPsi[i]=new ParticleGradient_t(nptcl);
    for (int i=0; i<nwtot; ++i)
      d2LogPsi[i]=new ParticleLaplacian_t(nptcl);
  }
}

/** evaluate everything before optimization */
void QMCCostFunctionOMP::checkConfigurations()
{
  RealType et_tot=0.0;
  RealType e2_tot=0.0;
#pragma omp parallel reduction(+:et_tot,e2_tot)
  {
    int ip = omp_get_thread_num();
    MCWalkerConfiguration& wRef(*wClones[ip]);
    if (RecordsOnNode[ip] ==0)
    {
      RecordsOnNode[ip]=new Matrix<Return_t>;
      RecordsOnNode[ip]->resize(wRef.numSamples(),SUM_INDEX_SIZE);
      if (needGrads)
      {
        DerivRecords[ip]=new Matrix<Return_t>;
        DerivRecords[ip]->resize(wRef.numSamples(),NumOptimizables);
        HDerivRecords[ip]=new Matrix<Return_t>;
        HDerivRecords[ip]->resize(wRef.numSamples(),NumOptimizables);
      }
    }
    else if (RecordsOnNode[ip]->size1()!=wRef.numSamples())
    {
      RecordsOnNode[ip]->resize(wRef.numSamples(),SUM_INDEX_SIZE);
      if (needGrads)
      {
        DerivRecords[ip]->resize(wRef.numSamples(),NumOptimizables);
        HDerivRecords[ip]->resize(wRef.numSamples(),NumOptimizables);
      }
    }
    QMCHamiltonianBase* nlpp = (includeNonlocalH =="no")?  0: hClones[ip]->getHamiltonian(includeNonlocalH);
    bool compute_nlpp=useNLPPDeriv && nlpp;
    //set the optimization mode for the trial wavefunction
    psiClones[ip]->startOptimization();
    //    synchronize the random number generator with the node
    (*MoverRng[ip]) = (*RngSaved[ip]);
    hClones[ip]->setRandomGenerator(MoverRng[ip]);
    //int nat = wRef.getTotalNum();
    Return_t e0=0.0;
    //       Return_t ef=0.0;
    Return_t e2=0.0;
    for (int iw=0, iwg=wPerNode[ip]; iw<wRef.numSamples(); ++iw,++iwg)
    {
      wRef.loadSample(wRef.R, iw);
      wRef.update(true);
      wRef.donePbyP();
      Return_t* restrict saved=(*RecordsOnNode[ip])[iw];
      psiClones[ip]->evaluateDeltaLog(wRef, saved[LOGPSI_FIXED], saved[LOGPSI_FREE], *dLogPsi[iwg], *d2LogPsi[iwg]);
      saved[REWEIGHT]=1.0;
      Return_t etmp;
      if (needGrads)
      {
        //allocate vector
        std::vector<Return_t> Dsaved(NumOptimizables,0.0);
        std::vector<Return_t> HDsaved(NumOptimizables,0.0);
        psiClones[ip]->evaluateDerivatives(wRef, OptVariablesForPsi, Dsaved, HDsaved);
        etmp =hClones[ip]->evaluateValueAndDerivatives(wRef,OptVariablesForPsi,Dsaved,HDsaved,compute_nlpp);
        copy(Dsaved.begin(),Dsaved.end(),(*DerivRecords[ip])[iw]);
        copy(HDsaved.begin(),HDsaved.end(),(*HDerivRecords[ip])[iw]);
      }
      else
        etmp= hClones[ip]->evaluate(wRef);

      e0 += saved[ENERGY_TOT] = saved[ENERGY_NEW] = etmp;
      e2 += etmp*etmp;
      saved[ENERGY_FIXED] = hClones[ip]->getLocalPotential();
      if(nlpp)
        saved[ENERGY_FIXED] -= nlpp->Value;
    }
    //add them all using reduction
    et_tot+=e0;
    e2_tot+=e2;
    // #pragma omp atomic
    //       eft_tot+=ef;
  }
  OptVariablesForPsi.setComputed();
  //     app_log() << "  VMC Efavg = " << eft_tot/static_cast<Return_t>(wPerNode[NumThreads]) << std::endl;
  //Need to sum over the processors
  std::vector<Return_t> etemp(3);
  etemp[0]=et_tot;
  etemp[1]=static_cast<Return_t>(wPerNode[NumThreads]);
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

#ifdef HAVE_LMY_ENGINE
/** evaluate everything before optimization */
void QMCCostFunctionOMP::engine_checkConfigurations(cqmc::engine::LMYEngine * EngineObj)
{
  RealType et_tot=0.0;
  RealType e2_tot=0.0;
#pragma omp parallel reduction(+:et_tot,e2_tot)
  {
    int ip = omp_get_thread_num();
    MCWalkerConfiguration& wRef(*wClones[ip]);
    if (RecordsOnNode[ip] ==0)
    {
      RecordsOnNode[ip]=new Matrix<Return_t>;
      RecordsOnNode[ip]->resize(wRef.numSamples(),SUM_INDEX_SIZE);
      if (needGrads)
      {
        DerivRecords[ip]=new Matrix<Return_t>;
        //DerivRecords[ip]->resize(wRef.numSamples(),NumOptimizables);
        HDerivRecords[ip]=new Matrix<Return_t>;
        //HDerivRecords[ip]->resize(wRef.numSamples(),NumOptimizables);
      }
    }
    else if (RecordsOnNode[ip]->size1()!=wRef.numSamples())
    {
      RecordsOnNode[ip]->resize(wRef.numSamples(),SUM_INDEX_SIZE);
      if (needGrads)
      {
        //DerivRecords[ip]->resize(wRef.numSamples(),NumOptimizables);
        //HDerivRecords[ip]->resize(wRef.numSamples(),NumOptimizables);
      }
    }
    QMCHamiltonianBase* nlpp = (includeNonlocalH =="no")?  0: hClones[ip]->getHamiltonian(includeNonlocalH.c_str());
    bool compute_nlpp=useNLPPDeriv && nlpp;
    //set the optimization mode for the trial wavefunction
    psiClones[ip]->startOptimization();
    //    synchronize the random number generator with the node
    (*MoverRng[ip]) = (*RngSaved[ip]);
    hClones[ip]->setRandomGenerator(MoverRng[ip]);
    //int nat = wRef.getTotalNum();
    Return_t e0=0.0;
    //       Return_t ef=0.0;
    Return_t e2=0.0;
    for (int iw=0, iwg=wPerNode[ip]; iw<wRef.numSamples(); ++iw,++iwg)
    {
      wRef.loadSample(wRef.R, iw);
      wRef.update(true);
      wRef.donePbyP();
      Return_t* restrict saved=(*RecordsOnNode[ip])[iw];
      psiClones[ip]->evaluateDeltaLog(wRef, saved[LOGPSI_FIXED], saved[LOGPSI_FREE], *dLogPsi[iwg], *d2LogPsi[iwg]);
      saved[REWEIGHT]=1.0;
      Return_t etmp;
      if (needGrads)
      {
        //allocate vector
        std::vector<Return_t> Dsaved(NumOptimizables,0.0);
        std::vector<Return_t> HDsaved(NumOptimizables,0.0);
        psiClones[ip]->evaluateDerivatives(wRef, OptVariablesForPsi, Dsaved, HDsaved);
        etmp =hClones[ip]->evaluateValueAndDerivatives(wRef,OptVariablesForPsi,Dsaved,HDsaved,compute_nlpp);
        //std::copy(Dsaved.begin(),Dsaved.end(),(*DerivRecords[ip])[iw]);
        //std::copy(HDsaved.begin(),HDsaved.end(),(*HDerivRecords[ip])[iw]);

        // add non-differentiated derivative vector
        std::vector<Return_t> der_rat_samp(NumOptimizables+1, 0.0);
        std::vector<Return_t> le_der_samp(NumOptimizables+1, 0.0);

        // dervative vectors
        der_rat_samp.at(0) = 1.0;
        for (int i = 0; i < Dsaved.size(); i++) 
          der_rat_samp.at(i+1) = Dsaved.at(i);

        // evaluate local energy
        etmp= hClones[ip]->evaluate(wRef);
        
        // energy dervivatives 
        le_der_samp.at(0) = etmp;
        for (int i = 0; i < HDsaved.size(); i++) 
          le_der_samp.at(i+1) = HDsaved.at(i) + etmp * Dsaved.at(i);
        
#ifdef HAVE_LMY_ENGINE
        // pass into engine
        EngineObj->take_sample(der_rat_samp, le_der_samp, le_der_samp, 1.0, saved[REWEIGHT]);
#endif

        //etmp= hClones[ip]->evaluate(wRef);
      }
      else
        etmp= hClones[ip]->evaluate(wRef);

      e0 += saved[ENERGY_TOT] = etmp;
      e2 += etmp*etmp;
      saved[ENERGY_FIXED] = hClones[ip]->getLocalPotential();
      if(nlpp)
        saved[ENERGY_FIXED] -= nlpp->Value;
    }

    //add them all using reduction
    et_tot+=e0;
    e2_tot+=e2;
    // #pragma omp atomic
    //       eft_tot+=ef;
  }
#ifdef HAVE_LMY_ENGINE
  // engine finish taking samples 
  EngineObj->sample_finish();
#endif

  if ( EngineObj->block_first() ) {
    OptVariablesForPsi.setComputed();
    app_log() << "calling setComputed function" << std::endl;
  }
  //     app_log() << "  VMC Efavg = " << eft_tot/static_cast<Return_t>(wPerNode[NumThreads]) << endl;
  //Need to sum over the processors
  std::vector<Return_t> etemp(3);
  etemp[0]=et_tot;
  etemp[1]=static_cast<Return_t>(wPerNode[NumThreads]);
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
}
#endif

void QMCCostFunctionOMP::resetPsi(bool final_reset)
{
  if (OptVariables.size() < OptVariablesForPsi.size())
    for (int i=0; i<equalVarMap.size(); ++i)
      OptVariablesForPsi[equalVarMap[i][0]]=OptVariables[equalVarMap[i][1]];
  else
    for (int i=0; i<OptVariables.size(); ++i)
      OptVariablesForPsi[i]=OptVariables[i];
  if (final_reset)
  {
    #pragma omp parallel for
    for (int i=0; i<psiClones.size(); ++i)
      psiClones[i]->stopOptimization();
  }
  //cout << "######### QMCCostFunctionOMP::resetPsi " << std::endl;
  //OptVariablesForPsi.print(std::cout);
  //cout << "-------------------------------------- " << std::endl;
  Psi.resetParameters(OptVariablesForPsi);
  for (int i=0; i<psiClones.size(); ++i)
    psiClones[i]->resetParameters(OptVariablesForPsi);
}

QMCCostFunctionOMP::Return_t QMCCostFunctionOMP::correlatedSampling(bool needGrad)
{
  for(int ip=0; ip<NumThreads; ++ip)
  {
    //    synchronize the random number generator with the node
    (*MoverRng[ip]) = (*RngSaved[ip]);
    hClones[ip]->setRandomGenerator(MoverRng[ip]);
  }

  const bool nlpp = (includeNonlocalH != "no");
  Return_t wgt_tot=0.0;
  Return_t wgt_tot2=0.0;
  Return_t NSm1 = 1.0/NumSamples;
  Return_t inv_n_samples=1.0/NumSamples;
#pragma omp parallel reduction(+:wgt_tot,wgt_tot2)
  {
    int ip = omp_get_thread_num();
    bool compute_nlpp=useNLPPDeriv && (includeNonlocalH != "no");
    bool compute_all_from_scratch=(includeNonlocalH != "no"); //true if we have nlpp

    MCWalkerConfiguration& wRef(*wClones[ip]);
    Return_t wgt_node=0.0, wgt_node2=0.0;
    for (int iw=0, iwg=wPerNode[ip]; iw<wRef.numSamples(); ++iw,++iwg)
    {
      wRef.loadSample(wRef.R, iw);
      wRef.update(true);
      if(nlpp) wRef.donePbyP(true);
      Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
      Return_t logpsi;
      logpsi=psiClones[ip]->evaluateDeltaLog(wRef,compute_all_from_scratch);
      wRef.G += *dLogPsi[iwg];
      wRef.L += *d2LogPsi[iwg];
      Return_t weight = saved[REWEIGHT] = vmc_or_dmc*(logpsi-saved[LOGPSI_FREE]);
      if (needGrad)
      {
        std::vector<Return_t> Dsaved(NumOptimizables,0);
        std::vector<Return_t> HDsaved(NumOptimizables,0);
        psiClones[ip]->evaluateDerivatives(wRef, OptVariablesForPsi, Dsaved, HDsaved);
        saved[ENERGY_NEW] =
          H_KE_Node[ip]->evaluateValueAndDerivatives(wRef,OptVariablesForPsi,Dsaved,HDsaved,compute_nlpp)
          +saved[ENERGY_FIXED];;
        for( int i=0; i<NumOptimizables; i++)
          if(OptVariablesForPsi.recompute(i))
          {
            (*DerivRecords[ip])(iw,i) = Dsaved[i];
            (*HDerivRecords[ip])(iw,i) = HDsaved[i];
          }
        //saved[ENERGY_NEW] = H_KE_Node[ip]->evaluate(wRef) + saved[ENERGY_FIXED];
      }
      else
        saved[ENERGY_NEW] = H_KE_Node[ip]->evaluate(wRef) + saved[ENERGY_FIXED];
      wgt_node+=inv_n_samples*weight;
      wgt_node2+=inv_n_samples*weight*weight;
    }
    wgt_tot += wgt_node;
    wgt_tot2 += wgt_node2;
  }
  //this is MPI barrier
  OHMMS::Controller->barrier();
  //collect the total weight for normalization and apply maximum weight
  myComm->allreduce(wgt_tot);
  myComm->allreduce(wgt_tot2);
  //    app_log()<<"Before Purge"<<wgt_tot<<" "<<wgt_tot2<< std::endl;
  Return_t wgtnorm = (wgt_tot==0)?0:wgt_tot;
  wgt_tot=0.0;
  for (int ip=0; ip<NumThreads; ip++)
  {
    int nw=wClones[ip]->numSamples();
    for (int iw=0; iw<nw; iw++)
    {
      Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
      saved[REWEIGHT] = std::min(std::exp(saved[REWEIGHT]-wgtnorm), std::numeric_limits<Return_t>::max()*(RealType)0.1 );
      wgt_tot+= inv_n_samples*saved[REWEIGHT];
    }
  }
  myComm->allreduce(wgt_tot);
  //    app_log()<<"During Purge"<<wgt_tot<<" "<< std::endl;
  wgtnorm =(wgt_tot==0)?1:1.0/wgt_tot;
  wgt_tot=0.0;
  for (int ip=0; ip<NumThreads; ip++)
  {
    int nw=wClones[ip]->numSamples();
    for (int iw=0; iw<nw; iw++)
    {
      Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
      saved[REWEIGHT] = std::min(saved[REWEIGHT]*wgtnorm,MaxWeight);
      wgt_tot+= inv_n_samples*saved[REWEIGHT];
    }
  }
  myComm->allreduce(wgt_tot);
  //    app_log()<<"After Purge"<<wgt_tot<<" "<< std::endl;
  for (int i=0; i<SumValue.size(); i++)
    SumValue[i]=0.0;
  CSWeight=wgt_tot=(wgt_tot==0)?1:1.0/wgt_tot;
  for (int ip=0; ip<NumThreads; ip++)
  {
    int nw=wClones[ip]->numSamples();
    for (int iw=0; iw<nw; iw++)
    {
      const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
      Return_t weight=saved[REWEIGHT]*wgt_tot;
      Return_t eloc_new=saved[ENERGY_NEW];
      Return_t delE=std::pow(std::abs(eloc_new-EtargetEff),PowerE);
      SumValue[SUM_E_BARE] += eloc_new;
      SumValue[SUM_ESQ_BARE] += eloc_new*eloc_new;
      SumValue[SUM_ABSE_BARE] += delE;
      SumValue[SUM_E_WGT] += eloc_new*saved[REWEIGHT];
      SumValue[SUM_ESQ_WGT] += eloc_new*eloc_new*saved[REWEIGHT];
      SumValue[SUM_ABSE_WGT] += delE*saved[REWEIGHT];
      SumValue[SUM_WGT] += saved[REWEIGHT];
      SumValue[SUM_WGTSQ] += saved[REWEIGHT]*saved[REWEIGHT];
    }
  }
  //collect everything
  myComm->allreduce(SumValue);
  //     for (int i=0; i<SumValue.size(); i++) std::cerr <<SumValue[i]<<"  ";
  //     std::cerr << std::endl;
  //     app_log()<<"After purge Energy Variance Weight "
  //      << SumValue[SUM_E_WGT]/SumValue[SUM_WGT] << " "
  //      << SumValue[SUM_ESQ_WGT]/SumValue[SUM_WGT] -(SumValue[SUM_E_WGT]/SumValue[SUM_WGT])*(SumValue[SUM_E_WGT]/SumValue[SUM_WGT]) << " "
  //      << SumValue[SUM_WGT]*SumValue[SUM_WGT]/SumValue[SUM_WGTSQ] << std::endl;
  return SumValue[SUM_WGT]*SumValue[SUM_WGT]/SumValue[SUM_WGTSQ];
}


QMCCostFunctionOMP::Return_t
QMCCostFunctionOMP::fillOverlapHamiltonianMatrices(Matrix<Return_t>& Left, Matrix<Return_t>& Right)
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

  //     resetPsi();
  //     Return_t NWE = NumWalkersEff=correlatedSampling(true);
  curAvg_w = SumValue[SUM_E_WGT]/SumValue[SUM_WGT];
  Return_t curAvg2_w = SumValue[SUM_ESQ_WGT]/SumValue[SUM_WGT];
  //    RealType H2_avg = 1.0/curAvg2_w;
  RealType H2_avg = 1.0/(curAvg_w*curAvg_w);
  //    RealType H2_avg = 1.0/std::sqrt(curAvg_w*curAvg_w*curAvg2_w);
  RealType V_avg = curAvg2_w - curAvg_w*curAvg_w;
  std::vector<Return_t> D_avg(NumParams(),0.0);
  Return_t wgtinv = 1.0/SumValue[SUM_WGT];
  for (int ip=0; ip<NumThreads; ip++)
  {
    int nw=wClones[ip]->numSamples();
    for (int iw=0; iw<nw; iw++)
    {
      const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
      Return_t weight=saved[REWEIGHT]*wgtinv;
      const Return_t* Dsaved= (*DerivRecords[ip])[iw];
      for (int pm=0; pm<NumParams(); pm++)
      {
        D_avg[pm]+= Dsaved[pm]*weight;
      }
    }
  }

  myComm->allreduce(D_avg);

  for (int ip=0; ip<NumThreads; ip++)
  {
    int nw=wClones[ip]->numSamples();
    for (int iw=0; iw<nw; iw++)
    {
      const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
      Return_t weight=saved[REWEIGHT]*wgtinv;
      Return_t eloc_new=saved[ENERGY_NEW];
      const Return_t* Dsaved= (*DerivRecords[ip])[iw];
      const Return_t* HDsaved= (*HDerivRecords[ip])[iw];
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
