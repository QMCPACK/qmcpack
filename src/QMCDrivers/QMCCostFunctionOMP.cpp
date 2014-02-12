//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
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
  app_log()<<" Using QMCCostFunctionOMP::QMCCostFunctionOMP"<<endl;
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


void QMCCostFunctionOMP::resetWalkers()
{
  //remove walkers of the clones
  for (int ip=0; ip<NumThreads; ++ip)
  {
    if(wClones[ip]->getActiveWalkers()>nVMCWalkers[ip])
      wClones[ip]->destroyWalkers(wClones[ip]->getActiveWalkers()-nVMCWalkers[ip]);
    else
      wClones[ip]->createWalkers(nVMCWalkers[ip]-wClones[ip]->getActiveWalkers());
    nVMCWalkers[ip]=0;
  }
}

void QMCCostFunctionOMP::GradCost(vector<Return_t>& PGradient, const vector<Return_t>& PM, Return_t FiniteDiff)
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
    vector<Return_t> EDtotals(NumOptimizables,0.0);
    vector<Return_t> EDtotals_w(NumOptimizables,0.0);
    vector<Return_t> E2Dtotals_w(NumOptimizables,0.0);
    vector<Return_t> URV(NumOptimizables,0.0);
    vector<Return_t> HD_avg(NumOptimizables,0.0);
    Return_t wgtinv = 1.0/SumValue[SUM_WGT];
    Return_t delE_bar=0;
    for (int ip=0; ip<NumThreads; ip++)
    {
      int nw=wClones[ip]->getActiveWalkers();
      for (int iw=0; iw<nw; iw++)
      {
        const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
        Return_t weight=saved[REWEIGHT]*wgtinv;
        Return_t eloc_new=saved[ENERGY_NEW];
        delE_bar += weight*std::pow(abs(eloc_new-EtargetEff),PowerE);
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
      int nw=wClones[ip]->getActiveWalkers();
      for (int iw=0; iw<nw; iw++)
      {
        const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
        Return_t weight=saved[REWEIGHT]*wgtinv;
        Return_t eloc_new=saved[ENERGY_NEW];
        Return_t delta_l = (eloc_new-curAvg_w);
        bool ltz(true);
        if (eloc_new-EtargetEff<0)
          ltz=false;
        Return_t delE=std::pow(abs(eloc_new-EtargetEff),PowerE);
        Return_t ddelE = PowerE*std::pow(abs(eloc_new-EtargetEff),PowerE-1);
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
      int nw=wClones[ip]->getActiveWalkers();
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
      if (std::fabs(w_var) > 1.0e-10)
        PGradient[j] += w_var*E2Dtotals_w[j];
      if (std::fabs(w_en)  > 1.0e-10)
        PGradient[j] += w_en*EDtotals_w[j];
      if (std::fabs(w_w)   > 1.0e-10)
        PGradient[j] += w_w*URV[j];
      if (std::fabs(w_abs) > 1.0e-10)
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


void QMCCostFunctionOMP::getConfigurations(const string& aroot)
{
  //makeClones(W,Psi,H);
  if (H_KE_Node.empty())
  {
    app_log() << "  QMCCostFunctionOMP is created with " << NumThreads << endl;
    //make H_KE_Node
    H_KE_Node.resize(NumThreads,0);
    RecordsOnNode.resize(NumThreads,0);
    DerivRecords.resize(NumThreads,0);
    HDerivRecords.resize(NumThreads,0);
  }
  
  app_log() << "   Loading configuration from MCWalkerConfiguration::SampleStack " << endl;
  app_log() << "    number of walkers before load " << W.getActiveWalkers() << endl;
  app_log()<<"  Using Nonlocal PP in Opt: "<<includeNonlocalH<<endl;
  OhmmsInfo::Log->turnoff();
  OhmmsInfo::Warn->turnoff();
  //#pragma omp parallel for
  for (int ip=0; ip<NumThreads; ++ip)
  {
    if (H_KE_Node[ip]==0)
    {
      H_KE_Node[ip]= new QMCHamiltonian;
      H_KE_Node[ip]->addOperator(hClones[ip]->getHamiltonian("Kinetic"),"Kinetic");
      if (includeNonlocalH.c_str()!="no")
      {
        QMCHamiltonianBase* a=hClones[ip]->getHamiltonian(includeNonlocalH.c_str());
        if(a)
        {
          app_log()<<" Found non-local Hamiltonian element named "<<includeNonlocalH<<endl;
          H_KE_Node[ip]->addOperator(a,includeNonlocalH.c_str());
        }
        else
          app_log()<<" Did not find non-local Hamiltonian element named "<<includeNonlocalH<<endl;
      }
    }
  }

  app_log() << "    number of walkers before load: ";
  nVMCWalkers.resize(NumThreads);
  for (int ip=0; ip<NumThreads; ++ip)
  {
    nVMCWalkers[ip]=wClones[ip]->getActiveWalkers();
    app_log() <<  wClones[ip]->getActiveWalkers() <<  " " ;
  }
  app_log() << endl;

#pragma omp parallel for
  for (int ip=0; ip<NumThreads; ++ip)
  {
    wClones[ip]->loadEnsemble();
  }

  //load walkers from SampleStack
  OhmmsInfo::Log->reset();
  OhmmsInfo::Warn->reset();
  app_log() << "    number of walkers after load: ";
  wPerNode[0]=0;
  for (int ip=0; ip<NumThreads; ++ip)
  {
    wPerNode[ip+1]=wPerNode[ip]+wClones[ip]->getActiveWalkers();
    app_log() <<  wClones[ip]->getActiveWalkers() <<  " " ;
  }
  app_log() << endl;
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
  //JNKIM TO JEREMY
  //for(int ip=1; ip<NumThreads;++ip)
  //{
  //  opt_variables_type dummy;
  //  psiClone->checkInVariables(dummy);
  //  PsiClone[ip]->checkOutVariables(OptVariablesForPsi);
  //}
}

/** evaluate everything before optimization */
void QMCCostFunctionOMP::checkConfigurations()
{
  /* mmorales:
     Since there are cases when memory is an issue (too many dets), the use of a buffer
     is decoupled from the use of includeNonlocalH in the cost function. Without a buffer,
     everything is recalculated.
  Options:
  - "yes" or "all"  : store everything
  - "minimum"       : store orbitals and inverses only, recalculate dets
  FIX FIX FIX: right now, there is no way to include the nonlocalH in the cost function
  without using a buffer because evaluateLog needs to be called to get the inverse of psiM
  This implies that isOptimizable must be set to true, which is risky. Fix this somehow
  */
  StoreDerivInfo=false;
  DerivStorageLevel=-1;
  if(usebuffer == "yes" || usebuffer == "all")
  {
    StoreDerivInfo=true;
    if(includeNonlocalH!="no")
      DerivStorageLevel=0;
    else
      DerivStorageLevel=1;
    app_log() <<"Using buffers for temporary storage in QMCCostFunction.\n" <<endl;
  }
  else if (usebuffer == "minimum")
  {
    StoreDerivInfo=true;
    // in this case the use of nonlocalH is irrelevant, since the same inf is enough for both cases
    DerivStorageLevel=2;
    app_log() <<"Using minimum storage for determinant evaluation. \n";
  }
  else
  {
    if(includeNonlocalH!="no")
    {
      APP_ABORT("Need to enable the use of includeNonlocalH=='name' without a buffer.");
    }
  }
  int numW = 0;
  for(int i=0; i<wClones.size(); i++)
    numW += wClones[i]->getActiveWalkers();
  app_log() <<"Memory usage: " <<endl;
  app_log() <<"Linear method (approx matrix usage: 4*N^2): " <<NumParams()*NumParams()*sizeof(QMCTraits::RealType)*4.0/1.0e6  <<" MB" <<endl; // assuming 4 matrices
  app_log() <<"Deriv,HDerivRecord:      " <<numW*NumOptimizables*sizeof(QMCTraits::RealType)*3.0/1.0e6 <<" MB" <<endl;
  if(StoreDerivInfo)
  {
    MCWalkerConfiguration& dummy(*wClones[0]);
    long memorb=0,meminv=0,memdets=0,memorbs_only=0;
    Psi.memoryUsage_DataForDerivatives(dummy,memorbs_only,memorb,meminv,memdets);
    memorbs_only*=sizeof(QMCTraits::RealType);
    memorb*=sizeof(QMCTraits::RealType);
    meminv*=sizeof(QMCTraits::RealType);
    memdets*=sizeof(QMCTraits::RealType);
    app_log() <<"Buffer memory cost:     MB/walker       MB/total " <<endl;
    app_log() <<"Orbitals only:           " <<memorbs_only/1.0e6 <<"      " <<memorbs_only*numW/1.0e6 <<endl;
    app_log() <<"Orbitals + dervs:        " <<memorb/1.0e6 <<"      " <<memorb*numW/1.0e6 <<endl;
    app_log() <<"Inverse:                 " <<meminv/1.0e6 <<"      " <<meminv*numW/1.0e6 <<endl;
    app_log() <<"Determinants:            " <<memdets/1.0e6 <<"      " <<memdets*numW/1.0e6 <<endl;
  }

  app_log().flush();

  RealType et_tot=0.0;
  RealType e2_tot=0.0;
#pragma omp parallel reduction(+:et_tot,e2_tot)
  {
    int ip = omp_get_thread_num();
    MCWalkerConfiguration& wRef(*wClones[ip]);
    if (RecordsOnNode[ip] ==0)
    {
      RecordsOnNode[ip]=new Matrix<Return_t>;
      RecordsOnNode[ip]->resize(wRef.getActiveWalkers(),SUM_INDEX_SIZE);
      if (needGrads)
      {
        DerivRecords[ip]=new Matrix<Return_t>;
        DerivRecords[ip]->resize(wRef.getActiveWalkers(),NumOptimizables);
        HDerivRecords[ip]=new Matrix<Return_t>;
        HDerivRecords[ip]->resize(wRef.getActiveWalkers(),NumOptimizables);
      }
    }
    else if (RecordsOnNode[ip]->size1()!=wRef.getActiveWalkers())
    {
      RecordsOnNode[ip]->resize(wRef.getActiveWalkers(),SUM_INDEX_SIZE);
      if (needGrads)
      {
        DerivRecords[ip]->resize(wRef.getActiveWalkers(),NumOptimizables);
        HDerivRecords[ip]->resize(wRef.getActiveWalkers(),NumOptimizables);
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
    //int totalElements=W.getTotalNum()*OHMMS_DIM;
    typedef MCWalkerConfiguration::Walker_t Walker_t;
    Return_t e0=0.0;
    //       Return_t ef=0.0;
    Return_t e2=0.0;
    for (int iw=0, iwg=wPerNode[ip]; iw<wRef.getActiveWalkers(); ++iw,++iwg)
    {
      ParticleSet::Walker_t& thisWalker(*wRef[iw]);
      wRef.R=thisWalker.R;
      wRef.update();
      Return_t* restrict saved=(*RecordsOnNode[ip])[iw];
      //          Return_t logpsi(0);
      //          psiClones[ip]->evaluateDeltaLog(wRef, saved[LOGPSI_FIXED], saved[LOGPSI_FREE], *dLogPsi[iwg],*d2LogPsi[iwg]);
      // buffer for MultiSlaterDet data
      //          if((usebuffer=="yes")||(includeNonlocalH=="yes"))
      if(StoreDerivInfo)
      {
        psiClones[ip]->registerDataForDerivatives(wRef, thisWalker.DataSetForDerivatives,DerivStorageLevel);
        psiClones[ip]->evaluateDeltaLog(wRef, saved[LOGPSI_FIXED], saved[LOGPSI_FREE], *dLogPsi[iwg], *d2LogPsi[iwg], thisWalker.DataSetForDerivatives);
        //            logpsi = saved[LOGPSI_FIXED] + saved[LOGPSI_FREE];
      }
      else
      {
        psiClones[ip]->evaluateDeltaLog(wRef, saved[LOGPSI_FIXED], saved[LOGPSI_FREE], *dLogPsi[iwg], *d2LogPsi[iwg]);
        //            logpsi = psiClones[ip]->evaluateLog(wRef);
      }
      //          if(includeNonlocalH!="no") logpsi = saved[LOGPSI_FIXED] + saved[LOGPSI_FREE];
      //if (includeNonlocalH!="no")
      //  saved[ENERGY_FIXED] = hClones[ip]->getLocalPotential() - (*(hClones[ip]->getHamiltonian(includeNonlocalH.c_str()))).Value;
      //else
      //  saved[ENERGY_FIXED] = hClones[ip]->getLocalPotential();
      //           ef += saved[ENERGY_FIXED];
      saved[REWEIGHT]=thisWalker.Weight=1.0;
      //          thisWalker.resetProperty(logpsi,psiClones[ip]->getPhase(),x);
      Return_t etmp;
      if (needGrads)
      {
        //allocate vector
        vector<Return_t> Dsaved(NumOptimizables,0.0);
        vector<Return_t> HDsaved(NumOptimizables,0.0);
        psiClones[ip]->evaluateDerivatives(wRef, OptVariablesForPsi, Dsaved, HDsaved);
        etmp =hClones[ip]->evaluateValueAndDerivatives(wRef,OptVariablesForPsi,Dsaved,HDsaved,compute_nlpp);
        std::copy(Dsaved.begin(),Dsaved.end(),(*DerivRecords[ip])[iw]);
        std::copy(HDsaved.begin(),HDsaved.end(),(*HDerivRecords[ip])[iw]);
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
  OptVariablesForPsi.setComputed();
  //     app_log() << "  VMC Efavg = " << eft_tot/static_cast<Return_t>(wPerNode[NumThreads]) << endl;
  //Need to sum over the processors
  vector<Return_t> etemp(3);
  etemp[0]=et_tot;
  etemp[1]=static_cast<Return_t>(wPerNode[NumThreads]);
  etemp[2]=e2_tot;
  myComm->allreduce(etemp);
  Etarget = static_cast<Return_t>(etemp[0]/etemp[1]);
  NumSamples = static_cast<int>(etemp[1]);
  app_log() << "  VMC Eavg = " << Etarget << endl;
  app_log() << "  VMC Evar = " << etemp[2]/etemp[1]-Etarget*Etarget << endl;
  app_log() << "  Total weights = " << etemp[1] << endl;

  app_log().flush();

  setTargetEnergy(Etarget);
  ReportCounter=0;
}

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
    for (int i=0; i<psiClones.size(); ++i)
      psiClones[i]->stopOptimization();
    #pragma omp parallel
    {
      int ip = omp_get_thread_num();
      MCWalkerConfiguration& wRef(*wClones[ip]);
      MCWalkerConfiguration::iterator it(wRef.begin());
      MCWalkerConfiguration::iterator it_end(wRef.end());
      for (; it!=it_end; ++it)
        (**it).DataSetForDerivatives.clear();
    }
    // is this correct with OMP?
    //       MCWalkerConfiguration::iterator it(W.begin());
    //       MCWalkerConfiguration::iterator it_end(W.end());
    //       for (; it!=it_end; ++it)
    //         (**it).DataSetForDerivatives.clear();
  }
  //cout << "######### QMCCostFunctionOMP::resetPsi " << endl;
  //OptVariablesForPsi.print(cout);
  //cout << "-------------------------------------- " << endl;
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

  Return_t wgt_tot=0.0;
  Return_t wgt_tot2=0.0;
  Return_t NSm1 = 1.0/NumSamples;
  Return_t inv_n_samples=1.0/NumSamples;
  typedef MCWalkerConfiguration::Walker_t Walker_t;
#pragma omp parallel reduction(+:wgt_tot,wgt_tot2)
  {
    int ip = omp_get_thread_num();
    bool compute_nlpp=useNLPPDeriv && (includeNonlocalH != "no");
    MCWalkerConfiguration& wRef(*wClones[ip]);
    Return_t wgt_node=0.0, wgt_node2=0.0;
    //int totalElements=W.getTotalNum()*OHMMS_DIM;
    MCWalkerConfiguration::iterator it(wRef.begin());
    MCWalkerConfiguration::iterator it_end(wRef.end());
    int iw=0,iwg=wPerNode[ip];
    for (; it!= it_end; ++it,++iw,++iwg)
    {
      ParticleSet::Walker_t& thisWalker(**it);
      wRef.R=thisWalker.R;
      wRef.update();
      Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
      // buffer for MultiSlaterDet data
      Return_t logpsi;
      //           Return_t logpsi_old = thisWalker.getPropertyBase()[LOGPSI];
      //          if((usebuffer=="yes")||(includeNonlocalH=="yes"))
      if(StoreDerivInfo)
      {
        Walker_t::Buffer_t& tbuffer=thisWalker.DataSetForDerivatives;
        logpsi=psiClones[ip]->evaluateDeltaLog(wRef,tbuffer);
        wRef.G += *dLogPsi[iwg];
        wRef.L += *d2LogPsi[iwg];
      }
      else
      {
        logpsi=psiClones[ip]->evaluateDeltaLog(wRef);
        wRef.G += *dLogPsi[iwg];
        wRef.L += *d2LogPsi[iwg];
        //             logpsi=psiClones[ip]->evaluateLog(wRef);
      }
      //          Return_t weight = std::exp(2.0*(logpsi-saved[LOGPSI_FREE]));
      Return_t weight = saved[REWEIGHT] = vmc_or_dmc*(logpsi-saved[LOGPSI_FREE])+std::log(thisWalker.Weight);
      //          if(std::isnan(weight)||std::isinf(weight)) weight=0;
      if (needGrad)
      {
        vector<Return_t> Dsaved(NumOptimizables,0);
        vector<Return_t> HDsaved(NumOptimizables,0);
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
  //    app_log()<<"Before Purge"<<wgt_tot<<" "<<wgt_tot2<<endl;
  Return_t wgtnorm = (wgt_tot==0)?0:wgt_tot;
  wgt_tot=0.0;
  for (int ip=0; ip<NumThreads; ip++)
  {
    int nw=wClones[ip]->getActiveWalkers();
    for (int iw=0; iw<nw; iw++)
    {
      Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
      saved[REWEIGHT] = std::min(std::exp(saved[REWEIGHT]-wgtnorm), numeric_limits<Return_t>::max()*0.1 );
      wgt_tot+= inv_n_samples*saved[REWEIGHT];
    }
  }
  myComm->allreduce(wgt_tot);
  //    app_log()<<"During Purge"<<wgt_tot<<" "<<endl;
  wgtnorm =(wgt_tot==0)?1:1.0/wgt_tot;
  wgt_tot=0.0;
  for (int ip=0; ip<NumThreads; ip++)
  {
    int nw=wClones[ip]->getActiveWalkers();
    for (int iw=0; iw<nw; iw++)
    {
      Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
      saved[REWEIGHT] = std::min(saved[REWEIGHT]*wgtnorm,MaxWeight);
      wgt_tot+= inv_n_samples*saved[REWEIGHT];
    }
  }
  myComm->allreduce(wgt_tot);
  //    app_log()<<"After Purge"<<wgt_tot<<" "<<endl;
  for (int i=0; i<SumValue.size(); i++)
    SumValue[i]=0.0;
  CSWeight=wgt_tot=(wgt_tot==0)?1:1.0/wgt_tot;
  for (int ip=0; ip<NumThreads; ip++)
  {
    int nw=wClones[ip]->getActiveWalkers();
    for (int iw=0; iw<nw; iw++)
    {
      const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
      Return_t weight=saved[REWEIGHT]*wgt_tot;
      Return_t eloc_new=saved[ENERGY_NEW];
      Return_t delE=std::pow(abs(eloc_new-EtargetEff),PowerE);
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
  //     for (int i=0; i<SumValue.size(); i++) cerr<<SumValue[i]<<"  ";
  //     cerr<<endl;
  //     app_log()<<"After purge Energy Variance Weight "
  //      << SumValue[SUM_E_WGT]/SumValue[SUM_WGT] << " "
  //      << SumValue[SUM_ESQ_WGT]/SumValue[SUM_WGT] -(SumValue[SUM_E_WGT]/SumValue[SUM_WGT])*(SumValue[SUM_E_WGT]/SumValue[SUM_WGT]) << " "
  //      << SumValue[SUM_WGT]*SumValue[SUM_WGT]/SumValue[SUM_WGTSQ] << endl;
  return SumValue[SUM_WGT]*SumValue[SUM_WGT]/SumValue[SUM_WGTSQ];
}

QMCCostFunctionOMP::Return_t QMCCostFunctionOMP::fillOverlapHamiltonianMatrices(Matrix<Return_t>& H2, Matrix<Return_t>& Hamiltonian, Matrix<Return_t>& Variance, Matrix<Return_t>& Overlap)
{
  //     resetPsi();
  //     Return_t NWE = NumWalkersEff=correlatedSampling(true);
  curAvg_w = SumValue[SUM_E_WGT]/SumValue[SUM_WGT];
  Return_t curAvg2_w = SumValue[SUM_ESQ_WGT]/SumValue[SUM_WGT];
  vector<Return_t> D_avg(NumParams(),0);
  Return_t wgtinv = 1.0/SumValue[SUM_WGT];
  for (int ip=0; ip<NumThreads; ip++)
  {
    int nw=wClones[ip]->getActiveWalkers();
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
  ///zero out matrices before we start
  for (int pm=0; pm<NumParams()+1; pm++)
  {
    for (int pm2=0; pm2<NumParams()+1; pm2++)
    {
      Overlap(pm,pm2)=0;
      Hamiltonian(pm,pm2)=0;
      H2(pm,pm2)=0;
      Variance(pm,pm2)=0;
    }
  }
  for (int ip=0; ip<NumThreads; ip++)
  {
    int nw=wClones[ip]->getActiveWalkers();
    for (int iw=0; iw<nw; iw++)
    {
      const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
      Return_t weight=saved[REWEIGHT]*wgtinv;
      Return_t eloc_new=saved[ENERGY_NEW];
      const Return_t* Dsaved= (*DerivRecords[ip])[iw];
      const Return_t* HDsaved= (*HDerivRecords[ip])[iw];
      for (int pm=0; pm<NumParams(); pm++)
      {
        Return_t wfe = (HDsaved[pm] + Dsaved[pm]*(eloc_new - curAvg_w) )*weight;
        Return_t wfm = (HDsaved[pm] - 2.0*Dsaved[pm]*(eloc_new - curAvg_w) )*weight;
        Return_t wfd = (Dsaved[pm]-D_avg[pm])*weight;
        H2(0,pm+1) += wfe*(eloc_new);
        H2(pm+1,0) += wfe*(eloc_new);
        Return_t vterm = HDsaved[pm]*(eloc_new-curAvg_w)+(eloc_new*eloc_new-curAvg2_w)*Dsaved[pm]-2.0*curAvg_w*Dsaved[pm]*(eloc_new - curAvg_w);
        Variance(0,pm+1) += vterm*weight;
        Variance(pm+1,0) += vterm*weight;
        Hamiltonian(0,pm+1) += wfe;
        Hamiltonian(pm+1,0) += wfd*(eloc_new-curAvg_w);
        for (int pm2=0; pm2<NumParams(); pm2++)
        {
          H2(pm+1,pm2+1) += wfe*(HDsaved[pm2]+ Dsaved[pm2]*(eloc_new - curAvg_w));
          Hamiltonian(pm+1,pm2+1) += wfd*(HDsaved[pm2]+ Dsaved[pm2]*(eloc_new-curAvg_w));
          Variance(pm+1,pm2+1) += wfm*(HDsaved[pm2] - 2.0*Dsaved[pm2]*(eloc_new - curAvg_w));
          Overlap(pm+1,pm2+1) += wfd*(Dsaved[pm2]-D_avg[pm2]);
        }
      }
    }
  }
  myComm->allreduce(Hamiltonian);
  myComm->allreduce(Overlap);
  myComm->allreduce(Variance);
  myComm->allreduce(H2);
  Hamiltonian(0,0) = curAvg_w;
  Overlap(0,0) = 1.0;
  H2(0,0) = curAvg2_w;
  Variance(0,0) = curAvg2_w - curAvg_w*curAvg_w;
  for (int pm=1; pm<NumParams()+1; pm++)
    for (int pm2=1; pm2<NumParams()+1; pm2++)
      Variance(pm,pm2) += Variance(0,0)*Overlap(pm,pm2);
  return 1.0;
}

QMCCostFunctionOMP::Return_t
QMCCostFunctionOMP::fillOverlapHamiltonianMatrices(Matrix<Return_t>& Left, Matrix<Return_t>& Right, Matrix<Return_t>& Overlap)
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
  Overlap=0.0;

  //     resetPsi();
  //     Return_t NWE = NumWalkersEff=correlatedSampling(true);
  curAvg_w = SumValue[SUM_E_WGT]/SumValue[SUM_WGT];
  Return_t curAvg2_w = SumValue[SUM_ESQ_WGT]/SumValue[SUM_WGT];
  //    RealType H2_avg = 1.0/curAvg2_w;
  RealType H2_avg = 1.0/(curAvg_w*curAvg_w);
  //    RealType H2_avg = 1.0/std::sqrt(curAvg_w*curAvg_w*curAvg2_w);
  RealType V_avg = curAvg2_w - curAvg_w*curAvg_w;
  vector<Return_t> D_avg(NumParams(),0.0);
  Return_t wgtinv = 1.0/SumValue[SUM_WGT];
  for (int ip=0; ip<NumThreads; ip++)
  {
    int nw=wClones[ip]->getActiveWalkers();
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
    int nw=wClones[ip]->getActiveWalkers();
    for (int iw=0; iw<nw; iw++)
    {
      const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
      Return_t weight=saved[REWEIGHT]*wgtinv;
      Return_t eloc_new=saved[ENERGY_NEW];
      const Return_t* Dsaved= (*DerivRecords[ip])[iw];
      const Return_t* HDsaved= (*HDerivRecords[ip])[iw];
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
          Overlap(pm+1,pm2+1) += ovlij;
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
  myComm->allreduce(Overlap);
  Left(0,0) = (1-b2)*curAvg_w + b2*V_avg;
  Overlap(0,0) = Right(0,0) = 1.0+b1*H2_avg*V_avg;
  if (GEVType=="H2")
    return H2_avg;

  return 1.0;
}
}
/***************************************************************************
* $RCSfile$   $Author: jnkim $
* $Revision: 1898 $   $Date: 2007-04-17 10:07:34 -0500 (Tue, 17 Apr 2007) $
* $Id: QMCCostFunctionOMP.cpp 1898 2007-04-17 15:07:34Z jnkim $
***************************************************************************/
