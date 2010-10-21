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
    CSWeight=(1.0);
    app_log()<<" Using QMCCostFunctionOMP::QMCCostFunctionOMP"<<endl;
  }


  /** Clean up the vector */
  QMCCostFunctionOMP::~QMCCostFunctionOMP()
  {
    delete_iter(H_KE_Node.begin(),H_KE_Node.end());
  }


  void QMCCostFunctionOMP::GradCost(vector<Return_t>& PGradient, vector<Return_t> PM, Return_t FiniteDiff)
  {
    if (FiniteDiff != 0)
      {
        QMCTraits::RealType dh=1.0/(2.0*FiniteDiff);
        for (int i=0; i<NumOptimizables ; i++)
          {
            for (int j=0; j<NumOptimizables; j++) OptVariables[j]=PM[j];
            OptVariables[i] = PM[i]+ FiniteDiff;
            QMCTraits::RealType CostPlus = this->Cost();
            OptVariables[i] = PM[i]- FiniteDiff;
            QMCTraits::RealType CostMinus = this->Cost();
            PGradient[i]= (CostPlus-CostMinus)*dh;
          }
      }
    else
      {
        for (int j=0; j<NumOptimizables; j++) OptVariables[j]=PM[j];
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
        Return_t delE_bar;
        for (int ip=0, wn=0; ip<NumThreads; ip++)
          {
            int nw=wClones[ip]->getActiveWalkers();
            for (int iw=0; iw<nw;iw++,wn++)
              {
                const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
                Return_t weight=saved[REWEIGHT]*wgtinv;
                
                Return_t eloc_new=saved[ENERGY_NEW];
                delE_bar += weight*std::pow(abs(eloc_new-EtargetEff),PowerE);
                const Return_t* Dsaved= (*DerivRecords[ip])[iw];
                const Return_t* HDsaved= (*HDerivRecords[ip])[iw];
                for (int pm=0; pm<NumOptimizables;pm++)
                  {
                    HD_avg[pm]+= HDsaved[pm];
                  }
              }
          }
          
        myComm->allreduce(HD_avg);
        myComm->allreduce(delE_bar);
        for (int pm=0; pm<NumOptimizables;pm++)  HD_avg[pm] *= 1.0/static_cast<Return_t>(NumSamples);

        for (int ip=0, wn=0; ip<NumThreads; ip++)
          {
            int nw=wClones[ip]->getActiveWalkers();
            for (int iw=0; iw<nw;iw++,wn++)
              {
                const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
                Return_t weight=saved[REWEIGHT]*wgtinv;
                Return_t eloc_new=saved[ENERGY_NEW];
                Return_t delta_l = (eloc_new-curAvg_w);
                bool ltz(true);
                if (eloc_new-EtargetEff<0) ltz=false;
                Return_t delE=std::pow(abs(eloc_new-EtargetEff),PowerE);
                Return_t ddelE = PowerE*std::pow(abs(eloc_new-EtargetEff),PowerE-1);
                const Return_t* Dsaved= (*DerivRecords[ip])[iw];
                const Return_t* HDsaved= (*HDerivRecords[ip])[iw];
                for (int pm=0; pm<NumOptimizables;pm++)
                  {
                    EDtotals_w[pm] += weight*(HDsaved[pm] + 2.0*Dsaved[pm]*delta_l);
                    URV[pm] += 2.0*(eloc_new*HDsaved[pm] - curAvg*HD_avg[pm]);
                    if (ltz) EDtotals[pm]+= weight*(2.0*Dsaved[pm]*(delE-delE_bar) + ddelE*HDsaved[pm]);
                    else EDtotals[pm] += weight*(2.0*Dsaved[pm]*(delE-delE_bar) - ddelE*HDsaved[pm]);
                  }
              }
          }
        myComm->allreduce(EDtotals);
        myComm->allreduce(EDtotals_w);
        myComm->allreduce(URV);
        Return_t smpinv=1.0/static_cast<Return_t>(NumSamples);

        for (int ip=0, wn=0; ip<NumThreads; ip++)
          {
            int nw=wClones[ip]->getActiveWalkers();
            for (int iw=0; iw<nw;iw++,wn++)
              {
                const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
                Return_t weight=saved[REWEIGHT]*wgtinv;
                Return_t eloc_new=saved[ENERGY_NEW];
                Return_t delta_l = (eloc_new-curAvg_w);
                Return_t sigma_l = delta_l*delta_l;
                const Return_t* Dsaved= (*DerivRecords[ip])[iw];
                const Return_t* HDsaved= (*HDerivRecords[ip])[iw];
                for (int pm=0; pm<NumOptimizables;pm++)
                  {
                    E2Dtotals_w[pm] += weight*2.0*(Dsaved[pm]*(sigma_l-curVar_w) + delta_l*(HDsaved[pm]-EDtotals_w[pm]));
                  }
              }
          }
        myComm->allreduce(E2Dtotals_w);

        for (int pm=0; pm<NumOptimizables;pm++)  URV[pm] *=smpinv;
        for (int j=0; j<NumOptimizables; j++)
          {
            PGradient[j] = 0.0;
            if (std::fabs(w_var) > 1.0e-10)   PGradient[j] += w_var*E2Dtotals_w[j];
            if (std::fabs(w_en)  > 1.0e-10)   PGradient[j] += w_en*EDtotals_w[j];
            if (std::fabs(w_w)   > 1.0e-10)   PGradient[j] += w_w*URV[j];
            if (std::fabs(w_abs) > 1.0e-10)   PGradient[j] += w_abs*EDtotals[j];
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


  void
  QMCCostFunctionOMP::getConfigurations(const string& aroot)
  {

    //makeClones(W,Psi,H);
    if (H_KE_Node.empty())
      {
        app_log() << "  QMCCostFunctionOMP is created with " << NumThreads << endl;
        //make H_KE_Node
        H_KE_Node.resize(NumThreads,0);
        RecordsOnNode.resize(NumThreads,0);
        //DerivRecords.resize(NumThreads );
        //HDerivRecords.resize(NumThreads );
        DerivRecords.resize(NumThreads,0);
        HDerivRecords.resize(NumThreads,0);
      }

    app_log() << "   Loading configuration from MCWalkerConfiguration::SampleStack " << endl;
    app_log() << "    number of walkers before load " << W.getActiveWalkers() << endl;

    OhmmsInfo::Log->turnoff();
    OhmmsInfo::Warn->turnoff();
    // #pragma omp parallel for
    for (int ip=0; ip<NumThreads; ++ip)
      {
        if (H_KE_Node[ip]==0)
          {
            H_KE_Node[ip]= new QMCHamiltonian;
            H_KE_Node[ip]->addOperator(hClones[ip]->getHamiltonian("Kinetic"),"Kinetic");
            if (( includeNonlocalH=="yes")&&(hClones[ip]->getHamiltonian("NonLocalECP")))
            {
              H_KE_Node[ip]->addOperator(hClones[ip]->getHamiltonian("NonLocalECP"),"NonLocalECP");
            }
          }
        wClones[ip]->loadEnsemble();
      }
    OhmmsInfo::Log->reset();
    OhmmsInfo::Warn->reset();

    app_log() << "    number of walkers after load: ";
    for (int ip=0; ip<NumThreads; ++ip)
      app_log() <<  wClones[ip]->getActiveWalkers() <<  " " ;
    app_log() << endl;
    FairDivideLow(W.getActiveWalkers()*NumThreads,NumThreads,wPerNode);
    
    app_log()<<"  Using Nonlocal PP in Opt: "<<includeNonlocalH<<endl;

    if (dLogPsi.size() != wPerNode[NumThreads])
      {
        delete_iter(dLogPsi.begin(),dLogPsi.end());
        delete_iter(d2LogPsi.begin(),d2LogPsi.end());
        int nptcl=W.getTotalNum();
        int nwtot=wPerNode[NumThreads];
        dLogPsi.resize(nwtot,0);
        d2LogPsi.resize(nwtot,0);
        for (int i=0; i<nwtot; ++i) dLogPsi[i]=new ParticleGradient_t(nptcl);
        for (int i=0; i<nwtot; ++i) d2LogPsi[i]=new ParticleLaplacian_t(nptcl);
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
  void
  QMCCostFunctionOMP::checkConfigurations()
  {

    RealType et_tot=0.0;
    RealType eft_tot=0.0;
    RealType e2_tot=0.0;
    RealType enl_tot=0.0;
    RealType TotLogPsi=0.0;

    //#pragma omp parallel reduction(+:et_tot,e2_tot)
#pragma omp parallel
    {
      int ip = omp_get_thread_num();
      MCWalkerConfiguration& wRef(*wClones[ip]);
      if (RecordsOnNode[ip] ==0)
        {
          RecordsOnNode[ip]=new Matrix<Return_t>;
          RecordsOnNode[ip]->resize(wRef.getActiveWalkers(),SUM_INDEX_SIZE);
          DerivRecords[ip]=new Matrix<Return_t>;
          DerivRecords[ip]->resize(wRef.getActiveWalkers(),NumOptimizables);
          HDerivRecords[ip]=new Matrix<Return_t>;
          HDerivRecords[ip]->resize(wRef.getActiveWalkers(),NumOptimizables);        
        }
       else if (RecordsOnNode[ip]->size1()!=wRef.getActiveWalkers())
        {
          RecordsOnNode[ip]->resize(wRef.getActiveWalkers(),SUM_INDEX_SIZE);
          DerivRecords[ip]->resize(wRef.getActiveWalkers(),NumOptimizables);
          HDerivRecords[ip]->resize(wRef.getActiveWalkers(),NumOptimizables);   
        }
      //set the optimization mode for the trial wavefunction
      psiClones[ip]->startOptimization();

      //int nat = wRef.getTotalNum();
      //int totalElements=W.getTotalNum()*OHMMS_DIM;
      typedef MCWalkerConfiguration::Walker_t Walker_t;
      Return_t e0=0.0;
      Return_t ef=0.0;
      Return_t e2=0.0;
      MCWalkerConfiguration::iterator it(wRef.begin());
      MCWalkerConfiguration::iterator it_end(wRef.end());
      int iw=0,iwg=wPerNode[ip];
      for (; it!=it_end; ++it,++iw,++iwg)
        {
          ParticleSet::Walker_t& thisWalker(**it);
          wRef.R=thisWalker.R;
          wRef.update();
          Return_t* restrict saved=(*RecordsOnNode[ip])[iw];
//           psiClones[ip]->evaluateDeltaLog(wRef, saved[LOGPSI_FIXED], saved[LOGPSI_FREE], *dLogPsi[iwg],*d2LogPsi[iwg]);

// buffer for MultiSlaterDet data 
          psiClones[ip]->registerDataForDerivatives(wRef, thisWalker.DataSetForDerivatives);
          //thisWalker.DataSetForDerivatives=tbuffer;
          psiClones[ip]->evaluateDeltaLog(wRef, saved[LOGPSI_FIXED], saved[LOGPSI_FREE], *dLogPsi[iwg], *d2LogPsi[iwg], thisWalker.DataSetForDerivatives);

          Return_t x= hClones[ip]->evaluate(wRef);
          e0 += saved[ENERGY_TOT] = x;
          e2 += x*x;
          if (includeNonlocalH=="yes")
            saved[ENERGY_FIXED] = hClones[ip]->getLocalPotential() - (*(hClones[ip]->getHamiltonian("NonLocalECP"))).Value;
          else 
            saved[ENERGY_FIXED] = hClones[ip]->getLocalPotential();
          ef += saved[ENERGY_FIXED];
          saved[REWEIGHT]=thisWalker.Weight=1.0;


          vector<Return_t> Dsaved(NumOptimizables);
          vector<Return_t> HDsaved(NumOptimizables);
          psiClones[ip]->evaluateDerivatives(wRef, OptVariablesForPsi, Dsaved, HDsaved);
          std::copy(Dsaved.begin(),Dsaved.end(),(*DerivRecords[ip])[iw]);
          std::copy(HDsaved.begin(),HDsaved.end(),(*HDerivRecords[ip])[iw]);
        }
        
      //add them all
#pragma omp atomic
      et_tot+=e0;
#pragma omp atomic
      e2_tot+=e2;
#pragma omp atomic
      eft_tot+=ef;
    }
    OptVariablesForPsi.setComputed();
    
    app_log() << "  VMC Efavg = " << eft_tot/static_cast<Return_t>(wPerNode[NumThreads]) << endl;
    
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

    setTargetEnergy(Etarget);

    ReportCounter=0;
  }

  void QMCCostFunctionOMP::resetPsi(bool final_reset)
  {
    if (OptVariables.size() < OptVariablesForPsi.size())
      {
        for (int i=0; i<equalVarMap.size(); ++i) OptVariablesForPsi[equalVarMap[i][0]]=OptVariables[equalVarMap[i][1]];
      }
    else
      for (int i=0; i<OptVariables.size(); ++i) OptVariablesForPsi[i]=OptVariables[i];

    if (final_reset) {
      for (int i=0; i<psiClones.size(); ++i)
        psiClones[i]->stopOptimization();
      
#pragma omp parallel
    {
      int ip = omp_get_thread_num();
      MCWalkerConfiguration& wRef(*wClones[ip]);
      
      MCWalkerConfiguration::iterator it(wRef.begin());
      MCWalkerConfiguration::iterator it_end(wRef.end());
      for (; it!=it_end; ++it) (**it).DataSetForDerivatives.clear();
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

    Return_t wgt_tot=0.0;
    Return_t wgt_tot2=0.0;
    Return_t NSm1 = 1.0/NumSamples;

    //#pragma omp parallel reduction(+:wgt_tot)
    typedef MCWalkerConfiguration::Walker_t Walker_t;
#pragma omp parallel
    {
      int ip = omp_get_thread_num();
      MCWalkerConfiguration& wRef(*wClones[ip]);
      Return_t eloc_new, wgt_node=0.0, wgt_node2=0.0;
      //int totalElements=W.getTotalNum()*OHMMS_DIM;
      MCWalkerConfiguration::iterator it(wRef.begin());
      MCWalkerConfiguration::iterator it_end(wRef.end());
      int iw=0,iwg=wPerNode[ip];
      for (; it!= it_end;++it,++iw,++iwg)
        {
          ParticleSet::Walker_t& thisWalker(**it);
          wRef.R=thisWalker.R;
          wRef.update();

          Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
          Walker_t::Buffer_t& tbuffer=thisWalker.DataSetForDerivatives;

          // buffer for MultiSlaterDet data
          Return_t logpsi=psiClones[ip]->evaluateDeltaLog(wRef,tbuffer);
          Return_t weight=saved[REWEIGHT] = std::exp(2.0*(logpsi-saved[LOGPSI_FREE])) ;
          wRef.G += *dLogPsi[iwg];
          wRef.L += *d2LogPsi[iwg];
          saved[ENERGY_NEW] = H_KE_Node[ip]->evaluate(wRef) + saved[ENERGY_FIXED];

          if (needGrad)
          {
            vector<Return_t> Dsaved(NumOptimizables);
            vector<Return_t> HDsaved(NumOptimizables);
            psiClones[ip]->evaluateDerivatives(wRef, OptVariablesForPsi, Dsaved, HDsaved);
            for( int i=0;i<NumOptimizables;i++)
              if(OptVariablesForPsi.recompute(i))
              {
               (*DerivRecords[ip])(iw,i) = Dsaved[i];
               (*HDerivRecords[ip])(iw,i) = HDsaved[i];
              }
          }
          
          wgt_node+=weight;
          wgt_node2+=weight*weight;
        }

#pragma omp atomic
      wgt_tot += wgt_node;
#pragma omp atomic
      wgt_tot2 += wgt_node2;
    }

    //this is MPI barrier
    //OHMMS::Controller->barrier();
    //collect the total weight for normalization and apply maximum weight
    myComm->allreduce(wgt_tot);
    myComm->allreduce(wgt_tot2);
    //app_log()<<"Before Purge"<<wgt_tot*wgt_tot/wgt_tot2<<endl;



    Return_t wgtnorm = (1.0*NumSamples)/wgt_tot;
    wgt_tot=0.0;
    for (int ip=0; ip<NumThreads; ip++)
      {
        int nw=wClones[ip]->getActiveWalkers();
        for (int iw=0; iw<nw;iw++)
          {
            Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
            saved[REWEIGHT] = std::min(saved[REWEIGHT]*wgtnorm,MaxWeight) ;
            wgt_tot+= saved[REWEIGHT];
          }
      }
    //app_log()<<"After Purge"<<wgt_tot*wgt_tot/new_wgt_tot2<<endl;
    myComm->allreduce(wgt_tot);

    for (int i=0; i<SumValue.size(); i++) SumValue[i]=0.0;
    CSWeight=wgt_tot=1.0/wgt_tot;
    for (int ip=0; ip<NumThreads; ip++)
      {
        int nw=wClones[ip]->getActiveWalkers();
        for (int iw=0; iw<nw;iw++)
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

//     app_log()<<"Energy After Purge   "<<SumValue[SUM_E_WGT]/SumValue[SUM_WGT]<<endl;
//     app_log()<<"Variance After Purge "<<SumValue[SUM_ESQ_WGT]/SumValue[SUM_WGT]<<endl;
//     app_log()<<"Weight After Purge   "<<SumValue[SUM_WGT]*SumValue[SUM_WGT]/SumValue[SUM_WGTSQ]<<endl;
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
    for (int ip=0, wn=0; ip<NumThreads; ip++)
      {
        int nw=wClones[ip]->getActiveWalkers();
        
        
        for (int iw=0; iw<nw;iw++,wn++)
          {
            const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
            Return_t weight=saved[REWEIGHT]*wgtinv;
            const Return_t* Dsaved= (*DerivRecords[ip])[iw];
            for (int pm=0; pm<NumParams();pm++)
              {
                D_avg[pm]+= Dsaved[pm]*weight;
              } 
          }
      }
    myComm->allreduce(D_avg);
    ///zero out matrices before we start
    for (int pm=0; pm<NumParams()+1;pm++)
      {
        for (int pm2=0; pm2<NumParams()+1;pm2++)
          {
            Overlap(pm,pm2)=0;
            Hamiltonian(pm,pm2)=0;
            H2(pm,pm2)=0;
            Variance(pm,pm2)=0;
          }
      }
      
    
    for (int ip=0, wn=0; ip<NumThreads; ip++)
      {
        int nw=wClones[ip]->getActiveWalkers();
        for (int iw=0; iw<nw;iw++,wn++)
          {
            const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
            Return_t weight=saved[REWEIGHT]*wgtinv;
            Return_t eloc_new=saved[ENERGY_NEW];
            const Return_t* Dsaved= (*DerivRecords[ip])[iw];
            const Return_t* HDsaved= (*HDerivRecords[ip])[iw];

              for (int pm=0; pm<NumParams();pm++)
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
                for (int pm2=0; pm2<NumParams();pm2++)
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
    for (int pm=1; pm<NumParams()+1;pm++)
      for (int pm2=1; pm2<NumParams()+1;pm2++)
        Variance(pm,pm2) += Variance(0,0)*Overlap(pm,pm2);
          
            
    return 1.0;
  }

}
/***************************************************************************
* $RCSfile$   $Author: jnkim $
* $Revision: 1898 $   $Date: 2007-04-17 10:07:34 -0500 (Tue, 17 Apr 2007) $
* $Id: QMCCostFunctionOMP.cpp 1898 2007-04-17 15:07:34Z jnkim $
***************************************************************************/
