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
    cout<<" Using QMCCostFunctionOMP::QMCCostFunctionOMP"<<endl;
  }


  /** Clean up the vector */
  QMCCostFunctionOMP::~QMCCostFunctionOMP()
  {
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
        NumWalkersEff=correlatedSampling();
        //Estimators::accumulate has been called by correlatedSampling


        curAvg_w = SumValue[SUM_E_WGT]/SumValue[SUM_WGT];
        Return_t curAvg2_w = curAvg_w*curAvg_w;
        curVar_w = SumValue[SUM_ESQ_WGT]/SumValue[SUM_WGT]-curAvg_w*curAvg_w;
        vector<Return_t> EDtotals(NumOptimizables,0.0);
        vector<Return_t> E2Dtotals(NumOptimizables,0.0);
        vector<Return_t> EDtotals_w(NumOptimizables,0.0);
        vector<Return_t> E2Dtotals_w(NumOptimizables,0.0);
        vector<Return_t> URV(NumOptimizables,0.0);
        vector<Return_t> HD_avg(NumOptimizables,0.0);
        Return_t e2(0.0);


        for (int ip=0, wn=0; ip<NumThreads; ip++)
          {
            int nw=wClones[ip]->getActiveWalkers();
            for (int iw=0; iw<nw;iw++,wn++)
              {
                const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
                Return_t weight=saved[REWEIGHT]/SumValue[SUM_WGT];
                Return_t eloc_new=saved[ENERGY_NEW];
                e2 += eloc_new*eloc_new*weight;
                vector<Return_t> Dsaved= (*TempDerivRecords[ip])[iw];
                vector<Return_t> HDsaved= (*TempHDerivRecords[ip])[iw];
                for (int pm=0; pm<NumOptimizables;pm++)
                  {
                    HD_avg[pm]+= HDsaved[pm];
                    Return_t val;
                    if (samplePsi2) val = (HDsaved[pm] + 2.0*(eloc_new - curAvg_w)*Dsaved[pm] );
                    else val = (HDsaved[pm] + (eloc_new - curAvg_w)*Dsaved[pm]);
                    EDtotals[pm] +=  val;
                    EDtotals_w[pm] += weight*val;
                  }
              }
          }
        myComm->allreduce(EDtotals);
        myComm->allreduce(EDtotals_w);
        myComm->allreduce(HD_avg);

        Return_t wgtinv=1.0/static_cast<Return_t>(NumSamples);
        curAvg = SumValue[SUM_E_BARE]*wgtinv;
        for (int pm=0; pm<NumOptimizables;pm++)  HD_avg[pm] *=wgtinv;
        Return_t weightFrac(0.0);
        for (int ip=0, wn=0; ip<NumThreads; ip++)
          {
            int nw=wClones[ip]->getActiveWalkers();
            for (int iw=0; iw<nw;iw++,wn++)
              {
                const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
                Return_t weight=saved[REWEIGHT]*CSWeight;
                Return_t eloc_new=saved[ENERGY_NEW];
                vector<Return_t> Dsaved= (*TempDerivRecords[ip])[iw];
                vector<Return_t> HDsaved= (*TempHDerivRecords[ip])[iw];
                for (int pm=0; pm<NumOptimizables;pm++)
                  {
                    URV[pm] += 2.0*(eloc_new*HDsaved[pm] - curAvg*HD_avg[pm]);
                    Return_t val;
                    if (samplePsi2) val=  2.0*((eloc_new-curAvg_w)*(eloc_new-curAvg_w)-curVar_w)*Dsaved[pm]
                                    + 2.0*(eloc_new-curAvg_w)*(HDsaved[pm] - EDtotals_w[pm]);
//                                     - 2.0*curAvg_w*EDtotals_w[pm];
                    else val = (eloc_new*eloc_new - e2)*Dsaved[pm]
                                    + 2.0*eloc_new*(HDsaved[pm])
                                    - 2.0*curAvg_w*EDtotals_w[pm];
                    E2Dtotals[pm] += val ;
                    E2Dtotals_w[pm] += val * weight;
                  }
              }
          }
        myComm->allreduce(E2Dtotals);
        myComm->allreduce(E2Dtotals_w);
        myComm->allreduce(URV);

        for (int pm=0; pm<NumOptimizables;pm++)  URV[pm] *=wgtinv;
	// app_log() << "Gradient:\n";
        for (int j=0; j<NumOptimizables; j++) {
	  PGradient[j] = 0.0;
	  if (std::fabs(w_var) > 1.0e-10)   PGradient[j] += w_var*E2Dtotals_w[j];
	  if (std::fabs(w_en)  > 1.0e-10)   PGradient[j] += w_en*EDtotals_w[j];
	  if (std::fabs(w_w)   > 1.0e-10)   PGradient[j] += w_w*URV[j];
	  
	  //	  PGradient[j] = (w_var*E2Dtotals_w[j] + w_en*EDtotals_w[j] + w_w*URV[j]) ;
	  // char buff[200];
	  // snprintf (buff, 200, "%10.5e\n", PGradient[j]);
	  // app_log() << buff;
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
        TempDerivRecords.resize(NumThreads,0);
        TempHDerivRecords.resize(NumThreads,0);
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
    RealType e2_tot=0.0;
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
          //DerivRecords[ip]=new vector<vector<Return_t> >(wRef.getActiveWalkers() , vector<Return_t>(0.0,NumOptimizables));
          //HDerivRecords[ip]=new vector<vector<Return_t> >(wRef.getActiveWalkers() , vector<Return_t>(0.0,NumOptimizables));
          TempDerivRecords[ip]=new vector<vector<Return_t> >(wRef.getActiveWalkers() , vector<Return_t>(NumOptimizables,0.0));
          TempHDerivRecords[ip]=new vector<vector<Return_t> >(wRef.getActiveWalkers() , vector<Return_t>(NumOptimizables,0.0));
        }

      //int nat = wRef.getTotalNum();
      //int totalElements=W.getTotalNum()*OHMMS_DIM;
      Return_t e0=0.0;
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
          //vector<Return_t> Dsaved= (*DerivRecords[ip])[iw];
          //vector<Return_t> HDsaved= (*HDerivRecords[ip])[iw];
          psiClones[ip]->evaluateDeltaLog(wRef, saved[LOGPSI_FIXED], saved[LOGPSI_FREE], *dLogPsi[iwg],*d2LogPsi[iwg]);
          Return_t x= hClones[ip]->evaluate(wRef);
          e0 += saved[ENERGY_TOT] = x;
          e2 += x*x;
          saved[ENERGY_FIXED] = hClones[ip]->getLocalPotential();
//           if (samplePsi2) thisWalker.Weight=1.0;
//           else thisWalker.Weight=std::exp(saved[LOGPSI_FIXED]+saved[LOGPSI_FREE]);
          thisWalker.Weight=1.0;


          vector<Return_t>* Dsaved= &((*TempDerivRecords[ip])[iw]);
          vector<Return_t>* HDsaved= &((*TempHDerivRecords[ip])[iw]);
          psiClones[ip]->evaluateDerivatives(wRef,saved[ENERGY_TOT]-saved[ENERGY_FIXED],OptVariablesForPsi,*Dsaved,*HDsaved);
          //for(int l=0; l<NumOptimizables; l++) (*HDsaved)[l] += saved[ENERGY_FIXED]* (*Dsaved)[l];
          //for(int l=0; l<NumOptimizables; l++) cout<<HDsaved[l]<<"  "<<Dsaved[l];
          //cout<<endl;
        }
      //add them all
#pragma omp atomic
      et_tot+=e0;
#pragma omp atomic
      e2_tot+=e2;
    }

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

  void QMCCostFunctionOMP::resetPsi()
  {
    if (OptVariables.size() < OptVariablesForPsi.size())
      {
        for (int i=0; i<equalVarMap.size(); ++i) OptVariablesForPsi[equalVarMap[i][0]]=OptVariables[equalVarMap[i][1]];
      }
    else
      for (int i=0; i<OptVariables.size(); ++i) OptVariablesForPsi[i]=OptVariables[i];

    //cout << "######### QMCCostFunctionOMP::resetPsi " << endl;
    //OptVariablesForPsi.print(cout);
    //cout << "-------------------------------------- " << endl;
    Psi.resetParameters(OptVariablesForPsi);

    for (int i=0; i<psiClones.size(); ++i)
      psiClones[i]->resetParameters(OptVariablesForPsi);

  }

  QMCCostFunctionOMP::Return_t QMCCostFunctionOMP::correlatedSampling()
  {

    Return_t wgt_tot=0.0;
    Return_t wgt_bare=0.0;
    Return_t wgt_bare2=0.0;
    Return_t NSm1 = 1.0/NumSamples;

    //#pragma omp parallel reduction(+:wgt_tot)
#pragma omp parallel
    {
      int ip = omp_get_thread_num();
      //bool usingWeight=UseWeight;
      MCWalkerConfiguration& wRef(*wClones[ip]);
      Return_t eloc_new, wgt_node=0.0;
//       , wgt_node_bare=0.0, wgt_bare2_node=0.0;
      //int totalElements=W.getTotalNum()*OHMMS_DIM;
      MCWalkerConfiguration::iterator it(wRef.begin());
      MCWalkerConfiguration::iterator it_end(wRef.end());
      int iw=0,iwg=wPerNode[ip];
      for (; it!= it_end;++it,++iw,++iwg)
        {
          ParticleSet::Walker_t& thisWalker(**it);
          wRef.R=thisWalker.R;
          wRef.update();
          Return_t logpsi=psiClones[ip]->evaluateDeltaLog(wRef);
          wRef.G += *dLogPsi[iwg];
          wRef.L += *d2LogPsi[iwg];

          Return_t* restrict saved = (*RecordsOnNode[ip])[iw];

          RealType KEtemp = H_KE_Node[ip]->evaluate(wRef);
          //if (KEtemp<MinKE) KEtemp=std::fabs(KEtemp-MinKE)+MinKE;
          eloc_new = KEtemp + saved[ENERGY_FIXED];
//           Return_t this_bare_weight=(1.0);
          Return_t weight;
          if (samplePsi2) weight = std::min( std::exp(2.0*(logpsi-saved[LOGPSI_FREE])),MaxWeight) ;
          else weight = std::min( std::exp( logpsi-saved[LOGPSI_FREE] ),MaxWeight) ;
          if (KEtemp<MinKE) weight=0.0;
          //  Return_t weight = usingWeight?std::exp(2.0*(logpsi-saved[LOGPSI_FREE])):1.0;
          saved[ENERGY_NEW]=eloc_new;
          saved[REWEIGHT]=weight;

          vector<Return_t>* Dsaved= &((*TempDerivRecords[ip])[iw]);
          vector<Return_t>* HDsaved= &((*TempHDerivRecords[ip])[iw]);
          psiClones[ip]->evaluateDerivatives(wRef,KEtemp,OptVariablesForPsi,*Dsaved,*HDsaved);
          //for(int l=0; l<NumOptimizables; l++) (*HDsaved)[l] += saved[ENERGY_FIXED]* (*Dsaved)[l];

          wgt_node+=weight;
//           wgt_node_bare+=this_bare_weight;
//           wgt_bare2_node+=this_bare_weight*this_bare_weight;
        }

#pragma omp atomic
      wgt_tot += wgt_node;
// #pragma omp atomic
//       wgt_bare += wgt_node_bare;
// #pragma omp atomic
//       wgt_bare2 +=wgt_bare2_node;
    }

    //this is MPI barrier
    //OHMMS::Controller->barrier();
    //collect the total weight for normalization and apply maximum weight
    myComm->allreduce(wgt_tot);
    myComm->allreduce(wgt_bare);
    myComm->allreduce(wgt_bare2);
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
//     app_log()<<"wgt_bare:"<<wgt_bare<<" wgt_bare2:"<<wgt_bare2<<"  SumValue[SUM_WGT]:"<<SumValue[SUM_WGT]<<"  SumValue[SUM_WGTSQ]:"<<SumValue[SUM_WGTSQ]<<endl;
//     return wgt_bare*wgt_bare/wgt_bare2;
//     return SumValue[SUM_WGT]*SumValue[SUM_WGT]/wgt_bare2;
    return SumValue[SUM_WGT]*SumValue[SUM_WGT]/SumValue[SUM_WGTSQ];
  }

}
/***************************************************************************
* $RCSfile$   $Author: jnkim $
* $Revision: 1898 $   $Date: 2007-04-17 10:07:34 -0500 (Tue, 17 Apr 2007) $
* $Id: QMCCostFunctionOMP.cpp 1898 2007-04-17 15:07:34Z jnkim $
***************************************************************************/
