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
                vector<Return_t> Dsaved= (*TempDerivRecords[ip])[iw];
                vector<Return_t> HDsaved= (*TempHDerivRecords[ip])[iw];
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
                vector<Return_t> Dsaved= (*TempDerivRecords[ip])[iw];
                vector<Return_t> HDsaved= (*TempHDerivRecords[ip])[iw];
                for (int pm=0; pm<NumOptimizables;pm++)
                {
                    EDtotals_w[pm] += weight*(HDsaved[pm] + 2.0*Dsaved[pm]*delta_l ); 
                    URV[pm] += 2.0*(eloc_new*HDsaved[pm] - curAvg*HD_avg[pm]);
                    if (ltz) EDtotals[pm]+= weight*( 2.0*Dsaved[pm]*(delE-delE_bar) + ddelE*HDsaved[pm]);
                    else EDtotals[pm] += weight*( 2.0*Dsaved[pm]*(delE-delE_bar) - ddelE*HDsaved[pm]);
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
                vector<Return_t> Dsaved= (*TempDerivRecords[ip])[iw];
                vector<Return_t> HDsaved= (*TempHDerivRecords[ip])[iw];
                for (int pm=0; pm<NumOptimizables;pm++)
                { 
                    E2Dtotals_w[pm] += weight*2.0*( Dsaved[pm]*(sigma_l-curVar_w) + delta_l*(HDsaved[pm]-EDtotals_w[pm]) ); 
                }
            }
        }
        myComm->allreduce(E2Dtotals_w); 
        
        for (int pm=0; pm<NumOptimizables;pm++)  URV[pm] *=smpinv; 
        for (int j=0; j<NumOptimizables; j++) {
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
        //set the optimization mode for the trial wavefunction
        psiClones[ip]->startOptimization();

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

void QMCCostFunctionOMP::resetPsi(bool final_reset)
{
    if (OptVariables.size() < OptVariablesForPsi.size())
    {
        for (int i=0; i<equalVarMap.size(); ++i) OptVariablesForPsi[equalVarMap[i][0]]=OptVariables[equalVarMap[i][1]];
    }
    else
        for (int i=0; i<OptVariables.size(); ++i) OptVariablesForPsi[i]=OptVariables[i];

    if (final_reset)
        for (int i=0; i<psiClones.size(); ++i)
            psiClones[i]->stopOptimization();

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
    Return_t wgt_tot2=0.0;
    Return_t NSm1 = 1.0/NumSamples;

    //#pragma omp parallel reduction(+:wgt_tot)
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
            Return_t logpsi=psiClones[ip]->evaluateDeltaLog(wRef);
            wRef.G += *dLogPsi[iwg];
            wRef.L += *d2LogPsi[iwg];

            Return_t* restrict saved = (*RecordsOnNode[ip])[iw];

            RealType KEtemp = H_KE_Node[ip]->evaluate(wRef);
            eloc_new = KEtemp + saved[ENERGY_FIXED];
            Return_t weight;
//           if (samplePsi2) weight = std::min( std::exp(2.0*(logpsi-saved[LOGPSI_FREE])),MaxWeight) ;
//           else weight = std::min( std::exp( logpsi-saved[LOGPSI_FREE] ),MaxWeight) ;
            if (samplePsi2) weight = std::exp(2.0*(logpsi-saved[LOGPSI_FREE])) ;
            else weight = std::exp( logpsi-saved[LOGPSI_FREE] ) ;
            if (KEtemp<MinKE) weight=0.0;
            //  Return_t weight = usingWeight?std::exp(2.0*(logpsi-saved[LOGPSI_FREE])):1.0;
            saved[ENERGY_NEW]=eloc_new;
            saved[REWEIGHT]=weight;

            vector<Return_t>* Dsaved= &((*TempDerivRecords[ip])[iw]);
            vector<Return_t>* HDsaved= &((*TempHDerivRecords[ip])[iw]);
            psiClones[ip]->evaluateDerivatives(wRef,KEtemp,OptVariablesForPsi,*Dsaved,*HDsaved);

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

    //app_log()<<"After After Purge"<<SumValue[SUM_WGT]*SumValue[SUM_WGT]/SumValue[SUM_WGTSQ]<<endl;

    return SumValue[SUM_WGT]*SumValue[SUM_WGT]/SumValue[SUM_WGTSQ];
}

QMCCostFunctionOMP::Return_t QMCCostFunctionOMP::fillOverlapHamiltonianMatrix(Matrix<Return_t>& Hamiltonian, Matrix<Return_t>& Overlap )
{

    resetPsi();
    Return_t NWE = NumWalkersEff=correlatedSampling();
    curAvg_w = SumValue[SUM_E_WGT]/SumValue[SUM_WGT];
    vector<Return_t> D_avg(NumParams(),0);
    Return_t wgtinv = 1.0/SumValue[SUM_WGT];
    for (int ip=0, wn=0; ip<NumThreads; ip++)
    {
        int nw=wClones[ip]->getActiveWalkers();
        for (int iw=0; iw<nw;iw++,wn++ )
        {
            const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
            Return_t weight=saved[REWEIGHT]*wgtinv;
            vector<Return_t> Dsaved= (*TempDerivRecords[ip])[iw];
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
        }
    }

    for (int ip=0, wn=0; ip<NumThreads; ip++)
    {
        int nw=wClones[ip]->getActiveWalkers();
        for (int iw=0; iw<nw;iw++,wn++ )
        {
            const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
            Return_t weight=saved[REWEIGHT]*wgtinv;
            Return_t eloc_new=saved[ENERGY_NEW];
//             Return_t ke_new=saved[ENERGY_NEW] - saved[ENERGY_FIXED];
            vector<Return_t> Dsaved= (*TempDerivRecords[ip])[iw];
            vector<Return_t> HDsaved= (*TempHDerivRecords[ip])[iw];

            for (int pm=0; pm<NumParams();pm++)
            {
                Return_t wfd = (Dsaved[pm]-D_avg[pm])*weight;
                Hamiltonian(0,pm+1) += weight*(HDsaved[pm] + Dsaved[pm]*(eloc_new-curAvg_w));
                Hamiltonian(pm+1,0) += wfd*(eloc_new-curAvg_w);
                for (int pm2=0; pm2<NumParams();pm2++)
                {
                    Hamiltonian(pm+1,pm2+1) += wfd*(HDsaved[pm2]+ Dsaved[pm2]*(eloc_new-curAvg_w) );
                    Overlap(pm+1,pm2+1) += wfd*(Dsaved[pm2]-D_avg[pm2]);
                }
            }
        }
    }
    myComm->allreduce(Hamiltonian);
    myComm->allreduce(Overlap);


    Overlap(0,0) = 1;
    Hamiltonian(0,0) = curAvg_w ;
    return NWE;
}



}
/***************************************************************************
* $RCSfile$   $Author: jnkim $
* $Revision: 1898 $   $Date: 2007-04-17 10:07:34 -0500 (Tue, 17 Apr 2007) $
* $Id: QMCCostFunctionOMP.cpp 1898 2007-04-17 15:07:34Z jnkim $
***************************************************************************/
