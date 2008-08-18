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


namespace qmcplusplus {

  QMCCostFunctionOMP::QMCCostFunctionOMP( MCWalkerConfiguration& w, 
      TrialWaveFunction& psi, QMCHamiltonian& h, HamiltonianPool& hpool):
    QMCCostFunctionBase(w,psi,h), CloneManager(hpool)
  { 
  }


  /** Clean up the vector */
  QMCCostFunctionOMP::~QMCCostFunctionOMP() {
  }

  /**  Perform the correlated sampling algorthim.
   */
  QMCCostFunctionOMP::Return_t QMCCostFunctionOMP::correlatedSampling() {

    Return_t wgt_tot=0.0;

//#pragma omp parallel reduction(+:wgt_tot)
#pragma omp parallel 
    {
      int ip = omp_get_thread_num();
      bool usingWeight=UseWeight;
      MCWalkerConfiguration& wRef(*wClones[ip]);
      ParticleSet::ParticleLaplacian_t pdL(wRef.getTotalNum());
      ParticleSet::ParticleGradient_t pdG(wRef.getTotalNum());
      typedef MCWalkerConfiguration::Walker_t Walker_t;
      MCWalkerConfiguration::iterator it(wRef.begin()); 
      MCWalkerConfiguration::iterator it_end(wRef.end()); 
      Return_t eloc_new, wgt_node=0.0;
      int totalElements=W.getTotalNum()*OHMMS_DIM;
      int iw=0;
      while(it != it_end) {
        Walker_t& thisWalker(**it);
        //rewind the buffer to get the data from buffer
        thisWalker.DataSet.rewind();
        //W is updated by thisWalker
        //BEGIN:MEMORY FIX
        wRef.R=thisWalker.R;
        wRef.update();
        //wRef.copyFromBuffer(thisWalker.DataSet);
        //END:MEMORY FIX

        thisWalker.DataSet.get(&(pdG[0][0]),&(pdG[0][0])+totalElements);
        thisWalker.DataSet.get(pdL.begin(),pdL.end());

        Return_t logpsi=psiClones[ip]->evaluateDeltaLog(wRef);
        wRef.G += pdG;
        wRef.L += pdL;

        Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
        eloc_new=H_KE_Node[ip]->evaluate(wRef)+saved[ENERGY_FIXED];
        Return_t weight = usingWeight?std::exp(2.0*(logpsi-saved[LOGPSI_FREE])):1.0;
        saved[ENERGY_NEW]=eloc_new;
        saved[REWEIGHT]=weight;
        wgt_node+=weight;

        ++it;
        ++iw;
      }
#pragma omp atomic
      wgt_tot += wgt_node;
    }

    //this is MPI barrier
    //OHMMS::Controller->barrier();
    //collect the total weight for normalization and apply maximum weight
    myComm->allreduce(wgt_tot);

    for(int i=0; i<SumValue.size(); i++) SumValue[i]=0.0;

    wgt_tot=1.0/wgt_tot;
    Return_t wgt_max=MaxWeight*wgt_tot;
    for(int ip=0; ip<NumThreads; ip++)
    {
      int nw=wClones[ip]->getActiveWalkers();
      for(int iw=0; iw<nw;iw++) {
        const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
        Return_t weight=saved[REWEIGHT]*wgt_tot;
        Return_t eloc_new=saved[ENERGY_NEW];

        weight = (weight>wgt_max)? wgt_max:weight;
        Return_t delE=std::pow(abs(eloc_new-EtargetEff),PowerE);
        SumValue[SUM_E_BARE] += eloc_new;
        SumValue[SUM_ESQ_BARE] += eloc_new*eloc_new;
        SumValue[SUM_ABSE_BARE] += delE;
        SumValue[SUM_E_WGT] += eloc_new*weight;
        SumValue[SUM_ESQ_WGT] += eloc_new*eloc_new*weight;
        SumValue[SUM_ABSE_WGT] += delE*weight;
        SumValue[SUM_WGT] += weight;
        SumValue[SUM_WGTSQ] += weight*weight;
      }
    }
    //collect everything
    myComm->allreduce(SumValue);

    return SumValue[SUM_WGT]*SumValue[SUM_WGT]/SumValue[SUM_WGTSQ];
  }

  void 
  QMCCostFunctionOMP::getConfigurations(const string& aroot) {

    //makeClones(W,Psi,H);
    if(H_KE_Node.empty())
    {
      app_log() << "  QMCCostFunctionOMP is created with " << NumThreads << endl;
      //make H_KE_Node
      H_KE_Node.resize(NumThreads,0);
      RecordsOnNode.resize(NumThreads,0);
    }

    app_log() << "   Loading configuration from MCWalkerConfiguration::SampleStack " << endl;
    app_log() << "    number of walkers before load " << W.getActiveWalkers() << endl;

#pragma omp parallel
    {
      int ip = omp_get_thread_num();
      if(H_KE_Node[ip]==0)
      {
        H_KE_Node[ip]= new QMCHamiltonian;
        H_KE_Node[ip]->addOperator(hClones[ip]->getHamiltonian("Kinetic"),"Kinetic");
      }
      wClones[ip]->destroyWalkers(wClones[ip]->begin(),wClones[ip]->end());
      wClones[ip]->loadEnsemble();
      wClones[ip]->clearEnsemble();
    }

    app_log() << "    number of walkers after load: ";
    for(int ip=0; ip<NumThreads; ++ip)
      app_log() <<  wClones[ip]->getActiveWalkers() <<  " " ;
    app_log() << endl;
    FairDivideLow(W.getActiveWalkers()*NumThreads,NumThreads,wPerNode);
  }

  /** evaluate everything before optimization */
  void 
  QMCCostFunctionOMP::checkConfigurations() {

    RealType et_tot=0.0;
    RealType e2_tot=0.0;

//#pragma omp parallel reduction(+:et_tot,e2_tot)
#pragma omp parallel
    {
      int ip = omp_get_thread_num();
      MCWalkerConfiguration& wRef(*wClones[ip]);
      ParticleSet::ParticleLaplacian_t pdL(wRef.getTotalNum());
      ParticleSet::ParticleGradient_t pdG(wRef.getTotalNum());
      int numLocWalkers=wRef.getActiveWalkers();

      if(RecordsOnNode[ip] ==0) RecordsOnNode[ip]=new Matrix<Return_t>;
      RecordsOnNode[ip]->resize(numLocWalkers,6);

      typedef MCWalkerConfiguration::Walker_t Walker_t;
      MCWalkerConfiguration::iterator it(wRef.begin()); 
      MCWalkerConfiguration::iterator it_end(wRef.end()); 
      int nat = wRef.getTotalNum();
      int totalElements=W.getTotalNum()*OHMMS_DIM;
      int iw=0;
      Return_t e0=0.0;
      Return_t e2=0.0;
      while(it != it_end) {

        Walker_t& thisWalker(**it);

        //clean-up DataSet to save re-used values
        thisWalker.DataSet.clear();
        //rewind the counter
        thisWalker.DataSet.rewind();
        //MCWalkerConfiguraton::registerData add distance-table data
        //BEGIN:MEMORY FIX
        wRef.R=thisWalker.R;
        wRef.update();
        //wRef.registerData(thisWalker.DataSet);
        //END:MEMORY FIX
        Return_t* restrict saved=(*RecordsOnNode[ip])[iw];
        psiClones[ip]->evaluateDeltaLog(wRef, saved[LOGPSI_FIXED], saved[LOGPSI_FREE], pdG, pdL);

        thisWalker.DataSet.add(&(pdG[0][0]),&(pdG[0][0])+totalElements);
        thisWalker.DataSet.add(pdL.first_address(),pdL.last_address());

        //e0 += saved[ENERGY_TOT] = hClones[ip]->evaluate(wRef);
        Return_t x= hClones[ip]->evaluate(wRef);
        e0 += saved[ENERGY_TOT] = x;
        e2 += x*x;
        saved[ENERGY_FIXED] = hClones[ip]->getLocalPotential();

        ++it;
        ++iw;
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
    if(OptVariables.size() < OptVariablesForPsi.size())
    {
      for(int i=0; i<equalVarMap.size(); ++i)
        OptVariablesForPsi[equalVarMap[i][0]]=OptVariables[equalVarMap[i][1]];
    }
    else
      for(int i=0; i<OptVariables.size(); ++i) OptVariablesForPsi[i]=OptVariables[i];

    //cout << "######### QMCCostFunctionOMP::resetPsi " << endl;
    //OptVariablesForPsi.print(cout);
    //cout << "-------------------------------------- " << endl;
    Psi.resetParameters(OptVariablesForPsi);

    for(int i=0; i<psiClones.size(); ++i)
      psiClones[i]->resetParameters(OptVariablesForPsi);

  }

}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1898 $   $Date: 2007-04-17 10:07:34 -0500 (Tue, 17 Apr 2007) $
 * $Id: QMCCostFunctionOMP.cpp 1898 2007-04-17 15:07:34Z jnkim $ 
 ***************************************************************************/
