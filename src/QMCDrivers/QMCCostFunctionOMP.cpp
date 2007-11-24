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

#pragma omp parallel reduction(+:wgt_tot)
    {
      int ip = omp_get_thread_num();
      bool usingWeight=UseWeight;
      MCWalkerConfiguration& wRef(*wClones[ip]);
      ParticleSet::ParticleLaplacian_t dl(wRef.getTotalNum());
      typedef MCWalkerConfiguration::Walker_t Walker_t;
      MCWalkerConfiguration::iterator it(wRef.begin()); 
      MCWalkerConfiguration::iterator it_end(wRef.end()); 
      Return_t eloc_new, wgt_node=0.0;
      int iw=0;
      while(it != it_end) {
        Walker_t& thisWalker(**it);

        //rewind the buffer to get the data from buffer
        thisWalker.DataSet.rewind();
        //W is updated by thisWalker
        wRef.copyFromBuffer(thisWalker.DataSet);

        Return_t logpsi=0.0;
        //copy dL from Buffer
        thisWalker.DataSet.get(dl.begin(),dl.end());
        logpsi=psiClones[ip]->evaluateDeltaLog(wRef);
        wRef.G += thisWalker.Drift;
        wRef.L += dl;

        Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
        eloc_new=H_KE_Node[ip]->evaluate(wRef)+saved[ENERGY_FIXED];
        Return_t weight = usingWeight?std::exp(2.0*(logpsi-saved[LOGPSI_FREE])):1.0;
        saved[ENERGY_NEW]=eloc_new;
        saved[REWEIGHT]=weight;
        wgt_node+=weight;

        ++it;
        ++iw;
      }
      //reduction
      wgt_tot += wgt_node;
    }

    //this is MPI barrier
    OHMMS::Controller->barrier();
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

// this is rather trivial operations
//#pragma omp parallel
//    {
//      int ip = omp_get_thread_num();
//      int pe=PowerE;
//      RealType e0=EtargetEff;
//      Return_t wgtnorm=1.0/wgt_tot;
//
//      vector<Return_t> sumLoc(SUM_INDEX_SIZE,0.0);
//      Return_t wgt_max=MaxWeight*wgtnorm;
//      int nw=wClones[ip]->getActiveWalkers();
//      for(int iw=0; iw<nw;iw++) {
//        const Return_t* restrict saved = (*RecordsOnNode[ip])[iw];
//        Return_t weight=saved[REWEIGHT]*wgtnorm;
//        Return_t eloc_new=saved[ENERGY_NEW];
//
//        weight = std::min(wgt_max,weight);
//        Return_t delE=std::pow(abs(eloc_new-e0),pe);
//        sumLoc[SUM_E_BARE] += eloc_new;
//        sumLoc[SUM_ESQ_BARE] += eloc_new*eloc_new;
//        sumLoc[SUM_ABSE_BARE] += delE;
//        sumLoc[SUM_E_WGT] += eloc_new*weight;
//        sumLoc[SUM_ESQ_WGT] += eloc_new*eloc_new*weight;
//        sumLoc[SUM_ABSE_WGT] += delE*weight;
//        sumLoc[SUM_WGT] += weight;
//        sumLoc[SUM_WGTSQ] += weight*weight;
//      }
//#pragma omp critical
//      {
//       for(int i=0; i<SUM_INDEX_SIZE; i++) SumValue[i]=sumLoc[i];
//      }
//    }

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
    W.loadEnsemble();
    app_log() << "    number of walkers after load " << W.getActiveWalkers() << endl;
    //if(aroot.size() && aroot != "invalid") {
    //  app_log() << "  Reading configurations from the previous qmc block" << endl;
    //  HDFWalkerInputCollect wReader(aroot);
    //  wReader.putSingle(W);
    //}

    //divide the walkers
    FairDivideLow(W.getActiveWalkers(),NumThreads,wPerNode);

#pragma omp parallel
    {
      int ip = omp_get_thread_num();
      if(H_KE_Node[ip]==0)
      {
        H_KE_Node[ip]= new QMCHamiltonian;
        H_KE_Node[ip]->addOperator(hClones[ip]->getHamiltonian("Kinetic"),"Kinetic");
      }
      if(ip)
      {
        wClones[ip]->createWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
      }
    }

    //remove walkers on other threads
    W.destroyWalkers(W.begin()+wPerNode[1],W.end());
  }

  /** evaluate everything before optimization */
  void 
  QMCCostFunctionOMP::checkConfigurations() {

    RealType et_tot=0.0;
#pragma omp parallel reduction(+:et_tot)
    {
      int ip = omp_get_thread_num();
      MCWalkerConfiguration& wRef(*wClones[ip]);
      ParticleSet::ParticleLaplacian_t dl(wRef.getTotalNum());
      int numLocWalkers=wRef.getActiveWalkers();

      if(RecordsOnNode[ip] ==0) RecordsOnNode[ip]=new Matrix<Return_t>;
      RecordsOnNode[ip]->resize(numLocWalkers,6);

      typedef MCWalkerConfiguration::Walker_t Walker_t;
      MCWalkerConfiguration::iterator it(wRef.begin()); 
      MCWalkerConfiguration::iterator it_end(wRef.end()); 
      int nat = wRef.getTotalNum();
      int iw=0;
      Return_t e0=0.0;
      while(it != it_end) {

        Walker_t& thisWalker(**it);

        //clean-up DataSet to save re-used values
        thisWalker.DataSet.clear();
        //rewind the counter
        thisWalker.DataSet.rewind();
        //MCWalkerConfiguraton::registerData add distance-table data
        wRef.registerData(thisWalker,thisWalker.DataSet);

        Return_t* restrict saved=(*RecordsOnNode[ip])[iw];
#if defined(QMC_COMPLEX)
        app_error() << " Optimization is not working with complex wavefunctions yet" << endl;
        app_error() << "  Needs to fix TrialWaveFunction::evaluateDeltaLog " << endl;
#else
        psiClones[ip]->evaluateDeltaLog(wRef, 
            saved[LOGPSI_FIXED], saved[LOGPSI_FREE], thisWalker.Drift, dl);
#endif
        thisWalker.DataSet.add(dl.first_address(),dl.last_address());
        e0 += saved[ENERGY_TOT] = hClones[ip]->evaluate(wRef);
        saved[ENERGY_FIXED] = hClones[ip]->getInvariantEnergy();

        ++it;
        ++iw;
      }
      //add them all
      et_tot+=e0;
    }

    //Need to sum over the processors
    vector<Return_t> etemp(2);
    etemp[0]=et_tot;
    etemp[1]=static_cast<Return_t>(wPerNode[NumThreads]);

    cout << "energy " << etemp[0] << " " << etemp[1] << endl;

    myComm->allreduce(etemp);
    Etarget = static_cast<Return_t>(etemp[0]/etemp[1]);
    NumSamples = static_cast<int>(etemp[1]);

    setTargetEnergy(Etarget);

    ReportCounter=0;
  }

  void QMCCostFunctionOMP::resetPsi()
  {

    OptimizableSetType::iterator oit(OptVariables.begin()), oit_end(OptVariables.end());
    while(oit != oit_end)
    {
      Return_t v=(*oit).second;
      OptVariablesForPsi[(*oit).first]=v;
      map<string,set<string>*>::iterator eit(equalConstraints.find((*oit).first));
      if(eit != equalConstraints.end())
      {
        set<string>::iterator f((*eit).second->begin()),l((*eit).second->end());
        while(f != l)
        {
          OptVariablesForPsi[(*f)]=v;
          ++f;
        }
      }
      ++oit;
    }

    //cout << "QMCCostFunctionOMP::resetPsi " <<endl;
    //oit=OptVariablesForPsi.begin();
    //oit_end=OptVariablesForPsi.end();
    //while(oit != oit_end)
    //{
    //  cout << (*oit).first << "=" << (*oit).second << " ";
    //  ++oit;
    //}
    //cout << endl;
    if(psiClones.empty())
      Psi.resetParameters(OptVariablesForPsi);
    else
      for(int i=0; i<NumThreads; i++) psiClones[i]->resetParameters(OptVariablesForPsi);
  }

  /** Reset the Wavefunction \f$ \Psi({\bf R},{{\bf \alpha_i}}) \f$
   *@return true always
   *
   * Reset from the old parameter set \f$ {{\bf \alpha_i}} \f$ to the new parameter
   * set \f$ {{\bf \alpha_{i+1}}}\f$  
   */
  bool
  QMCCostFunctionOMP::resetWaveFunctions() {

    //loop over all the unique id's
    for(int i=0; i<IDtag.size(); i++)
    {
      OptVariables[IDtag[i]]=OptParams[i];
    }
    resetPsi();
    return true;
  }
}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1898 $   $Date: 2007-04-17 10:07:34 -0500 (Tue, 17 Apr 2007) $
 * $Id: QMCCostFunctionOMP.cpp 1898 2007-04-17 15:07:34Z jnkim $ 
 ***************************************************************************/
