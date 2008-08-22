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
#include "QMCDrivers/QMCCostFunctionSingle.h"
#include "Particle/MCWalkerConfiguration.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Particle/HDFWalkerInputCollect.h"
#include "Message/CommOperators.h"
#include "Utilities/Timer.h"
//#define QMCCOSTFUNCTION_DEBUG


namespace qmcplusplus {

  QMCCostFunctionSingle::QMCCostFunctionSingle(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h):
    QMCCostFunctionBase(w,psi,h)
  { 
  }

  /** Clean up the vector */
  QMCCostFunctionSingle::~QMCCostFunctionSingle() {
  }

  /**  Perform the correlated sampling algorthim.
   */
  QMCCostFunctionSingle::Return_t QMCCostFunctionSingle::correlatedSampling() {

    typedef MCWalkerConfiguration::Walker_t Walker_t;
    //Return_t eloc_new;
    //int totalElements=W.getTotalNum()*OHMMS_DIM;

    Return_t wgt_tot=0.0;
    MCWalkerConfiguration::iterator it(W.begin()); 
    MCWalkerConfiguration::iterator it_end(W.end()); 
    int iw=0;
    for(; it!=it_end;++it,++iw)
    {
      Walker_t& thisWalker(**it);
      Return_t* restrict saved = Records[iw];
      W.R=(*it)->R;
      W.update();
      Return_t logpsi=Psi.evaluateDeltaLog(W);
      W.G += *dLogPsi[iw];
      W.L += *d2LogPsi[iw];

      Return_t eloc_new=H_KE.evaluate(W)+saved[ENERGY_FIXED];
      Return_t weight = UseWeight?std::exp(2.0*(logpsi-saved[LOGPSI_FREE])):1.0;

      saved[ENERGY_NEW]=eloc_new;
      saved[REWEIGHT]=weight;
      wgt_tot+=weight;
    }

    //collect the total weight for normalization and apply maximum weight
    myComm->allreduce(wgt_tot);

    for(int i=0; i<SumValue.size(); i++) SumValue[i]=0.0;
    wgt_tot=1.0/wgt_tot;

    Return_t wgt_max=MaxWeight*wgt_tot;
    int nw=W.getActiveWalkers();
    for(iw=0; iw<nw;iw++) {
      Return_t* restrict saved = Records[iw];
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

    //collect everything
    myComm->allreduce(SumValue);

    return SumValue[SUM_WGT]*SumValue[SUM_WGT]/SumValue[SUM_WGTSQ];
  }

  void 
  QMCCostFunctionSingle::getConfigurations(const string& aroot) {
    app_log() << "   Loading configuration from MCWalkerConfiguration::SampleStack " << endl;
    app_log() << "    number of walkers before load " << W.getActiveWalkers() << endl;
    Timer t1;
    W.destroyWalkers(W.begin(),W.end());
    W.loadEnsemble();
    W.clearEnsemble();
    app_log() << "    Loading time = " << t1.elapsed() << endl;
    app_log() << "    number of walkers after load " << W.getActiveWalkers() << endl;
    //if(aroot.size() && aroot != "invalid") {
    //  app_log() << "  Reading configurations from the previous qmc block" << endl;
    //  HDFWalkerInputCollect wReader(aroot);
    //  wReader.putSingle(W);
    //}
    //
    if(dLogPsi.size() != W.getActiveWalkers())
    {
      delete_iter(dLogPsi.begin(),dLogPsi.end());
      delete_iter(d2LogPsi.begin(),d2LogPsi.end());
      int nptcl=W.getTotalNum();
      int nwtot=W.getActiveWalkers();
      dLogPsi.resize(nwtot,0);
      d2LogPsi.resize(nwtot,0);
      for(int i=0; i<nwtot; ++i) dLogPsi[i]=new ParticleGradient_t(nptcl);
      for(int i=0; i<nwtot; ++i) d2LogPsi[i]=new ParticleLaplacian_t(nptcl);
    }

#if defined(QMCCOSTFUNCTION_DEBUG)
    if(debug_stream) delete debug_stream;
    char fname[16];
    sprintf(fname,"optdebug.p%d",OHMMS::Controller->mycontext());
    debug_stream = new ofstream(fname);
    debug_stream->setf(ios::scientific, ios::floatfield);
    debug_stream->precision(8);

    *debug_stream << "Initial : " << endl;
    for(int i=0; i<OptParams.size(); i++) 
      *debug_stream << " " << IDtag[i] << "=" << OptParams[i];
    *debug_stream << endl;
#endif
  }

 // /** evaluate everything before optimization */
 // void 
 // QMCCostFunctionSingle::getConfigurations(vector<string>& ConfigFile, 
 //     int partid, int nparts) {
 //   if(ConfigFile.size()) {

 //     app_log() << "  Reading configurations from mcwalkerset " << endl;

 //     W.destroyWalkers(W.begin(),W.end());
 //     for(int i=0; i<ConfigFile.size(); i++) {
 //       //JNKIM: needs to change to HDFWalkerInputCollect
 //       //HDFWalkerInput0 wReader(ConfigFile[i],partid,nparts);
 //       HDFWalkerInputCollect wReader(ConfigFile[i]);
 //       wReader.putSingle(W);
 //       //wReader.put(W,-1);
 //     }

 //     //remove input files
 //     ConfigFile.erase(ConfigFile.begin(),ConfigFile.end());
 //   }
 // }

  /** evaluate everything before optimization */
  void 
  QMCCostFunctionSingle::checkConfigurations() {

    //dG.resize(W.getTotalNum());
    //dL.resize(W.getTotalNum());
    int numLocWalkers=W.getActiveWalkers();
    Records.resize(numLocWalkers,6);

    typedef MCWalkerConfiguration::Walker_t Walker_t;
    MCWalkerConfiguration::iterator it(W.begin()); 
    MCWalkerConfiguration::iterator it_end(W.end()); 
    //int totalElements=W.getTotalNum()*OHMMS_DIM;
    Etarget=0.0;
    Return_t e2sum=0.0;
    for(int iw=0; it!=it_end; ++it,++iw)
    {
      Walker_t& thisWalker(**it);
      W.R=thisWalker.R;
      W.update();
      Return_t*  saved=Records[iw];
      Psi.evaluateDeltaLog(W, saved[LOGPSI_FIXED], saved[LOGPSI_FREE], *dLogPsi[iw], *d2LogPsi[iw]);
      Return_t e=H.evaluate(W);
      e2sum += e*e;
      Etarget += saved[ENERGY_TOT] = e;
      saved[ENERGY_FIXED] = H.getLocalPotential();
    }

    //Need to sum over the processors
    vector<Return_t> etemp(3);
    etemp[0]=Etarget;
    etemp[1]=static_cast<Return_t>(numLocWalkers);
    etemp[2]=e2sum;

    myComm->allreduce(etemp);
    Etarget = static_cast<Return_t>(etemp[0]/etemp[1]);
    NumSamples = static_cast<int>(etemp[1]);

    app_log() << "  VMC Eavg = " << Etarget << endl;
    app_log() << "  VMC Evar = " << etemp[2]/etemp[1]-Etarget*Etarget << endl;
    app_log() << "  Total weights = " << etemp[1] << endl;

    setTargetEnergy(Etarget);

    ReportCounter=0;
  }

  void QMCCostFunctionSingle::resetPsi()
  {
    if(OptVariables.size() < OptVariablesForPsi.size())
    {
      for(int i=0; i<equalVarMap.size(); ++i)
        OptVariablesForPsi[equalVarMap[i][0]]=OptVariables[equalVarMap[i][1]];
    }
    else
      for(int i=0; i<OptVariables.size(); ++i) OptVariablesForPsi[i]=OptVariables[i];
  }

  ///** Reset the Wavefunction \f$ \Psi({\bf R},{{\bf \alpha_i}}) \f$
  // *@return true always
  // *
  // * Reset from the old parameter set \f$ {{\bf \alpha_i}} \f$ to the new parameter
  // * set \f$ {{\bf \alpha_{i+1}}}\f$  
  // */
  //bool
  //QMCCostFunctionSingle::resetWaveFunctions() {

  //  resetPsi();
  //  return true;
  //}
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
