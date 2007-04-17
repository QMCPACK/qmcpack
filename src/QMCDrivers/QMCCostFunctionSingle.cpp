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
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
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
    MCWalkerConfiguration::iterator it(W.begin()); 
    MCWalkerConfiguration::iterator it_end(W.end()); 
    Return_t eloc_new;
    int iw=0;
    Return_t wgt_tot=0.0;

    while(it != it_end) {
      Walker_t& thisWalker(**it);
      Return_t* restrict saved = Records[iw];

      //rewind the buffer to get the data from buffer
      thisWalker.DataSet.rewind();
      //W is updated by thisWalker
      W.copyFromBuffer(thisWalker.DataSet);

      Return_t logpsi=0.0;
      //copy dL from Buffer
      thisWalker.DataSet.get(dL.begin(),dL.end());
      logpsi=Psi.evaluateDeltaLog(W);
      W.G += thisWalker.Drift;
      W.L += dL;

      eloc_new=H_KE.evaluate(W)+saved[ENERGY_FIXED];
      Return_t weight = UseWeight?std::exp(2.0*(logpsi-saved[LOGPSI_FREE])):1.0;

      saved[ENERGY_NEW]=eloc_new;
      saved[REWEIGHT]=weight;
      wgt_tot+=weight;
      ++it;
      ++iw;
    }

    OHMMS::Controller->barrier();
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
    if(aroot.size() && aroot != "invalid") {
      app_log() << "  Reading configurations from the previous qmc block" << endl;
      HDFWalkerInputCollect wReader(aroot);
      wReader.putSingle(W);
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

    dL.resize(W.getTotalNum());
    int numLocWalkers=W.getActiveWalkers();
    Records.resize(numLocWalkers,6);

    typedef MCWalkerConfiguration::Walker_t Walker_t;
    MCWalkerConfiguration::iterator it(W.begin()); 
    MCWalkerConfiguration::iterator it_end(W.end()); 
    int nat = W.getTotalNum();
    int iw=0;
    Etarget=0.0;
    while(it != it_end) {

      Walker_t& thisWalker(**it);

      //clean-up DataSet to save re-used values
      thisWalker.DataSet.clear();
      //rewind the counter
      thisWalker.DataSet.rewind();
      //MCWalkerConfiguraton::registerData add distance-table data
      W.registerData(thisWalker,thisWalker.DataSet);

      Return_t*  saved=Records[iw];
#if defined(QMC_COMPLEX)
      app_error() << " Optimization is not working with complex wavefunctions yet" << endl;
      app_error() << "  Needs to fix TrialWaveFunction::evaluateDeltaLog " << endl;
#else
      Psi.evaluateDeltaLog(W, saved[LOGPSI_FIXED], saved[LOGPSI_FREE], thisWalker.Drift, dL);
#endif
      thisWalker.DataSet.add(dL.first_address(),dL.last_address());
      Etarget += saved[ENERGY_TOT] = H.evaluate(W);
      saved[ENERGY_FIXED] = H.getInvariantEnergy();

      ++it;
      ++iw;
    }

    //Need to sum over the processors
    vector<Return_t> etemp(2);
    etemp[0]=Etarget;
    etemp[1]=static_cast<Return_t>(numLocWalkers);
    myComm->allreduce(etemp);
    Etarget = static_cast<Return_t>(etemp[0]/etemp[1]);
    NumSamples = static_cast<int>(etemp[1]);

    setTargetEnergy(Etarget);

    ReportCounter=0;
  }

  void QMCCostFunctionSingle::resetPsi()
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

    //cout << "QMCCostFunctionSingle::resetPsi " <<endl;
    //oit=OptVariablesForPsi.begin();
    //oit_end=OptVariablesForPsi.end();
    //while(oit != oit_end)
    //{
    //  cout << (*oit).first << "=" << (*oit).second << " ";
    //  ++oit;
    //}
    //cout << endl;
    Psi.resetParameters(OptVariablesForPsi);
  }

  /** Reset the Wavefunction \f$ \Psi({\bf R},{{\bf \alpha_i}}) \f$
   *@return true always
   *
   * Reset from the old parameter set \f$ {{\bf \alpha_i}} \f$ to the new parameter
   * set \f$ {{\bf \alpha_{i+1}}}\f$  
   */
  bool
  QMCCostFunctionSingle::resetWaveFunctions() {

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
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
