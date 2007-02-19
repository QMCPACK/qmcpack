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
#ifndef QMCPLUSPLUS_COSTFUNCTION_H
#define QMCPLUSPLUS_COSTFUNCTION_H

#include <deque>
#include "Configuration.h"
#include "Optimize/OptimizeBase.h"
#include "Optimize/VarList.h"
#include "QMCHamiltonians/QMCHamiltonian.h"

namespace qmcplusplus {

  class MCWalkerConfiguration;
  class TrialWaveFunction;

  /** @ingroup QMCDrivers
   * @brief Implements wave-function optimization
   *
   * Optimization by correlated sampling method with configurations 
   * generated from VMC.
   */

  class QMCCostFunction: public CostFunctionBase<OHMMS_PRECISION>
  {
  public:

    typedef VarRegistry<Return_t> OptimizableSetType;

    enum FieldIndex_OPT {LOGPSI_FIXED=0, LOGPSI_FREE=1, ENERGY_TOT=2, ENERGY_FIXED=3, ENERGY_NEW=4, REWEIGHT=5};
    enum SumIndex_OPT {SUM_E_BARE, SUM_ESQ_BARE, SUM_ABSE_BARE,
      SUM_E_WGT, SUM_ESQ_WGT, SUM_ABSE_WGT, SUM_WGT, SUM_WGTSQ};


    ///Constructor.
    QMCCostFunction(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);
    
    ///Destructor
    ~QMCCostFunction();

    ///process xml node
    bool put(xmlNodePtr cur);
    ///assign optimization parameter i
    Return_t& Params(int i) { return OptParams[i]; }
    ///return optimization parameter i
    Return_t Params(int i) const { return OptParams[i]; }
    ///return the cost value for CGMinimization
    Return_t Cost();
    ///return the number of optimizable parameters
    int NumParams() { return OptParams.size(); }
    ///dump the current parameters and other report
    void Report();
    ///report  parameters at the end
    void reportParameters(std::ostream& os);
    ///evaluate the local energies of the samples
    Return_t correlatedSampling();

    void setWaveFunctionNode(xmlNodePtr cur) { m_wfPtr=cur; }

    void getConfigurations(vector<string>& ConfigFile, int partid, int nparts);
    void getConfigurations(const string& aroot);

    void checkConfigurations();

    void setTargetEnergy(Return_t et);

    void setRootName(const string& aroot) { RootName=aroot;}

    void setStream(ostream* os) { msg_stream = os;}

    /** set myComm 
     * @param c Communicate*
     */
    void setCommunicator(Communicate* c)
    {
      myComm=c;
    }

  protected:

    ///walker ensemble
    MCWalkerConfiguration& W;

    ///trial function
    TrialWaveFunction& Psi;

    ///Hamiltonian
    QMCHamiltonian& H;

    ///boolean to turn on/off the psi^2/psi^2_old for correlated sampling
    bool UseWeight;
    ///bollean to turn on/off the use of anonymous buffer
    bool needBuffering;
    ///bollean to turn on/off the use of anonymous buffer for the ratio
    bool hamiltonianNeedRatio;
    /** |E-E_T|^PowerE is used for the cost function
     *
     * default PowerE=1
     */
    int PowerE;
    ///number of times cost function evaluated
    int NumCostCalls;
    ///total number of samples to use in correlated sampling
    int NumSamples;
    ///total number of optimizable variables
    int NumOptimizables;
    ///counter for output
    int ReportCounter;
    ///weights for energy and variance in the cost function
    Return_t w_en, w_var, w_abs;
    ///value of the cost function
    Return_t CostValue;
    ///target energy
    Return_t Etarget;
    ///real target energy with the Correlation Factor
    Return_t EtargetEff;
    ///effective number of walkers 
    Return_t NumWalkersEff;
    ///fraction of the number of walkers below which the costfunction becomes invalid
    Return_t MinNumWalkers;
    ///maximum weight beyond which the weight is set to 1
    Return_t MaxWeight;

    ///current Average
    Return_t curAvg;
    ///current Variance
    Return_t curVar;
    ///current weighted average (correlated sampling) 
    Return_t curAvg_w;
    ///current weighted variance (correlated sampling)
    Return_t curVar_w;
    ///current variance of SUM_ABSE_WGT/SUM_WGT
    Return_t curVar_abs;

    /** Rescaling factor to correct the target energy Etarget=(1+CorrelationFactor)*Etarget
     *
     * default CorrelationFactor=0.0;
     */
    Return_t CorrelationFactor;
    ///storage for previous values of the cost function
    std::deque<Return_t> costList;
    ///storage for previous sets of parameters
    std::deque<vector<Return_t> > paramList;
    ///parameters to be optimized
    vector<Return_t> OptParams;
    ///ID tag for each optimizable parameter
    vector<string> IDtag;  
    ///list of optimizables
    OptimizableSetType OptVariables;
    ///communicator
    Communicate* myComm;
    ///stream to which progress is sent
    ostream* msg_stream;
    ///xml node to be dumped
    xmlNodePtr m_wfPtr;
    ///document node to be dumped
    xmlDocPtr m_doc_out;
    ///parameters to be updated
    std::map<string,xmlNodePtr> paramNodes;
    ///attributes to be updated
    std::map<string,pair<xmlNodePtr,string> > attribNodes;
    ///string for the file root
    string RootName;
    ///Hamiltonians that depend on the optimization: KE
    QMCHamiltonian H_KE;

    /** Sum of energies and weights for averages
     *
     * SumValues[k] where k is one of SumIndex_opt
     */
    TinyVector<Return_t,8> SumValue;
    /** Saved properties of all the walkers
     *
     * Records(iw,field_id) returns the field_id value of the iw-th walker
     * field_id is one of FieldIndex_opt
     */
    Matrix<Return_t> Records;
    /** Fixed  Laplacian , \f$\nabla^2\ln\Psi\f$, components */
    ParticleSet::ParticleLaplacian_t dL;

    bool resetWaveFunctions();
    bool checkParameters();
    bool putOptParams();  
    void updateXmlNodes();
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
