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
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_COSTFUNCTIONBASE_H
#define QMCPLUSPLUS_COSTFUNCTIONBASE_H

#include "Configuration.h"
#include "Optimize/OptimizeBase.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Message/MPIObjectBase.h"
#include <deque>
#include <set>

#ifdef HAVE_LMY_ENGINE
//#include "Eigen/Dense"
#include "formic/utils/matrix.h"
#include "formic/utils/lmyengine/engine.h"
#endif

namespace qmcplusplus
{

class MCWalkerConfiguration;

/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC.
 */

class QMCCostFunctionBase: public CostFunctionBase<QMCTraits::RealType>, public MPIObjectBase
{
public:

  enum FieldIndex_OPT {LOGPSI_FIXED=0, LOGPSI_FREE=1, ENERGY_TOT=2, ENERGY_FIXED=3, ENERGY_NEW=4, REWEIGHT=5};
  enum SumIndex_OPT {SUM_E_BARE=0, SUM_ESQ_BARE, SUM_ABSE_BARE,
                     SUM_E_WGT, SUM_ESQ_WGT, SUM_ABSE_WGT, SUM_WGT, SUM_WGTSQ,
                     SUM_INDEX_SIZE
                    };


  ///Constructor.
  QMCCostFunctionBase(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);

  ///Destructor
  virtual ~QMCCostFunctionBase();

  ///process xml node
  bool put(xmlNodePtr cur);
  void resetCostFunction(std::vector<xmlNodePtr>& cset);
  ///assign optimization parameter i
  Return_t& Params(int i)
  {
    return OptVariables[i];
  }
  ///return optimization parameter i
  Return_t Params(int i) const
  {
    return OptVariables[i];
  }
  int getType(int i)
  {
    return OptVariables.getType(i);
  }
  ///return the cost value for CGMinimization
  Return_t Cost(bool needGrad=true);

  ///return the cost value for CGMinimization
  Return_t computedCost();
  void printEstimates();
  ///return the gradient of cost value for CGMinimization
  virtual void GradCost(std::vector<Return_t>& PGradient, const std::vector<Return_t>& PM, Return_t FiniteDiff=0) {};
  ///return the number of optimizable parameters
  inline int NumParams()
  {
    return OptVariables.size();
  }
  ///return the number of optimizable parameters
  inline int getNumSamples()
  {
    return NumSamples;
  }
  inline void setNumSamples(int newNumSamples)
  {  
    NumSamples = newNumSamples;     
  }
  ///reset the wavefunction
  virtual void resetPsi(bool final_reset=false)=0;

  inline void getParameterTypes(std::vector<int>& types)
  {
    return OptVariablesForPsi.getParameterTypeList(types);
  }

  ///dump the current parameters and other report
  void Report();
  ///report  parameters at the end
  void reportParameters();

  ///return the counter which keeps track of optimization steps
  inline int getReportCounter() const
  {
    return ReportCounter;
  }

  void setWaveFunctionNode(xmlNodePtr cur)
  {
    m_wfPtr=cur;
  }

  void recordParametersToPsi(Return_t e, Return_t v)
  {
    Psi.coefficientHistory.addParams(OptVariables,e,v);
  }
  void getAvgParameters(int N)
  {
    OptVariables=Psi.coefficientHistory.getAvgCoefficients(N);
  }
  //void getConfigurations(std::vector<std::string>& ConfigFile, int partid, int nparts);

  void setTargetEnergy(Return_t et);

  void setRootName(const std::string& aroot)
  {
    RootName=aroot;
  }

  void setStream(std::ostream* os)
  {
    msg_stream = os;
  }

  void addCoefficients(xmlXPathContextPtr acontext, const char* cname);

  /** implement the virtual function
   * @param x0 current parameters
   * @param gr gradients or conjugate gradients
   * @param dl return the displacelement to minimize the cost function
   * @param val_proj projected cost
   *
   * If successful, any optimization object updates the parameters by x0 + dl*gr
   * and proceeds with a new step.
   */
  bool lineoptimization(const std::vector<Return_t>& x0, const std::vector<Return_t>& gr, Return_t val0,
                        Return_t& dl, Return_t& val_proj, Return_t& lambda_max);

  virtual Return_t fillOverlapHamiltonianMatrices(Matrix<Return_t>& Left, Matrix<Return_t>& Right)=0;

#ifdef HAVE_LMY_ENGINE
  Return_t LMYEngineCost(const bool needDeriv, cqmc::engine::LMYEngine * EngineObj);
#endif

  virtual void getConfigurations(const std::string& aroot)=0;

  virtual void checkConfigurations()=0;

#ifdef HAVE_LMY_ENGINE
  virtual void engine_checkConfigurations(cqmc::engine::LMYEngine * EngineObj)=0;
#endif

  void setRng(std::vector<RandomGenerator_t*>& r);

  inline bool getneedGrads() const
  {
    return needGrads;
  }

  inline void setneedGrads(bool tf)
  {
    needGrads=tf;
  }
  inline void setDMC()
  {
    vmc_or_dmc=1.0;
  }

protected:

  ///walker ensemble
  MCWalkerConfiguration& W;

  ///trial function
  TrialWaveFunction& Psi;

  ///Hamiltonian
  QMCHamiltonian& H;

  /////boolean to turn on/off the psi^2/psi^2_old for correlated sampling
  //bool UseWeight;
  ///bollean to turn on/off the use of anonymous buffer
  bool needBuffering;
  ///bollean to turn on/off the use of anonymous buffer for the ratio
  bool hamiltonianNeedRatio;
  ///if true, do not write the *.opt.#.xml
  bool Write2OneXml;
  ///if true, use analytic derivatives for the non-local potential component
  bool useNLPPDeriv;
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
  Return_t w_en, w_var, w_abs, w_w;
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
  ///threshold to remove configurations from sample with |Psi_old| < SmallWeight
  Return_t SmallWeight;
  Return_t w_beta;
  std::string GEVType;
  Return_t vmc_or_dmc;
  bool needGrads;
  ///whether we are targeting an excited state
  std::string targetExcitedStr;
  ///whether we are targeting an excited state
  bool targetExcited;
  ///the shift to use when targeting an excited state
  double omega_shift;
  /** Rescaling factor to correct the target energy Etarget=(1+CorrelationFactor)*Etarget
   *
   * default CorrelationFactor=0.0;
   */
  Return_t CorrelationFactor;
  ///list of optimizables
  opt_variables_type OptVariables;
  /** full list of optimizables
   *
   * The size of OptVariablesForPsi is equal to or larger than
   * that of OptVariables due to the dependent variables.
   * This is used for TrialWaveFunction::resetParameters and
   * is normally the same as OptVariables.
   */
  opt_variables_type OptVariablesForPsi;
  /** index mapping for <equal> constraints
   *
   * - equalVarMap[i][0] : index in OptVariablesForPsi
   * - equalVarMap[i][1] : index in OptVariables
   */
  std::vector<TinyVector<int,2> > equalVarMap;
  /** index mapping for <negate> constraints
   *
   * - negateVarMap[i][0] : index in OptVariablesForPsi
   * - negateVarMap[i][1] : index in OptVariables
   */
  ///index mapping for <negative> constraints
  std::vector<TinyVector<int,2> > negateVarMap;
  ///stream to which progress is sent
  std::ostream* msg_stream;
  ///xml node to be dumped
  xmlNodePtr m_wfPtr;
  ///document node to be dumped
  xmlDocPtr m_doc_out;
  ///parameters to be updated
  std::map<std::string,xmlNodePtr> paramNodes;
  ///coefficients to be updated
  std::map<std::string,xmlNodePtr> coeffNodes;
  ///attributes to be updated
  std::map<std::string,std::pair<xmlNodePtr,std::string> > attribNodes;
  ///string for the file root
  std::string RootName;
  ///Hamiltonians that depend on the optimization: KE
  QMCHamiltonian H_KE;

  ///Random number generators
  std::vector<RandomGenerator_t*> RngSaved,MoverRng;
  std::string includeNonlocalH;


  /** Sum of energies and weights for averages
   *
   * SumValues[k] where k is one of SumIndex_opt
   */
  std::vector<Return_t> SumValue;
  /** Saved properties of all the walkers
   *
   * Records(iw,field_id) returns the field_id value of the iw-th walker
   * field_id is one of FieldIndex_opt
   */
  Matrix<Return_t> Records;
  ///** Saved derivative properties and Hderivative properties of all the walkers
  //*/
  //vector<std::vector<vector<Return_t> >* > DerivRecords;
  //vector<std::vector<vector<Return_t> >* > HDerivRecords;

  typedef ParticleSet::ParticleGradient_t ParticleGradient_t;
  typedef ParticleSet::ParticleLaplacian_t ParticleLaplacian_t;
  ///** Fixed  Gradients , \f$\nabla\ln\Psi\f$, components */
  std::vector<ParticleGradient_t*> dLogPsi;
  ///** Fixed  Laplacian , \f$\nabla^2\ln\Psi\f$, components */
  std::vector<ParticleLaplacian_t*> d2LogPsi;
  ///stream for debug
  std::ostream* debug_stream;

  bool checkParameters();
  void updateXmlNodes();


  virtual Return_t correlatedSampling(bool needGrad=true)=0;

  #ifdef HAVE_LMY_ENGINE
  virtual Return_t LMYEngineCost_detail(cqmc::engine::LMYEngine * EngineObj)
  {
    APP_ABORT("NOT IMPLEMENTED");
    return 0;
  }
  #endif

};
}
#endif
