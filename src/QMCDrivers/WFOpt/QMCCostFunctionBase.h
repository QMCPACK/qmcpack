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

#include <deque>
#include <set>
#include "Configuration.h"
#include "Optimize/OptimizeBase.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "Message/MPIObjectBase.h"

#ifdef HAVE_LMY_ENGINE
//#include "Eigen/Dense"
#include "formic/utils/matrix.h"
#include "formic/utils/lmyengine/engine.h"
#endif

namespace qmcplusplus
{
class MCWalkerConfiguration;
class DescentEngine;

/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC.
 */

class QMCCostFunctionBase : public CostFunctionBase<QMCTraits::RealType>, public MPIObjectBase
{
public:
  enum FieldIndex_OPT
  {
    LOGPSI_FIXED = 0,
    LOGPSI_FREE  = 1,
    ENERGY_TOT   = 2,
    ENERGY_FIXED = 3,
    ENERGY_NEW   = 4,
    REWEIGHT     = 5
  };
  enum SumIndex_OPT
  {
    SUM_E_BARE = 0,
    SUM_ESQ_BARE,
    SUM_ABSE_BARE,
    SUM_E_WGT,
    SUM_ESQ_WGT,
    SUM_ABSE_WGT,
    SUM_WGT,
    SUM_WGTSQ,
    SUM_INDEX_SIZE
  };


  ///Constructor.
  QMCCostFunctionBase(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, Communicate* comm);

  ///Destructor
  ~QMCCostFunctionBase() override;

  ///process xml node
  bool put(xmlNodePtr cur);
  void resetCostFunction(std::vector<xmlNodePtr>& cset);
  ///Save opt parameters to HDF5
  bool reportH5;
  bool CI_Opt;
  ///Path and name of the HDF5 prefix where CI coeffs are saved
  std::string newh5;
  ///assign optimization parameter i
  Return_t& Params(int i) override { return OptVariables[i]; }
  ///return optimization parameter i
  Return_t Params(int i) const override { return OptVariables[i]; }
  int getType(int i) const { return OptVariables.getType(i); }
  ///return the cost value for CGMinimization
  Return_rt Cost(bool needGrad = true) override;

  ///return the cost value for CGMinimization
  Return_rt computedCost();
  void printEstimates();
  ///return the gradient of cost value for CGMinimization
  void GradCost(std::vector<Return_rt>& PGradient,
                const std::vector<Return_rt>& PM,
                Return_rt FiniteDiff = 0) override{};
  ///return the number of optimizable parameters
  inline int getNumParams() const override { return OptVariables.size(); }
  ///return the global number of samples
  inline int getNumSamples() const { return NumSamples; }
  inline void setNumSamples(int newNumSamples) { NumSamples = newNumSamples; }
  ///reset the wavefunction
  virtual void resetPsi(bool final_reset = false) = 0;

  inline void getParameterTypes(std::vector<int>& types) { return OptVariablesForPsi.getParameterTypeList(types); }

  ///dump the current parameters and other report
  void Report() override;
  ///report  parameters at the end
  void reportParameters();

  ///report  parameters in HDF5 at the end
  void reportParametersH5();
  ///return the counter which keeps track of optimization steps
  inline int getReportCounter() const { return ReportCounter; }

  void setWaveFunctionNode(xmlNodePtr cur) { m_wfPtr = cur; }

  void setTargetEnergy(Return_rt et);

  void setRootName(const std::string& aroot) { RootName = aroot; }

  void setStream(std::ostream* os) { msg_stream = os; }

  void addCoefficients(xmlXPathContextPtr acontext, const char* cname);

  void printCJParams(xmlNodePtr cur, std::string& rname);

  void addCJParams(xmlXPathContextPtr acontext, const char* cname);

  virtual Return_rt fillOverlapHamiltonianMatrices(Matrix<Return_rt>& Left, Matrix<Return_rt>& Right) = 0;

#ifdef HAVE_LMY_ENGINE
  Return_rt LMYEngineCost(const bool needDeriv, cqmc::engine::LMYEngine<Return_t>* EngineObj);
#endif

  virtual void getConfigurations(const std::string& aroot) = 0;

  virtual void checkConfigurations() = 0;

#ifdef HAVE_LMY_ENGINE
  virtual void engine_checkConfigurations(cqmc::engine::LMYEngine<Return_t>* EngineObj,
                                          DescentEngine& descentEngineObj,
                                          const std::string& MinMethod) = 0;

#endif

  void setRng(RefVector<RandomGenerator> r);

  inline bool getneedGrads() const { return needGrads; }

  inline void setneedGrads(bool tf) { needGrads = tf; }
  inline void setDMC() { vmc_or_dmc = 1.0; }

  inline std::string getParamName(int i) const override { return OptVariables.name(i); }

  inline const opt_variables_type& getOptVariables() const { return OptVariables; }

  /// return variance after checkConfigurations
  inline Return_rt getVariance() const
  {
    return SumValue[SUM_ESQ_WGT] / SumValue[SUM_WGT] -
        (SumValue[SUM_E_WGT] / SumValue[SUM_WGT]) * (SumValue[SUM_E_WGT] / SumValue[SUM_WGT]);
  }

protected:
  ///walker ensemble
  MCWalkerConfiguration& W;

  ///trial function
  TrialWaveFunction& Psi;

  ///Hamiltonian
  QMCHamiltonian& H;

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
  /// global number of samples to use in correlated sampling
  int NumSamples;
  ///total number of optimizable variables
  int NumOptimizables;
  ///counter for output
  int ReportCounter;
  ///weights for energy and variance in the cost function
  Return_rt w_en, w_var, w_abs, w_w;
  ///value of the cost function
  Return_rt CostValue;
  ///target energy
  Return_rt Etarget;
  ///real target energy with the Correlation Factor
  Return_rt EtargetEff;
  ///effective number of walkers
  Return_rt NumWalkersEff;
  ///fraction of the number of walkers below which the costfunction becomes invalid
  Return_rt MinNumWalkers;
  ///maximum weight beyond which the weight is set to 1
  Return_rt MaxWeight;
  ///current Average
  Return_rt curAvg;
  ///current Variance
  Return_rt curVar;
  ///current weighted average (correlated sampling)
  Return_rt curAvg_w;
  ///current weighted variance (correlated sampling)
  Return_rt curVar_w;
  ///current variance of SUM_ABSE_WGT/SUM_WGT
  Return_rt curVar_abs;

  Return_rt w_beta;
  std::string GEVType;
  Return_rt vmc_or_dmc;
  bool needGrads;
  ///whether we are targeting an excited state
  std::string targetExcitedStr;
  ///whether we are targeting an excited state
  bool targetExcited;
  ///the shift to use when targeting an excited state
  double omega_shift;

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
  // unchanged initial checked-in variables
  opt_variables_type InitVariables;
  /** index mapping for <equal> constraints
   *
   * - equalVarMap[i][0] : index in OptVariablesForPsi
   * - equalVarMap[i][1] : index in OptVariables
   */
  std::vector<TinyVector<int, 2>> equalVarMap;
  /** index mapping for <negate> constraints
   *
   * - negateVarMap[i][0] : index in OptVariablesForPsi
   * - negateVarMap[i][1] : index in OptVariables
   */
  ///index mapping for <negative> constraints
  std::vector<TinyVector<int, 2>> negateVarMap;
  ///stream to which progress is sent
  std::ostream* msg_stream;
  ///xml node to be dumped
  xmlNodePtr m_wfPtr;
  ///document node to be dumped
  xmlDocPtr m_doc_out;
  ///parameters to be updated
  std::map<std::string, xmlNodePtr> paramNodes;
  ///coefficients to be updated
  std::map<std::string, xmlNodePtr> coeffNodes;
  ///attributes to be updated
  std::map<std::string, std::pair<xmlNodePtr, std::string>> attribNodes;
  ///string for the file root
  std::string RootName;

  ///Random number generators
  UPtrVector<RandomGenerator> RngSaved;
  std::vector<RandomGenerator*> MoverRng;
  std::string includeNonlocalH;


  /** Sum of energies and weights for averages
   *
   * SumValues[k] where k is one of SumIndex_opt
   */
  std::vector<Return_rt> SumValue;
  /** Saved properties of all the walkers
   *
   * Records(iw,field_id) returns the field_id value of the iw-th walker
   * field_id is one of FieldIndex_opt
   */
  Matrix<Return_rt> Records;
  ///** Saved derivative properties and Hderivative properties of all the walkers
  //*/
  //vector<std::vector<vector<Return_t> >* > DerivRecords;
  //vector<std::vector<vector<Return_t> >* > HDerivRecords;

  using ParticleGradient  = ParticleSet::ParticleGradient;
  using ParticleLaplacian = ParticleSet::ParticleLaplacian;
  ///** Fixed  Gradients , \f$\nabla\ln\Psi\f$, components */
  std::vector<ParticleGradient*> dLogPsi;
  ///** Fixed  Laplacian , \f$\nabla^2\ln\Psi\f$, components */
  std::vector<ParticleLaplacian*> d2LogPsi;
  ///stream for debug
  std::ostream* debug_stream;

  bool checkParameters();
  void updateXmlNodes();

  /// Flag on whether the variational parameter override is output to the new wavefunction
  bool do_override_output;

  virtual Return_rt correlatedSampling(bool needGrad = true) = 0;

#ifdef HAVE_LMY_ENGINE
  virtual Return_rt LMYEngineCost_detail(cqmc::engine::LMYEngine<Return_t>* EngineObj)
  {
    APP_ABORT("NOT IMPLEMENTED");
    return 0;
  }
#endif
};
} // namespace qmcplusplus
#endif
