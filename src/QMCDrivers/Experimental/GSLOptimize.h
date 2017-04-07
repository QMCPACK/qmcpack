//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_GSL_OPTIMIZE_H
#define QMCPLUSPLUS_GSL_OPTIMIZE_H

#include <deque>
#include "Configuration.h"
#include "QMCDrivers/QMCDriver.h"
#include "Optimize/Minimize.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus
{

/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 * @authors Jeongnim Kim, Jordan Vincent and Simone Chiesa
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC. Modify VMC_OPT to handle various cost functions.
 * This class is only for internal use but should not be distributed to the public domain.
 */
class GSLOptimize: public QMCDriver, public CostFunctionBase<OHMMS_PRECISION>
{
public:

  enum FieldIndex_OPT {LOGPSI_FIXED=0, LOGPSI_FREE=1, ENERGY_TOT=2, ENERGY_FIXED=3, ENERGY_NEW=4, REWEIGHT=5};
  enum SumIndex_OPT {SUM_E_BARE, SUM_ESQ_BARE, SUM_ABSE_BARE,
                     SUM_E_WGT, SUM_ESQ_WGT, SUM_ABSE_WGT, SUM_WGT, SUM_WGTSQ
                    };
  //enum {ENERGY, ENERGYSQ};

  ///Constructor.
  GSLOptimize(MCWalkerConfiguration& w, TrialWaveFunction& psi,
              QMCHamiltonian& h);

  ///Destructor
  ~GSLOptimize();

  ///Run the Optimization algorithm.
  bool run();
  ///process xml node
  bool put(xmlNodePtr cur);
  ///assign optimization parameter i
  Return_t& Params(int i)
  {
    return OptParams[i];
  }
  ///return optimization parameter i
  Return_t Params(int i) const
  {
    return OptParams[i];
  }
  ///return the cost value for CGMinimization
  Return_t Cost();

  ///evaluate the local energies of the samples
  RealType correlatedSampling();

  ///return the number of optimizable parameters
  int NumParams()
  {
    return OptParams.size();
  }
  ///add a configuration file to the list of files
  void addConfiguration(const std::string& a);

  void setWaveFunctionNode(xmlNodePtr cur)
  {
    m_wfPtr=cur;
  }

  void Report();

private:

  ///boolean to turn on/off the psi^2/psi^2_old for correlated sampling
  bool UseWeight;
  ///bollean to turn on/off the use of anonymous buffer
  bool needBuffering;
  ///bollean to turn on/off the use of anonymous buffer for the ratio
  bool hamiltonianNeedRatio;
  ///index to denote the partition id
  int PartID;
  ///total number of partitions that will share a set of configuratons
  int NumParts;
  /** |E-E_T|^PowerE is used for the cost function
   *
   * default PowerE=1
   */
  int PowerE;
  ///number of times cost function evaluated
  int NumCostCalls;
  ///total number of samples to use in correlated sampling
  int NumSamples;
  ///conjugate gradient tolerance
  RealType cg_tolerance;
  ///conjugate gradient stepsize
  RealType cg_stepsize;
  ///conjugate gradient epsilon
  RealType cg_epsilon;
  ///weights for energy and variance in the cost function
  RealType w_en, w_var, w_abs;
  ///value of the cost function
  RealType CostValue;
  ///target energy
  RealType Etarget;
  ///real target energy with the Correlation Factor
  RealType EtargetEff;
  ///minimum fraction of the effective walkers, default 0.9
  RealType MinNumWalkers;
  /** Rescaling factor to correct the target energy Etarget=(1+CorrelationFactor)*Etarget
   *
   * default CorrelationFactor=0.0;
   */
  RealType CorrelationFactor;
  ///xml node to be dumped
  xmlNodePtr m_wfPtr;
  ///document node to be dumped
  xmlDocPtr m_doc_out;
  ///parameters to be updated
  std::vector<xmlNodePtr> m_param_out;
  ///storage for previous values of the cost function
  std::deque<Return_t> costList;
  ///storage for previous sets of parameters
  std::deque<std::vector<Return_t> > paramList;
  ///parameters to be optimized
  std::vector<Return_t> OptParams;
  ///ID tag for each optimizable parameter
  std::vector<std::string> IDtag;
  ///method for optimization, default conjugate gradient
  std::string optmethod;
  ///list of files storing configurations
  std::vector<std::string> ConfigFile;
  ///Hamiltonians that depend on the optimization: KE
  QMCHamiltonian H_KE;

  /** Sum of energies and weights for averages
   *
   * SumValues[k] where k is one of SumIndex_opt
   */
  TinyVector<RealType,8> SumValue;
  /** Saved properties of all the walkers
   *
   * Records(iw,field_id) returns the field_id value of the iw-th walker
   * field_id is one of FieldIndex_opt
   */
  Matrix<ValueType> Records;
  /** Fixed  Laplacian , \f$\nabla^2\ln\Psi\f$, components */
  ParticleSet::ParticleLaplacian_t dL;

  void checkConfigurations();
  bool resetWaveFunctions();
  bool checkParameters();
  bool putOptParams();

  ///Copy Constructor (disabled).
  GSLOptimize(const GSLOptimize& a): QMCDriver(a) { }
  ///Copy operator (disabled).
  GSLOptimize& operator=(const GSLOptimize&)
  {
    return *this;
  }
};
}
#endif
