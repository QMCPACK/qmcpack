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
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_COSTFUNCTION_H
#define QMCPLUSPLUS_COSTFUNCTION_H

#include "QMCDrivers/WFOpt/QMCCostFunctionBase.h"
#include "QMCDrivers/CloneManager.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"

namespace qmcplusplus
{
/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC running on a single thread.
 */
class QMCCostFunctionBatched : public QMCCostFunctionBase, public CloneManager
{
public:
  ///Constructor.
  QMCCostFunctionBatched(MCWalkerConfiguration& w,
                         TrialWaveFunction& psi,
                         QMCHamiltonian& h,
                         SampleStack& samples,
                         Communicate* comm);

  ///Destructor
  ~QMCCostFunctionBatched();

  void getConfigurations(const std::string& aroot);
  void checkConfigurations();
#ifdef HAVE_LMY_ENGINE
  void engine_checkConfigurations(cqmc::engine::LMYEngine<Return_t>* EngineObj,
                                  DescentEngine& descentEngineObj,
                                  const std::string& MinMethod);
#endif


  void resetPsi(bool final_reset = false);
  void GradCost(std::vector<Return_t>& PGradient, const std::vector<Return_t>& PM, Return_rt FiniteDiff = 0);
  Return_rt fillOverlapHamiltonianMatrices(Matrix<Return_rt>& Left, Matrix<Return_rt>& Right);

protected:
  std::vector<QMCHamiltonian*> H_KE_Node;
  std::vector<Matrix<Return_rt>*> RecordsOnNode;

  /** Temp derivative properties and Hderivative properties of all the walkers
  */
  std::vector<Matrix<Return_rt>*> DerivRecords;
  std::vector<Matrix<Return_rt>*> HDerivRecords;
  Return_rt CSWeight;

  Return_rt correlatedSampling(bool needGrad = true);

  SampleStack& samples_;

  std::vector<Return_rt> log_psi_fixed_;
  std::vector<Return_rt> log_psi_opt_;

  std::vector<TrialWaveFunction*> wf_ptr_list_;
  std::vector<ParticleSet*> p_ptr_list_;
  std::vector<QMCHamiltonian*> h_ptr_list_;
  std::vector<QMCHamiltonian*> h0_ptr_list_;


#ifdef HAVE_LMY_ENGINE
  int total_samples();
  Return_rt LMYEngineCost_detail(cqmc::engine::LMYEngine<Return_t>* EngineObj);
#endif
};
} // namespace qmcplusplus
#endif
