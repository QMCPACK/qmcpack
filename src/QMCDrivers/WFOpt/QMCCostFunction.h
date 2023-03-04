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
#include "HamiltonianRef.h"

namespace qmcplusplus
{
/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC running on a single thread.
 */
class QMCCostFunction : public QMCCostFunctionBase, public CloneManager
{
public:
  ///Constructor.
  QMCCostFunction(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, Communicate* comm);

  ///Destructor
  ~QMCCostFunction() override;

  void getConfigurations(const std::string& aroot) override;
  void checkConfigurations(EngineHandle& handle) override;
#ifdef HAVE_LMY_ENGINE
  void engine_checkConfigurations(cqmc::engine::LMYEngine<Return_t>* EngineObj,
                                  DescentEngine& descentEngineObj,
                                  const std::string& MinMethod) override;
#endif


  void resetPsi(bool final_reset = false) override;
  void GradCost(std::vector<Return_rt>& PGradient, const std::vector<Return_rt>& PM, Return_rt FiniteDiff = 0) override;
  Return_rt fillOverlapHamiltonianMatrices(Matrix<Return_rt>& Left, Matrix<Return_rt>& Right) override;

protected:
  std::vector<std::unique_ptr<HamiltonianRef>> H_KE_Node;
  std::vector<Matrix<Return_rt>*> RecordsOnNode;

  /** Temp derivative properties and Hderivative properties of all the walkers
  */
  std::vector<Matrix<Return_rt>*> DerivRecords;
  std::vector<Matrix<Return_rt>*> HDerivRecords;
  Return_rt CSWeight;

  EffectiveWeight correlatedSampling(bool needGrad = true) override;

#ifdef HAVE_LMY_ENGINE
  size_t total_samples();
  Return_rt LMYEngineCost_detail(cqmc::engine::LMYEngine<Return_t>* EngineObj) override;
#endif

  NewTimer& fill_timer_;
};
} // namespace qmcplusplus
#endif
