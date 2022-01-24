//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_COSTFUNCTION_CUDA_H
#define QMCPLUSPLUS_COSTFUNCTION_CUDA_H

#include "QMCDrivers/WFOpt/QMCCostFunctionBase.h"
#include "QMCDrivers/CloneManager.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "type_traits/CUDATypes.h"
#include "HamiltonianRef.h"

namespace qmcplusplus
{
/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC running on a single thread.
 */
class QMCCostFunctionCUDA : public QMCCostFunctionBase, public CloneManager
{
public:
  using Walker_t = MCWalkerConfiguration::Walker_t;

  ///Constructor.
  QMCCostFunctionCUDA(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, Communicate* comm);

  ///Destructor
  ~QMCCostFunctionCUDA() override;

  void getConfigurations(const std::string& aroot) override;
  void checkConfigurations() override;
  void GradCost(std::vector<Return_rt>& PGradient, const std::vector<Return_rt>& PM, Return_rt FiniteDiff = 0) override;
  Return_rt fillOverlapHamiltonianMatrices(Matrix<Return_rt>& Left, Matrix<Return_rt>& Right) override;

protected:
  using CTS = CUDAGlobalTypes;
  Matrix<Return_rt> Records;
  using RealMatrix_t = TrialWaveFunction::RealMatrix_t;
  using ValueMatrix  = TrialWaveFunction::ValueMatrix;
  using GradMatrix   = TrialWaveFunction::GradMatrix;

  ///Hamiltonians that depend on the optimization: KE
  HamiltonianRef H_KE;
  /** Temp derivative properties and Hderivative properties of all the walkers
  */
  std::vector<std::vector<Return_rt>> TempDerivRecords;
  std::vector<std::vector<Return_rt>> TempHDerivRecords;
  RealMatrix_t LogPsi_Derivs, LocE_Derivs;
  ValueMatrix d2logPsi_opt, d2logPsi_fixed;
  GradMatrix dlogPsi_opt, dlogPsi_fixed;

  std::vector<Matrix<Return_rt>*> RecordsOnNode;
  std::vector<Matrix<Return_rt>*> DerivRecords;
  std::vector<Matrix<Return_rt>*> HDerivRecords;

  Return_rt CSWeight;
  void resetPsi(bool final_reset = false) override;
  Return_rt correlatedSampling(bool needDerivs) override;
};
} // namespace qmcplusplus
#endif
