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

#include "QMCDrivers/QMCCostFunctionBase.h"
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
class QMCCostFunctionCUDA: public QMCCostFunctionBase, public CloneManager
{
public:
  typedef MCWalkerConfiguration::Walker_t Walker_t;

  ///Constructor.
  QMCCostFunctionCUDA( MCWalkerConfiguration& w, TrialWaveFunction& psi,
                       QMCHamiltonian& h, HamiltonianPool& hpool);

  ///Destructor
  ~QMCCostFunctionCUDA();

  void getConfigurations(const std::string& aroot);
  void checkConfigurations();
  void GradCost(std::vector<Return_t>& PGradient, const std::vector<Return_t>& PM, Return_t FiniteDiff=0);
  Return_t fillOverlapHamiltonianMatrices(Matrix<Return_t>& Left, Matrix<Return_t>& Right);

protected:
  Matrix<Return_t> Records;
  typedef TrialWaveFunction::RealMatrix_t  RealMatrix_t;
  typedef TrialWaveFunction::ValueMatrix_t ValueMatrix_t;
  typedef TrialWaveFunction::GradMatrix_t  GradMatrix_t;
  /** Temp derivative properties and Hderivative properties of all the walkers
  */
  std::vector<std::vector<Return_t> >  TempDerivRecords;
  std::vector<std::vector<Return_t> >  TempHDerivRecords;
  RealMatrix_t LogPsi_Derivs, LocE_Derivs;
  ValueMatrix_t d2logPsi_opt, d2logPsi_fixed;
  GradMatrix_t   dlogPsi_opt,  dlogPsi_fixed;

  std::vector<Matrix<Return_t>*> RecordsOnNode;
  std::vector<Matrix<Return_t>* > DerivRecords;
  std::vector<Matrix<Return_t>* > HDerivRecords;

  Return_t CSWeight;
  void resetPsi(bool final_reset=false);
  Return_t correlatedSampling(bool needDerivs);
};
}
#endif
