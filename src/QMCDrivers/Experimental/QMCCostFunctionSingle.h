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
    
    
#ifndef QMCPLUSPLUS_COSTFUNCTIONSINGLE_H
#define QMCPLUSPLUS_COSTFUNCTIONSINGLE_H

#include "QMCDrivers/QMCCostFunctionBase.h"

namespace qmcplusplus
{

/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC running on a single thread.
 */
class QMCCostFunctionSingle: public QMCCostFunctionBase
{
public:

  ///Constructor.
  QMCCostFunctionSingle(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h);

  ///Destructor
  ~QMCCostFunctionSingle();

  void getConfigurations(const std::string& aroot);
  void checkConfigurations();
  void GradCost(std::vector<QMCTraits::RealType>& PGradient, const std::vector<QMCTraits::RealType>& PM, QMCTraits::RealType FiniteDiff=0) ;
  Return_t fillOverlapHamiltonianMatrices(Matrix<Return_t>& H2, Matrix<Return_t>& Hamiltonian, Matrix<Return_t>& Variance, Matrix<Return_t>& Overlap);
protected:
  std::vector<std::vector<Return_t> > TempDerivRecords;
  std::vector<std::vector<Return_t> > TempHDerivRecords;
  Return_t CSWeight;
  void resetPsi(bool final_reset=false);
  Return_t correlatedSampling(bool needGrad=true);
};
}
#endif
