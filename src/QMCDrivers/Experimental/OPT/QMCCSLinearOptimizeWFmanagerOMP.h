//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_LINOPTWFMANOMP_H
#define QMCPLUSPLUS_LINOPTWFMANOMP_H

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
class QMCCSLinearOptimizeWFmanagerOMP: public QMCCostFunctionBase, public CloneManager
{
public:

  ///Constructor.
  QMCCSLinearOptimizeWFmanagerOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                                  QMCHamiltonian& h, HamiltonianPool& hpool);

  ///Destructor
  ~QMCCSLinearOptimizeWFmanagerOMP();

  void getConfigurations(const std::string& aroot) {};
  void checkConfigurations() {};
  void startOptimization();
  void resetPsi(bool final_reset=false);
  void GradCost(std::vector<Return_t>& PGradient, const std::vector<Return_t>& PM, Return_t FiniteDiff=0) {};
  Return_t fillOverlapHamiltonianMatrices(Matrix<Return_t>& H2, Matrix<Return_t>& Hamiltonian
                                          , Matrix<Return_t>& Variance, Matrix<Return_t>& Overlap)
  {
    return 0.0;
  }
protected:
  Return_t correlatedSampling(bool needGrad=true)
  {
    return 0.0;
  }
};
}
#endif
