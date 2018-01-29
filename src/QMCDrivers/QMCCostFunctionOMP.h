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
    
    


#ifndef QMCPLUSPLUS_COSTFUNCTIONOMP_H
#define QMCPLUSPLUS_COSTFUNCTIONOMP_H

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
class QMCCostFunctionOMP: public QMCCostFunctionBase, public CloneManager
{
public:

  ///Constructor.
  QMCCostFunctionOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                     QMCHamiltonian& h, HamiltonianPool& hpool);

  ///Destructor
  ~QMCCostFunctionOMP();

  void getConfigurations(const std::string& aroot);
  void checkConfigurations();
#ifdef HAVE_LMY_ENGINE
  void engine_checkConfigurations(cqmc::engine::LMYEngine * EngineObj);
#endif
  void resetPsi(bool final_reset=false);
  void GradCost(std::vector<Return_t>& PGradient, const std::vector<Return_t>& PM, Return_t FiniteDiff=0);
  Return_t fillOverlapHamiltonianMatrices(Matrix<Return_t>& Left, Matrix<Return_t>& Right);

protected:
  std::vector<QMCHamiltonian*> H_KE_Node;
  std::vector<Matrix<Return_t>*> RecordsOnNode;

  /** Temp derivative properties and Hderivative properties of all the walkers
  */
  std::vector<Matrix<Return_t>* > DerivRecords;
  std::vector<Matrix<Return_t>* > HDerivRecords;
  Return_t CSWeight;

  Return_t correlatedSampling(bool needGrad=true);

  #ifdef HAVE_LMY_ENGINE
  int total_samples();
  Return_t LMYEngineCost_detail(cqmc::engine::LMYEngine * EngineObj);
  #endif

};
}
#endif
