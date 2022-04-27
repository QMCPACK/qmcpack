//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file QMCLinearOptimize.h
 * @brief Definition of QMCDriver which performs VMC and optimization.
 */
#ifndef QMCPLUSPLUS_QMCLINEAROPTIMIZATION_VMCSINGLE_H
#define QMCPLUSPLUS_QMCLINEAROPTIMIZATION_VMCSINGLE_H

#include <memory>
#include "QMCDrivers/QMCDriver.h"
#include "Optimize/OptimizeBase.h"
#include "QMCWaveFunctions/WaveFunctionPool.h"
#include "Numerics/LinearFit.h"
#ifdef HAVE_LMY_ENGINE
#include "formic/utils/matrix.h"
#include "formic/utils/lmyengine/engine.h"
#include "QMCDrivers/Optimizers/DescentEngine.h"
#endif
#include "LinearMethod.h"


namespace qmcplusplus
{
/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC.
 */

class QMCLinearOptimize : public QMCDriver, public LinearMethod
{
public:
  ///Constructor.
  QMCLinearOptimize(MCWalkerConfiguration& w,
                    TrialWaveFunction& psi,
                    QMCHamiltonian& h,
                    Communicate* comm,
                    const std::string& QMC_driver_type = "QMCLinearOptimize");

  ///Destructor
  ~QMCLinearOptimize() override = default;

  ///Run the Optimization algorithm.
  bool run() override = 0;
  ///process xml node
  bool put(xmlNodePtr cur) override;
  void setWaveFunctionNode(xmlNodePtr cur) { wfNode = cur; }

  std::vector<RealType> optdir, optparam;
  ///total number of VMC walkers
  int NumOfVMCWalkers;
  ///Number of iterations maximum before generating new configurations.
  int Max_iterations;
  ///target cost function to optimize
  std::unique_ptr<QMCCostFunctionBase> optTarget;
  ///vmc engine
  std::unique_ptr<QMCDriver> vmcEngine;
  ///xml node to be dumped
  xmlNodePtr wfNode;

  RealType param_tol;

  ///common operation to start optimization, used by the derived classes
  void start();
#ifdef HAVE_LMY_ENGINE
  void engine_start(cqmc::engine::LMYEngine<ValueType>* EngineObj,
                    DescentEngine& descentEngineObj,
                    std::string MinMethod);
#endif
  ///common operation to finish optimization, used by the derived classes
  void finish();
  void generateSamples();

  QMCRunType getRunType() override { return QMCRunType::LINEAR_OPTIMIZE; }
  Timer t1;
};
} // namespace qmcplusplus
#endif
