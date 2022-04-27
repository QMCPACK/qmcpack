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


namespace qmcplusplus
{
///forward declaration of a cost function
class QMCCostFunctionBase;
class HamiltonianPool;


/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC.
 */

class QMCLinearOptimize : public QMCDriver
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
  ///index to denote the partition id
  int PartID;
  ///total number of partitions that will share a set of configuratons
  int NumParts;
  ///total number of VMC walkers
  int NumOfVMCWalkers;
  ///Number of iterations maximum before generating new configurations.
  int Max_iterations;
  ///target cost function to optimize
  std::unique_ptr<QMCCostFunctionBase> optTarget;
  ///Dimension of matrix and number of parameters
  int N, numParams;
  ///vmc engine
  std::unique_ptr<QMCDriver> vmcEngine;
  ///xml node to be dumped
  xmlNodePtr wfNode;
  ///xml node for optimizer
  xmlNodePtr optNode;

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
  //asymmetric generalized EV
  RealType getLowestEigenvector(Matrix<RealType>& A, Matrix<RealType>& B, std::vector<RealType>& ev);
  //asymmetric EV
  RealType getLowestEigenvector(Matrix<RealType>& A, std::vector<RealType>& ev);
  void getNonLinearRange(int& first, int& last);
  bool nonLinearRescale(std::vector<RealType>& dP, Matrix<RealType>& S);
  RealType getNonLinearRescale(std::vector<RealType>& dP, Matrix<RealType>& S);
  void generateSamples();

  QMCRunType getRunType() override { return QMCRunType::LINEAR_OPTIMIZE; }
  NewTimer& generate_samples_timer_;
  NewTimer& initialize_timer_;
  NewTimer& eigenvalue_timer_;
  NewTimer& line_min_timer_;
  NewTimer& cost_function_timer_;
  Timer t1;
};
} // namespace qmcplusplus
#endif
