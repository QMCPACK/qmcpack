//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


/** @file QMCLinearOptimizeBatched.h
 * @brief Definition of QMCDriver which performs VMC and optimization.
 */
#ifndef QMCPLUSPLUS_QMCLINEAROPTIMIZATION_BATCHED_H
#define QMCPLUSPLUS_QMCLINEAROPTIMIZATION_BATCHED_H

#include <memory>
#include "QMCDrivers/QMCDriverNew.h"
#include "QMCDrivers/QMCDriverInput.h"
#include "QMCDrivers/VMC/VMCDriverInput.h"
#include "QMCDrivers/VMC/VMCBatched.h"
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

class QMCLinearOptimizeBatched : public QMCDriverNew
{
public:
  ///Constructor.
  QMCLinearOptimizeBatched(const ProjectData& project_data,
                           MCWalkerConfiguration& w,
                           QMCDriverInput&& qmcdriver_input,
                           VMCDriverInput&& vmcdriver_input,
                           MCPopulation&& population,
                           SampleStack& samples,
                           Communicate* comm,
                           const std::string& QMC_driver_type = "QMCLinearOptimizeBatched");

  ///Destructor
  ~QMCLinearOptimizeBatched() override = default;

  ///Run the Optimization algorithm.
  bool run() override = 0;

  void setWaveFunctionNode(xmlNodePtr cur) { wfNode = cur; }

  // ------------------------------------
  // Used by legacy linear method algos

  std::vector<RealType> optdir, optparm;

  ///Number of iterations maximum before generating new configurations.
  int Max_iterations;

  RealType param_tol;
  //-------------------------------------

  ///target cost function to optimize
  std::unique_ptr<QMCCostFunctionBase> optTarget;
  ///vmc engine
  std::unique_ptr<VMCBatched> vmcEngine;
  ///xml node to be dumped
  xmlNodePtr wfNode;
  ///xml node for optimizer
  xmlNodePtr optNode;

  VMCDriverInput vmcdriver_input_;
  SampleStack& samples_;

  QMCRunType getRunType() override { return QMCRunType::LINEAR_OPTIMIZE; }

  ParameterSet m_param;

};
} // namespace qmcplusplus
#endif
