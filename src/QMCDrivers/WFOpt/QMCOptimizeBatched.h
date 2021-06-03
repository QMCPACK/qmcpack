//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file QMCOptimizeBatched.h
 * @brief Definition of QMCDriver which performs VMC (using batched driver) and optimization.
 */
#ifndef QMCPLUSPLUS_QMCOPTIMIZATION_VMCBATCHED_H
#define QMCPLUSPLUS_QMCOPTIMIZATION_VMCBATCHED_H

#include <memory>
#include "QMCDrivers/QMCDriver.h"
#include "Optimize/OptimizeBase.h"
#include "QMCDrivers/QMCDriverNew.h"
#include "QMCDrivers/QMCDriverInput.h"
#include "QMCDrivers/VMC/VMCDriverInput.h"
#include "QMCDrivers/VMC/VMCBatched.h"
#include "QMCDrivers/WFOpt/WFOptDriverInput.h"


namespace qmcplusplus
{
///forward declaration of a cost function
class QMCCostFunctionBase;
class HamiltonianPool;
class MCPopulation;

/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC.
 */

class QMCOptimizeBatched : public QMCDriverNew
{
public:
  ///Constructor.
  QMCOptimizeBatched(const ProjectData& project_data,
                     MCWalkerConfiguration& w,
                     QMCDriverInput&& qmcdriver_input,
                     VMCDriverInput&& vmcdriver_input,
                     WFOptDriverInput&& wfoptdriver_input,
                     MCPopulation&& population,
                     SampleStack& samples,
                     Communicate* comm);

  ///Destructor
  ~QMCOptimizeBatched();

  ///Run the Optimization algorithm.
  bool run();
  ///process xml node
  void process(xmlNodePtr cur);
  ///add a configuration file to the list of files
  void addConfiguration(const std::string& a);

  void setWaveFunctionNode(xmlNodePtr cur) { wfNode = cur; }
  QMCRunType getRunType() { return QMCRunType::OPTIMIZE_BATCH; }

private:
  ///index to denote the partition id
  int PartID;
  ///total number of partitions that will share a set of configuratons
  int NumParts;
  ///total number of VMC walkers
  int NumOfVMCWalkers;
  ///target cost function to optimize
  //QMCCostFunction* optTarget;
  std::unique_ptr<QMCCostFunctionBase> optTarget;
  ///solver
  MinimizerBase<RealType>* optSolver;
  ///vmc engine
  VMCBatched* vmcEngine;
  ///xml node to be dumped
  xmlNodePtr wfNode;
  ///xml node for optimizer
  xmlNodePtr optNode;
  ///method for optimization, default conjugate gradient
  std::string optmethod;
  ///list of files storing configurations
  std::vector<std::string> ConfigFile;
  ///Copy Constructor (disabled).
  QMCOptimizeBatched(const QMCOptimizeBatched&) = delete;
  ///Copy operator (disabled).
  QMCOptimizeBatched& operator=(const QMCOptimizeBatched&) = delete;

  void generateSamples();

  /// VMC-specific driver input
  VMCDriverInput vmcdriver_input_;

  WFOptDriverInput wfoptdriver_input_;

  /// Samples to use in optimizer
  SampleStack& samples_;

  // Need to keep this around, unfortunately, since QMCCostFunctionBatched uses QMCCostFunctionBase,
  // which still takes an MCWalkerConfiguration in the constructor.
  MCWalkerConfiguration& W;
};
} // namespace qmcplusplus
#endif
