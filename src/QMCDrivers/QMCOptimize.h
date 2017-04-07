//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


/** @file QMCOptimize.h
 * @brief Definition of QMCDriver which performs VMC and optimization.
 */
#ifndef QMCPLUSPLUS_QMCOPTIMIZATION_VMCSINGLE_H
#define QMCPLUSPLUS_QMCOPTIMIZATION_VMCSINGLE_H

#include "QMCDrivers/QMCDriver.h"
#include "Optimize/OptimizeBase.h"

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

class QMCOptimize: public QMCDriver
{
public:

  ///Constructor.
  QMCOptimize(MCWalkerConfiguration& w, TrialWaveFunction& psi,
              QMCHamiltonian& h, HamiltonianPool& hpool, WaveFunctionPool& ppool);

  ///Destructor
  ~QMCOptimize();

  ///Run the Optimization algorithm.
  bool run();
  ///process xml node
  bool put(xmlNodePtr cur);
  ///add a configuration file to the list of files
  void addConfiguration(const std::string& a);

  void setWaveFunctionNode(xmlNodePtr cur)
  {
    wfNode=cur;
  }

private:

  ///index to denote the partition id
  int PartID;
  ///total number of partitions that will share a set of configuratons
  int NumParts;
  ///total number of Warmup Blocks
  int WarmupBlocks;
  ///total number of Warmup Blocks
  int NumOfVMCWalkers;
  ///yes/no applicable only first time
  std::string SkipSampleGeneration;
  ///need to know HamiltonianPool to use OMP
  HamiltonianPool& hamPool;
  ///target cost function to optimize
  //QMCCostFunction* optTarget;
  QMCCostFunctionBase* optTarget;
  ///solver
  MinimizerBase<RealType>* optSolver;
  ///vmc engine
  QMCDriver* vmcEngine;
  ///xml node to be dumped
  xmlNodePtr wfNode;
  ///xml node for optimizer
  xmlNodePtr optNode;
  ///method for optimization, default conjugate gradient
  std::string optmethod;
  ///list of files storing configurations
  std::vector<std::string> ConfigFile;
  ///Copy Constructor (disabled).
  QMCOptimize(const QMCOptimize& a): QMCDriver(a),hamPool(a.hamPool) { }
  ///Copy operator (disabled).
  QMCOptimize& operator=(const QMCOptimize&)
  {
    return *this;
  }

  void generateSamples();
};
}
#endif
