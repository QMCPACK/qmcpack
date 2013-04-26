//////////////////////////////////////////////////////////////////
// (c) Copyright 2007- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
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
  void addConfiguration(const string& a);

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
  string SkipSampleGeneration;
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
  string optmethod;
  ///list of files storing configurations
  vector<string> ConfigFile;
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
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 757 $   $Date: 2005-10-31 10:10:28 -0600 (Mon, 31 Oct 2005) $
 * $Id: QMCOptimize.h 757 2005-10-31 16:10:28Z jnkim $
 ***************************************************************************/
