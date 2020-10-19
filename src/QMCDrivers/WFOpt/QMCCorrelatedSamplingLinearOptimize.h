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


/** @file QMCCorrelatedSamplingLinearOptimize.h
 * @brief Definition of QMCDriver which performs VMC and optimization.
 */
#ifndef QMCPLUSPLUS_QMCCSLINEAROPTIMIZATION_H
#define QMCPLUSPLUS_QMCCSLINEAROPTIMIZATION_H

#include "QMCDrivers/WFOpt/QMCLinearOptimize.h"
#include "QMCDrivers/VMC/VMCLinearOpt.h"
#include "Optimize/NRCOptimization.h"
#if defined(QMC_CUDA)
#include "QMCDrivers/VMC/VMC_CUDA.h"
#endif

namespace qmcplusplus
{
/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC.
 */

class QMCCorrelatedSamplingLinearOptimize : public QMCLinearOptimize, private NRCOptimization<QMCTraits::RealType>
{
public:
  ///Constructor.
  QMCCorrelatedSamplingLinearOptimize(MCWalkerConfiguration& w,
                                      TrialWaveFunction& psi,
                                      QMCHamiltonian& h,
                                      Communicate* comm);

  ///Destructor
  ~QMCCorrelatedSamplingLinearOptimize();

  ///Run the Optimization algorithm.
  bool run();
  ///process xml node
  bool put(xmlNodePtr cur);
  RealType Func(RealType dl);
  virtual QMCRunType getRunType() { return QMCRunType::CS_LINEAR_OPTIMIZE; }

private:
  inline bool ValidCostFunction(bool valid)
  {
    if (!valid)
      app_log() << " Cost Function is Invalid. If this frequently, try reducing the step size of the line minimization "
                   "or reduce the number of cycles. "
                << std::endl;
    return valid;
  }

#if defined(QMC_CUDA)
  VMCcuda* vmcCSEngine;
#else
  VMCLinearOpt* vmcCSEngine;
#endif

  int NumOfVMCWalkers;
  ///Number of iterations maximum before generating new configurations.
  int Max_iterations;
  int nstabilizers;
  RealType stabilizerScale, bigChange, exp0, exp1, stepsize, savedQuadstep;
  std::string MinMethod, GEVtype, StabilizerMethod, GEVSplit;
  RealType w_beta;
  /// number of previous steps to orthogonalize to.
  int eigCG;
  /// total number of cg steps per iterations
  int TotalCGSteps;
};
} // namespace qmcplusplus
#endif
