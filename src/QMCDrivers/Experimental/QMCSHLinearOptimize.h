//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file QMCSHLinearOptimize.h
 * @brief Definition of QMCDriver which performs VMC and optimization.
 */
#ifndef QMCPLUSPLUS_QMCSHOPTIMIZATION_VMCSINGLE_H
#define QMCPLUSPLUS_QMCSHOPTIMIZATION_VMCSINGLE_H

#include "QMCDrivers/QMCLinearOptimize.h"
#include "QMCDrivers/DMC/DMCOMPOPT.h"
#include "Optimize/NRCOptimization.h"

namespace qmcplusplus
{
/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC.
 */

class QMCSHLinearOptimize: public QMCLinearOptimize, private NRCOptimization<QMCTraits::RealType>
{
public:

  ///Constructor.
  QMCSHLinearOptimize(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                      QMCHamiltonian& h, HamiltonianPool& hpool, WaveFunctionPool& ppool);

  ///Destructor
  ~QMCSHLinearOptimize();

  ///Run the Optimization algorithm.
  bool run();
  ///process xml node
  bool put(xmlNodePtr cur);
  RealType Func(Return_t dl);

private:
  inline bool ValidCostFunction(bool valid)
  {
    if (!valid)
      app_log()<<" Cost Function is Invalid. If this frequently, try reducing the step size of the line minimization or reduce the number of cycles. " << std::endl;
    return valid;
  }
  DMCOMPOPT* dmcEngine;
  std::string MinMethod;
  RealType bigChange, w_beta, stepsize;

};
}
#endif
