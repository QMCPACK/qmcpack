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
/** @file QMCFixedSampleLinearOptimize.h
 * @brief Definition of QMCDriver which performs VMC and optimization.
 */
#ifndef QMCPLUSPLUS_QMCFSLINEAROPTIMIZATION_VMCSINGLE_H
#define QMCPLUSPLUS_QMCFSLINEAROPTIMIZATION_VMCSINGLE_H

#include "QMCDrivers/QMCLinearOptimize.h"
#include "Optimize/NRCOptimization.h"

namespace qmcplusplus
{
/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC.
 */

class QMCFixedSampleLinearOptimize: public QMCLinearOptimize, private NRCOptimization<QMCTraits::RealType>
{
public:

  ///Constructor.
  QMCFixedSampleLinearOptimize(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                               QMCHamiltonian& h, HamiltonianPool& hpool, WaveFunctionPool& ppool);

  ///Destructor
  ~QMCFixedSampleLinearOptimize();

  ///Run the Optimization algorithm.
  bool run();
  ///process xml node
  bool put(xmlNodePtr cur);
  RealType Func(Return_t dl);

private:
  inline bool ValidCostFunction(bool valid)
  {
    if (!valid)
      app_log()<<" Cost Function is Invalid. If this frequently, try reducing the step size of the line minimization or reduce the number of cycles. " <<endl;
    return valid;
  }

  int NumOfVMCWalkers;
  ///Number of iterations maximum before generating new configurations.
  int Max_iterations;
  int nstabilizers;
  RealType stabilizerScale, bigChange, exp0, exp1, stepsize, savedQuadstep;
  string MinMethod, GEVtype, StabilizerMethod, GEVSplit;
  RealType w_beta;
  /// number of previous steps to orthogonalize to.
  int eigCG;
  /// total number of cg steps per iterations
  int  TotalCGSteps;
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 757 $   $Date: 2005-10-31 10:10:28 -0600 (Mon, 31 Oct 2005) $
 * $Id: QMCFixedSampleLinearOptimize.h 757 2005-10-31 16:10:28Z jnkim $
 ***************************************************************************/
