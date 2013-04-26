//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_LINOPTWFMANOMP_H
#define QMCPLUSPLUS_LINOPTWFMANOMP_H

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
class QMCCSLinearOptimizeWFmanagerOMP: public QMCCostFunctionBase, public CloneManager
{
public:

  ///Constructor.
  QMCCSLinearOptimizeWFmanagerOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                                  QMCHamiltonian& h, HamiltonianPool& hpool);

  ///Destructor
  ~QMCCSLinearOptimizeWFmanagerOMP();

  void getConfigurations(const string& aroot) {};
  void checkConfigurations() {};
  void startOptimization();
  void resetPsi(bool final_reset=false);
  void GradCost(vector<Return_t>& PGradient, const vector<Return_t>& PM, Return_t FiniteDiff=0) {};
  Return_t fillOverlapHamiltonianMatrices(Matrix<Return_t>& H2, Matrix<Return_t>& Hamiltonian
                                          , Matrix<Return_t>& Variance, Matrix<Return_t>& Overlap)
  {
    return 0.0;
  }
protected:
  Return_t correlatedSampling(bool needGrad=true)
  {
    return 0.0;
  }
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1804 $   $Date: 2007-02-24 14:49:09 -0600 (Sat, 24 Feb 2007) $
 * $Id: QMCCSLinearOptimizeWFmanagerOMP.h 1804 2007-02-24 20:49:09Z jnkim $
 ***************************************************************************/
