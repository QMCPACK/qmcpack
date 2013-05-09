//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Ken Esler
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: kpesler@gmail.com
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_COSTFUNCTION_CUDA_H
#define QMCPLUSPLUS_COSTFUNCTION_CUDA_H

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
class QMCCostFunctionCUDA: public QMCCostFunctionBase, public CloneManager
{
public:
  typedef MCWalkerConfiguration::Walker_t Walker_t;

  ///Constructor.
  QMCCostFunctionCUDA( MCWalkerConfiguration& w, TrialWaveFunction& psi,
                       QMCHamiltonian& h, HamiltonianPool& hpool);

  ///Destructor
  ~QMCCostFunctionCUDA();

  void getConfigurations(const string& aroot);
  void checkConfigurations();
  void resetWalkers();   
  void GradCost(vector<Return_t>& PGradient, const vector<Return_t>& PM, Return_t FiniteDiff=0);
  Return_t fillOverlapHamiltonianMatrices
  (Matrix<Return_t>& H2, Matrix<Return_t>& Hamiltonian, Matrix<Return_t>& Variance, Matrix<Return_t>& Overlap);
  Return_t fillOverlapHamiltonianMatrices(Matrix<Return_t>& Left, Matrix<Return_t>& Right);
  Return_t fillOverlapHamiltonianMatrices(Matrix<Return_t>& Left, Matrix<Return_t>& Right, Matrix<Return_t>& Overlap);

protected:
  Matrix<Return_t> Records;
  typedef TrialWaveFunction::ValueMatrix_t ValueMatrix_t;
  typedef TrialWaveFunction::GradMatrix_t  GradMatrix_t;
  /** Temp derivative properties and Hderivative properties of all the walkers
  */
  vector<vector<Return_t> >  TempDerivRecords;
  vector<vector<Return_t> >  TempHDerivRecords;
  ValueMatrix_t LogPsi_Derivs, LocE_Derivs;
  ValueMatrix_t d2logPsi_opt, d2logPsi_fixed;
  GradMatrix_t   dlogPsi_opt,  dlogPsi_fixed;

  vector<Matrix<Return_t>*> RecordsOnNode;
  vector<Matrix<Return_t>* > DerivRecords;
  vector<Matrix<Return_t>* > HDerivRecords;

  ///number of vmc walkers
  int nVMCWalkers;
  Return_t CSWeight;
  void resetPsi(bool final_reset=false);
  Return_t correlatedSampling(bool needDerivs);
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1804 $   $Date: 2007-02-24 14:49:09 -0600 (Sat, 24 Feb 2007) $
 * $Id: QMCCostFunctionCUDA.h 1804 2007-02-24 20:49:09Z jnkim $
 ***************************************************************************/
