//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_COSTFUNCTIONOMP_H
#define QMCPLUSPLUS_COSTFUNCTIONOMP_H

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
class QMCCostFunctionOMP: public QMCCostFunctionBase, public CloneManager
{
public:

  ///Constructor.
  QMCCostFunctionOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                     QMCHamiltonian& h, HamiltonianPool& hpool);

  ///Destructor
  ~QMCCostFunctionOMP();

  void getConfigurations(const std::string& aroot);
  void checkConfigurations();
  void resetPsi(bool final_reset=false);
  void resetWalkers();   
  void GradCost(std::vector<Return_t>& PGradient, const std::vector<Return_t>& PM, Return_t FiniteDiff=0);
  Return_t fillOverlapHamiltonianMatrices(Matrix<Return_t>& H2, Matrix<Return_t>& Hamiltonian, Matrix<Return_t>& Variance, Matrix<Return_t>& Overlap);
  Return_t fillOverlapHamiltonianMatrices(Matrix<Return_t>& Left, Matrix<Return_t>& Right, Matrix<Return_t>& Overlap);

protected:
  std::vector<QMCHamiltonian*> H_KE_Node;
  std::vector<Matrix<Return_t>*> RecordsOnNode;

  /** Temp derivative properties and Hderivative properties of all the walkers
  */
  std::vector<Matrix<Return_t>* > DerivRecords;
  std::vector<Matrix<Return_t>* > HDerivRecords;
  Return_t CSWeight;

  ///vmc walkers to clean up
  std::vector<int> nVMCWalkers;
  Return_t correlatedSampling(bool needGrad=true);
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1804 $   $Date: 2007-02-24 14:49:09 -0600 (Sat, 24 Feb 2007) $
 * $Id: QMCCostFunctionOMP.h 1804 2007-02-24 20:49:09Z jnkim $
 ***************************************************************************/
