//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file QMCChooseBestParameters.h
 * @brief Definition of QMCDriver which performs VMC and optimization.
 */
#ifndef QMCPLUSPLUS_QMCCHOSEBESTPARAMETERS_VMCSINGLE_H
#define QMCPLUSPLUS_QMCCHOSEBESTPARAMETERS_VMCSINGLE_H

#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
#include "QMCDrivers/QMCCSLinearOptimizeWFmanagerOMP.h"

namespace qmcplusplus
{

/** @ingroup QMCDrivers
 * @brief Implements wave-function optimization
 *
 * Optimization by correlated sampling method with configurations
 * generated from VMC.
 */
class HamiltonianPool;

class QMCChooseBestParameters: public QMCDriver, public CloneManager
{
public:

  ///Constructor.
  QMCChooseBestParameters(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h, HamiltonianPool& hpool, WaveFunctionPool& ppool);

  ///Destructor
  ~QMCChooseBestParameters();

  ///Run the algorithm.
  bool run();
  ///process xml node
  bool put(xmlNodePtr cur);
  void setWaveFunctionNode(xmlNodePtr cur)
  {
    wfNode=cur;
  };

private:
  QMCCSLinearOptimizeWFmanagerOMP optTarget;
  int naverage;
  xmlNodePtr wfNode;
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 757 $   $Date: 2005-10-31 10:10:28 -0600 (Mon, 31 Oct 2005) $
 * $Id: QMCChooseBestParameters.h 757 2005-10-31 16:10:28Z jnkim $
 ***************************************************************************/
