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
/** @file QMCChooseBestParameters.h
 * @brief Definition of QMCDriver which performs VMC and optimization.
 */
#ifndef QMCPLUSPLUS_QMCCHOSEBESTPARAMETERS_VMCSINGLE_H
#define QMCPLUSPLUS_QMCCHOSEBESTPARAMETERS_VMCSINGLE_H

#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"

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

    private:
      ///total number of Warmup Blocks
      int WarmupBlocks;
      ///total number of Warmup Blocks
      int NumOfVMCWalkers;
      ///need to know HamiltonianPool to use OMP
      HamiltonianPool& hamPool;
      ///vmc engine
      QMCDriver* vmcEngine;
      ///pointer to trial wave function
      TrialWaveFunction* WF;
      ///cost for choosing best wf 0<=alpha<=1 0=energy, 1=variance
      RealType alpha;
      ///Copy Constructor (disabled).
      QMCChooseBestParameters(const QMCChooseBestParameters& a): QMCDriver(a),CloneManager(a.hamPool),hamPool(a.hamPool){ }
      ///Copy operator (disabled).
      QMCChooseBestParameters& operator=(const QMCChooseBestParameters&)
      {
        return *this;
      }
    };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 757 $   $Date: 2005-10-31 10:10:28 -0600 (Mon, 31 Oct 2005) $
 * $Id: QMCChooseBestParameters.h 757 2005-10-31 16:10:28Z jnkim $
 ***************************************************************************/
