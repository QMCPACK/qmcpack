//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file DMCEstimatorManager.h
 * @brief Manager class of estimators for DMC
 */
#ifndef QMCPLUSPLUS_DMC_SCALAR_ESTIMATORMANAGER_H
#define QMCPLUSPLUS_DMC_SCALAR_ESTIMATORMANAGER_H

#include "Estimators/ScalarEstimatorManager.h"
#include "Estimators/DMCEnergyEstimator.h"

namespace qmcplusplus {

  class DMCEstimatorManager: public ScalarEstimatorManager {

  public:

    DMCEstimatorManager(QMCHamiltonian& h);
    ~DMCEstimatorManager();

    ///process xml tag associated with estimators
    bool put(xmlNodePtr cur);

    void accumulate(MCWalkerConfiguration& W);
  
    DMCEnergyEstimator<RealType> energyEstimator;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
