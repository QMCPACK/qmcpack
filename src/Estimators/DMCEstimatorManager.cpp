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
/** @file DMCEstimatorManager.cpp
 */
#include "Estimators/DMCEstimatorManager.h"

namespace qmcplusplus {

  DMCEstimatorManager::DMCEstimatorManager(QMCHamiltonian& h): 
    ScalarEstimatorManager(h) {
    add(&energyEstimator,"elocal");
  }

  DMCEstimatorManager::~DMCEstimatorManager() { }

  bool DMCEstimatorManager::put(xmlNodePtr cur) {
    return true;
  }

  void DMCEstimatorManager::accumulate(MCWalkerConfiguration& W) {
    energyEstimator.collect(W.begin(),W.end());
    BinSize++;
    MyData[WEIGHT_INDEX]+=energyEstimator.getWeight();
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
