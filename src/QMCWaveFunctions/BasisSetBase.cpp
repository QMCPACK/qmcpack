//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "QMCWaveFunctions/BasisSetBase.h"
#include "Numerics/MatrixOperators.h"

namespace qmcplusplus {

  BasisSetBase::BasisSetBase():BasisSetSize(0) { }
  BasisSetBase::~BasisSetBase() { }

  void 
  BasisSetBase::resize() {
    if(BasisSetSize) {
      Phi.resize(BasisSetSize);
      dPhi.resize(BasisSetSize);
      d2Phi.resize(BasisSetSize);
      Temp.resize(BasisSetSize,MAXINDEX);
    } else {
      app_error() << "  BasisSetBase::BasisSetSize == 0" << endl;
    }
  }
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

