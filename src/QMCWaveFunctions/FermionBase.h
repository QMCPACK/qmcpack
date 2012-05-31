//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#ifndef QMCPLUSPLUS_FERMIONBASE_H
#define QMCPLUSPLUS_FERMIONBASE_H

#include "QMCWaveFunctions/SPOSetBase.h"

namespace qmcplusplus {

  /** base class for Fermion
   */
  struct FermionBase
  {
    map<string,SPOSetBasePtr> spo_map;
  };

}
#endif

/***************************************************************************
 * $RCSfile$   $Author: jmcminis $
 * $Revision: 5443 $   $Date: 2012-03-21 20:14:56 -0500 (Wed, 21 Mar 2012) $
 * $Id: FermionBase.h 5443 2012-03-22 01:14:56Z jmcminis $ 
 ***************************************************************************/

