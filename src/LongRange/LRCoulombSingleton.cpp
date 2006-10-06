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
/** @file LRHandlerTemp.h
 * @brief Define a LRHandler with two template parameters
 */
#include "LongRange/LRCoulombSingleton.h"

namespace qmcplusplus {
  //initialization of the static data
  LRCoulombSingleton::LRHandlerType* LRCoulombSingleton::CoulombHandler=0;

  LRCoulombSingleton::LRHandlerType*
    LRCoulombSingleton::getHandler(ParticleSet& ref) {
      if(CoulombHandler ==0) {
        CoulombHandler=new LRHandlerType(ref);
        CoulombHandler->initBreakup(ref);
      }
      return CoulombHandler;
    }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
