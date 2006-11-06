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
#include "LongRange/LRJastrowSingleton.h"

namespace qmcplusplus {
  //initialization of the static data
  LRJastrowSingleton::LRHandlerType* LRJastrowSingleton::JastrowHandler=0;

  LRJastrowSingleton::LRHandlerType*
    LRJastrowSingleton::getHandler(ParticleSet& ref) {
      if(JastrowHandler ==0) {
        JastrowHandler=new LRHandlerType(ref);
        JastrowHandler->initBreakup(ref);
      }
      return JastrowHandler;
    }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
