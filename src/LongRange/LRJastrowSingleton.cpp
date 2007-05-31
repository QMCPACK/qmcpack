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
    LRJastrowSingleton::getHandler(ParticleSet& ref, double kc) {
      if(JastrowHandler ==0) {
        app_log() << "  LRJastrowSingleton::getHanlder " << endl;
        JastrowHandler=new LRHandlerType(ref,kc);
        JastrowHandler->initBreakup(ref);
        return JastrowHandler;
      }
      else
      {
        app_log() << "  Copy JastrowHandler. " << endl;
        return new LRHandlerType(*JastrowHandler,ref);
      }
    }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
