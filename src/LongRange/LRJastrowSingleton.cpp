//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file LRHandlerTemp.h
 * @brief Define a LRHandler with two template parameters
 */
#include "LongRange/LRJastrowSingleton.h"

namespace qmcplusplus
{
//initialization of the static data
LRJastrowSingleton::LRHandlerType* LRJastrowSingleton::JastrowHandler=0;

LRJastrowSingleton::LRHandlerType*
LRJastrowSingleton::getHandler(ParticleSet& ref, double kc)
{
  if(JastrowHandler ==0)
  {
    app_log() << "  LRJastrowSingleton::getHanlder " << std::endl;
    JastrowHandler=new LRHandlerType(ref,kc);
    JastrowHandler->initBreakup(ref);
    return JastrowHandler;
  }
  else
  {
    app_log() << "  Copy JastrowHandler. " << std::endl;
    return new LRHandlerType(*JastrowHandler,ref);
  }
}
}
