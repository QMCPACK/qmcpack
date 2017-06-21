//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_WALKER_CONTROL_FACTORY_H
#define QMCPLUSPLUS_WALKER_CONTROL_FACTORY_H

#include "QMCDrivers/WalkerControlBase.h"

namespace qmcplusplus
{

/** function to create WalkerControlBase or its derived class
 * @param current number of walkers
 */
WalkerControlBase* createWalkerController(int nwtot,
    Communicate* comm, xmlNodePtr cur, bool reconfig=false);
}
#endif

