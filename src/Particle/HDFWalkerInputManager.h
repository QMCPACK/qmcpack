//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//                    Jeremy McMinnis, jmcminis@gmail.com, Navar Inc.
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Inc.
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_WALKER_INPUT_MANAGER_H
#define QMCPLUSPLUS_WALKER_INPUT_MANAGER_H

#include "Particle/MCWalkerConfiguration.h"
#include "OhmmsData/HDFAttribIO.h"
#include <queue>

class Communicate;

namespace qmcplusplus
{

class HDFWalkerInputManager
{

  MCWalkerConfiguration& targetW;
  Communicate* myComm;
  std::string CurrentFileRoot;

public:

  HDFWalkerInputManager(MCWalkerConfiguration& w, Communicate* c);
  ~HDFWalkerInputManager();
  bool put(xmlNodePtr cur);
  //bool put(std::vector<xmlNodePtr>& mset, int pid);
  //bool put(std::vector<xmlNodePtr>& mset, Communicate* comm);
  std::string getFileRoot()
  {
    return CurrentFileRoot;
  }

  void rewind(const std::string& h5root, int blocks);

};
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
