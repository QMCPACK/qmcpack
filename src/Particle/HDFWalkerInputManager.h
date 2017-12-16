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
