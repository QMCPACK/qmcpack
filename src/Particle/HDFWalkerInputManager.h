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

#include "Particle/WalkerConfigurations.h"
#include <queue>

class Communicate;

namespace qmcplusplus
{
class HDFWalkerInputManager
{
  /// reference to the list of walker configurations to be read from file
  WalkerConfigurations& wc_list_;
  /// number of particles
  const size_t num_ptcls_;
  Communicate* myComm;
  std::string CurrentFileRoot;

public:
  HDFWalkerInputManager(WalkerConfigurations& w, size_t num_ptcls, Communicate* c);
  ~HDFWalkerInputManager();
  bool put(xmlNodePtr cur);
  //bool put(std::vector<xmlNodePtr>& mset, int pid);
  //bool put(std::vector<xmlNodePtr>& mset, Communicate* comm);
  std::string getFileRoot() { return CurrentFileRoot; }
};
} // namespace qmcplusplus

#endif
